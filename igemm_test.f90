module igemm
    implicit none

    integer, parameter :: MC = 384
    integer, parameter :: KC = 384
    integer, parameter :: NC = 4096
    integer, parameter :: MR = 4
    integer, parameter :: NR = 4

contains

    subroutine pack_MRxk(k, A, incRowA, incColA, buffer)
        integer, intent(in) :: k, incRowA, incColA
        double precision, intent(in) :: A(incRowA, *)
        double precision, intent(out) :: buffer(MR, KC)
        integer :: i, j

        do j = 1, k
            do i = 1, MR
                buffer(i, j) = A((i - 1)*incRowA + 1)
            end do
        end do
    end subroutine pack_MRxk

    subroutine pack_A(mc, kc, A, incRowA, incColA, buffer)
        integer, intent(in) :: mc, kc, incRowA, incColA
        double precision, intent(in) :: A(incRowA, *)
        double precision, intent(out) :: buffer(MC, KC)
        integer :: mp, _mr, i, j

        mp = mc/MR
        _mr = mc%MR

        do i = 1, mp
           call pack_MRxk(kc, A, incRowA, incColA, buffer((i - 1)*kc*MR + 1, 1))
            A = A + MR*incRowA
        end do

        if (_mr > 0) then
            do j = 1, kc
                do i = 1, _mr
                    buffer(i, j) = A((i - 1)*incRowA + 1)
                end do
                do i = _mr + 1, MR
                    buffer(i, j) = 0.0
                end do
                A = A + incColA
            end do
        end if
    end subroutine pack_A

    subroutine pack_kxNR(k, B, incRowB, incColB, buffer)
        integer, intent(in) :: k, incRowB, incColB
        double precision, intent(in) :: B(*, incColB)
        double precision, intent(out) :: buffer(KC, NR)
        integer :: i, j

        do i = 1, k
            do j = 1, NR
                buffer(i, j) = B(i, (j - 1)*incColB + 1)
            end do
        end do
    end subroutine pack_kxNR

    subroutine pack_B(kc, nc, B, incRowB, incColB, buffer)
        integer, intent(in) :: kc, nc, incRowB, incColB
        double precision, intent(in) :: B(*, incColB)
        double precision, intent(out) :: buffer(KC, NC)
        integer :: np, _nr, i, j

        np = nc/NR
        _nr = nc%NR

        do j = 1, np
           call pack_kxNR(kc, B, incRowB, incColB, buffer(1, (j - 1)*kc*NR + 1))
            B = B + NR*incColB
        end do

        if (_nr > 0) then
            do i = 1, kc
                do j = 1, _nr
                    buffer(i, j) = B(i, (j - 1)*incColB + 1)
                end do
                do j = _nr + 1, NR
                    buffer(i, j) = 0.0
                end do
                B = B + incRowB
            end do
        end if
    end subroutine pack_B

    subroutine dgemm_micro_kernel(kc, alpha, A, B, beta, C, incRowC, incColC)
        integer, intent(in) :: kc, incRowC, incColC
        double precision, intent(in) :: alpha, A(*), B(*)
        double precision, intent(inout) :: beta, C(incRowC, *)
        double precision :: AB(MR, NR)
        integer :: i, j, l

        AB = 0.0

        do l = 1, kc
            do j = 1, NR
                do i = 1, MR
                    AB(i, j) = AB(i, j) + A((i - 1)*MR + 1)*B((j - 1)*kc + i)
                end do
            end do
            A = A + MR
            B = B + NR
        end do

        if (beta == 0.0) then
            do j = 1, NR
                do i = 1, MR
                    C((i - 1)*incRowC + (j - 1)*incColC + 1) = 0.0
                end do
            end do
        elseif (beta /= 1.0) then
            do j = 1, NR
                do i = 1, MR
                    C((i - 1)*incRowC + (j - 1)*incColC + 1) = beta*C((i - 1)*incRowC + (j - 1)*incColC + 1)
                end do
            end do
        end if

        if (alpha == 1.0) then
            do j = 1, NR
                do i = 1, MR
                    C((i - 1)*incRowC + (j - 1)*incColC + 1) = C((i - 1)*incRowC + (j - 1)*incColC + 1) + AB(i, j)
                end do
            end do
        else
            do j = 1, NR
                do i = 1, MR
                    C((i - 1)*incRowC + (j - 1)*incColC + 1) = C((i - 1)*incRowC + (j - 1)*incColC + 1) + alpha*AB(i, j)
                end do
            end do
        end if
    end subroutine dgemm_micro_kernel

    subroutine dgemm_macro_kernel(mc, nc, kc, alpha, beta, C, incRowC, incColC)
        integer, intent(in) :: mc, nc, kc, incRowC, incColC
        double precision, intent(in) :: alpha, beta
        double precision, intent(inout) :: C(incRowC, *)
        double precision :: _A(MC, KC), _B(KC, NC), _C(MR, NR)
        integer :: mp, np, _mr, _nr, i, j, mr, nr

        mp = (mc + MR - 1)/MR
        np = (nc + NR - 1)/NR
        _mr = mod(mc, MR)
        _nr = mod(nc, NR)

        do j = 1, np
            nr = merge(NR, _nr, j /= np .or. _nr == 0)
            do i = 1, mp
                mr = merge(MR, _mr, i /= mp .or. _mr == 0)

                if (mr == MR .and. nr == NR) then
                  call dgemm_micro_kernel(kc, alpha, _A((i - 1)*kc*MR + 1, 1), &
                                            _B(1, (j - 1)*kc*NR + 1), beta, &
            C((i - 1)*MR*incRowC + (j - 1)*NR*incColC + 1, 1), incRowC, incColC)
                else
                  call dgemm_micro_kernel(kc, alpha, _A((i - 1)*kc*MR + 1, 1), &
                                            _B(1, (j - 1)*kc*NR + 1), 0.0, &
                                            _C, 1, MR)
                    call dgescal(mr, nr, beta, C((i - 1)*MR*incRowC + (j - 1)*NR*incColC + 1, 1), incRowC, incColC)
                    call dgeaxpy(mr, nr, 1.0, _C, 1, MR, &
            C((i - 1)*MR*incRowC + (j - 1)*NR*incColC + 1, 1), incRowC, incColC)
                end if
            end do
        end do
    end subroutine dgemm_macro_kernel

 subroutine dgemm_nn(m, n, k, alpha, A, incRowA, incColA, B, incRowB, incColB, &
                        beta, C, incRowC, incColC)
        integer, intent(in) :: m, n, k, incRowA, incColA, incRowB, incColB, incRowC, incColC
        double precision, intent(in) :: alpha, beta
        double precision, intent(in) :: A(incRowA, *), B(incRowB, *)
        double precision, intent(inout) :: C(incRowC, *)
        integer :: mb, nb, kb, _mc, _nc, _kc, mc, nc, kc, i, j, l

        double precision :: _A(MC, KC), _B(KC, NC), _C(MR, NR)
        double precision :: _beta

        mb = (m + MC - 1)/MC
        nb = (n + NC - 1)/NC
        kb = (k + KC - 1)/KC
        _mc = mod(m, MC)
        _nc = mod(n, NC)
        _kc = mod(k, KC)

        if (alpha == 0.0 .or. k == 0) then
            call dgescal(m, n, beta, C, incRowC, incColC)
            return
        end if

        do j = 1, nb
            nc = merge(NC, _nc, j /= nb .or. _nc == 0)

            do l = 1, kb
                kc = merge(KC, _kc, l /= kb .or. _kc == 0)
                _beta = merge(beta, 1.0, l == 1)

        call pack_B(kc, nc, B((l - 1)*KC*incRowB + 1, (j - 1)*NC*incColB + 1), &
                            incRowB, incColB, _B)

                do i = 1, mb
                    mc = merge(MC, _mc, i /= mb .or. _mc == 0)

        call pack_A(mc, kc, A((i - 1)*MC*incRowA + (l - 1)*KC*incColA + 1, 1), &
                                incRowA, incColA, _A)

                    call dgemm_macro_kernel(mc, nc, kc, alpha, _beta, &
                            C((i - 1)*MC*incRowC + (j - 1)*NC*incColC + 1, 1), &
                                            incRowC, incColC)
                end do
            end do
        end do
    end subroutine dgemm_nn

    subroutine dgeaxpy(m, n, alpha, X, incRowX, incColX, Y, incRowY, incColY)
        integer, intent(in) :: m, n, incRowX, incColX, incRowY, incColY
        double precision, intent(in) :: alpha
        double precision, intent(in) :: X(incRowX, *)
        double precision, intent(inout) :: Y(incRowY, *)
        integer :: i, j

        if (alpha /= 1.0) then
            do j = 1, n
                do i = 1, m
                    Y((i - 1)*incRowY + (j - 1)*incColY + 1) = Y((i - 1)*incRowY + (j - 1)*incColY + 1) + &
                                  alpha*X((i - 1)*incRowX + (j - 1)*incColX + 1)
                end do
            end do
        else
            do j = 1, n
                do i = 1, m
                    Y((i - 1)*incRowY + (j - 1)*incColY + 1) = Y((i - 1)*incRowY + (j - 1)*incColY + 1) + &
                                        X((i - 1)*incRowX + (j - 1)*incColX + 1)
                end do
            end do
        end if
    end subroutine dgeaxpy

    subroutine dgescal(m, n, alpha, X, incRowX, incColX)
        integer, intent(in) :: m, n, incRowX, incColX
        double precision, intent(in) :: alpha
        double precision, intent(inout) :: X(incRowX, *)
        integer :: i, j

        if (alpha /= 0.0) then
            do j = 1, n
                do i = 1, m
                    X((i - 1)*incRowX + (j - 1)*incColX + 1) = alpha*X((i - 1)*incRowX + (j - 1)*incColX + 1)
                end do
            end do
        else
            do j = 1, n
                do i = 1, m
                    X((i - 1)*incRowX + (j - 1)*incColX + 1) = 0.0
                end do
            end do
        end if
    end subroutine dgescal

end module igemm

