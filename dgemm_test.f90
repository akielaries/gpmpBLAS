!/*************************************************************************
! *
! *  Project
! *                         _____ _____  __  __ _____
! *                        / ____|  __ \|  \/  |  __ \
! *  ___  _ __   ___ _ __ | |  __| |__) | \  / | |__) |
! * / _ \| '_ \ / _ \ '_ \| | |_ |  ___/| |\/| |  ___/
! *| (_) | |_) |  __/ | | | |__| | |    | |  | | |
! * \___/| .__/ \___|_| |_|\_____|_|    |_|  |_|_|
! *      | |
! *      |_|
! *
! * Copyright (C) Akiel Aries, <akiel@akiel.org>, et al.
! *
! * This software is licensed as described in the file LICENSE, which
! * you should have received as part of this distribution. The terms
! * among other details are referenced in the official documentation
! * seen here : https://akielaries.github.io/openGPMP/ along with
! * important files seen in this project.
! *
! * You may opt to use, copy, modify, merge, publish, distribute
! * and/or sell copies of the Software, and permit persons to whom
! * the Software is furnished to do so, under the terms of the
! * LICENSE file. As this is an Open Source effort, all implementations
! * must be of the same methodology.
! *
! *
! *
! * This software is distributed on an AS IS basis, WITHOUT
! * WARRANTY OF ANY KIND, either express or implied.
! *
! ************************************************************************/
! dgemm.f90

MODULE DGEMM
    implicit none

    INTEGER, PARAMETER :: MC = 384
    INTEGER, PARAMETER :: KC = 384
    INTEGER, PARAMETER :: NC = 4096
    INTEGER, PARAMETER :: MR = 4
    INTEGER, PARAMETER :: NR = 4

contains

    ! packs panels from matrix A without padding
    SUBROUTINE pack_MRxk(k, A, incRowA, incColA, buffer)
        INTEGER, INTENT(IN) :: k, incRowA, incColA
        DOUBLE PRECISION, INTENT(IN) :: A(incRowA, *)
        DOUBLE PRECISION, INTENT(OUT) :: buffer(MR, KC)

        INTEGER :: i, j

        DO j = 1, k
            DO i = 1, MR
                !buffer(i,j) = A((i-1)*incRowA + 1)
                buffer(i, j) = A(i, (j - 1)*incColA + 1)

            END DO
        END DO
    END SUBROUTINE pack_MRxk

    ! packs panel from A with padding if needed
    SUBROUTINE pack_A(mc, kc, A, incRowA, incColA, buffer)
        INTEGER, INTENT(IN) :: mc, kc, incRowA, incColA
        DOUBLE PRECISION, INTENT(IN) :: A(mc, kc)
        DOUBLE PRECISION, INTENT(OUT) :: buffer(mc*kc)

        INTEGER :: mp, mr, i, j

        DO i = 1, mp
            CALL pack_MRxk(kc, A((i - 1)*MR + 1, 1), incRowA, incColA, buffer((i - 1)*kc*MR + 1:i*kc*MR))

        END DO

        IF (mr > 0) THEN
            DO j = 1, kc
                buffer((mp*kc + j - 1)*MR + 1:(mp*kc + j)*MR) = A((mp*MR + 1):(mp*MR + mr), j)
                buffer((mp*kc + j)*MR + 1:(mp*kc + j + 1)*MR) = 0.0
            END DO

        END IF

    END SUBROUTINE pack_A

    ! packs panels from B without padding
    SUBROUTINE pack_kxNR(k, B, incRowB, incColB, buffer)
        INTEGER, INTENT(IN) :: k, incRowB, incColB
        DOUBLE PRECISION, INTENT(IN) :: B(incRowB, *)
        DOUBLE PRECISION, INTENT(OUT) :: buffer(MR, KC)
        INTEGER :: i, j

        DO j = 1, k
            DO i = 1, MR

                buffer(i, j) = B(i, (j - 1)*incColB + 1)

            END DO
        END DO

    END SUBROUTINE pack_kxNR

    ! packs panel from B with padding if needed
    SUBROUTINE pack_B(kc, nc, B, incRowB, incColB, buffer)
        INTEGER, INTENT(IN) :: kc, nc, incRowB, incColB
        DOUBLE PRECISION, INTENT(IN) :: B(kc, nc)
        DOUBLE PRECISION, INTENT(OUT) :: buffer(kc*nc)

        INTEGER :: np, nr, i, j

        DO i = 1, np
            CALL pack_kxNR(kc, B((i - 1)*NR + 1, 1), incRowB, incColB, buffer((i - 1)*kc*NR + 1:i*kc*NR))

        END DO

        IF (nr > 0) THEN
            DO j = 1, kc
                buffer((np*kc + j - 1)*NR + 1:(np*nc + j)*NR) = B((np*NR + 1):(np*NR + nr), j)
                buffer((np*kc + j)*NR + 1:(np*kc + j + 1)*NR) = 0.0
            END DO

        END IF

    END SUBROUTINE pack_B

    ! micro kernel for multiplying panels from A and B
    SUBROUTINE micro_kernel(kc, alpha, A, B, beta, C, incRowC, incColC)
        INTEGER, INTENT(IN) :: kc, incRowC, incColC
        DOUBLE PRECISION, INTENT(IN) :: alpha, beta, A(MR), B(NR)
        DOUBLE PRECISION, INTENT(INOUT) :: C(MR, NR)

        DOUBLE PRECISION :: AB(MR*NR)
        INTEGER :: i, j, l

        ! Compute AB = A*B
        AB = 0.0
        DO l = 1, kc
            DO j = 1, NR
                DO i = 1, MR
                    AB(i, j) = AB(i, j) + A(i)*B(j)
                END DO
            END DO
            A = A + MR
            B = B + NR
        END DO

    END SUBROUTINE micro_kernel

END MODULE DGEMM

