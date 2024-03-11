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
        DOUBLE PRECISION, INTENT(IN) :: A(incRowA,*)
        DOUBLE PRECISION, INTENT(OUT) :: buffer(MR, KC)
        INTEGER :: i, j

        DO j = 1, k
            DO i = 1, MR
                !buffer(i,j) = A((i-1)*incRowA + 1)
                buffer(i, j) = A(i, (j - 1) * incColA + 1)

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
            call pack_MRxk(kc, A((i-1)*MR+1, 1), incRowA, incColA, buffer((i-1)*kc*MR+1:i*kc*MR))

            END DO

        if (mr > 0) THEN
            do j = 1, kc
                buffer((mp*kc+j-1)*MR+1:(mp*kc+j)*MR) = A((mp*MR+1):(mp*MR+mr), j)
                buffer((mp*kc+j)*MR+1:(mp*kc+j+1)*MR) = 0.0
            end do
        end if



    END SUBROUTINE pack_A


END MODULE DGEMM 


