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

module DGEMM 
    implicit none

    integer, parameter :: MC = 384
    integer, parameter :: KC = 384
    integer, parameter :: NC = 4096
    integer, parameter :: MR = 4
    integer, parameter :: NR = 4

contains

    subroutine pack_MRxk(k, A, incRowA, incColA, buffer)
        integer, intent(in) :: k, incRowA, incColA
        double precision, intent(in) :: A(incRowA,*)
        double precision, intent(out) :: buffer(MR, KC)
        integer :: i, j

        do j = 1, k
            do i = 1, MR
                !buffer(i,j) = A((i-1)*incRowA + 1)
                buffer(i, j) = A(i, (j - 1) * incColA + 1)

            end do
        end do
    end subroutine pack_MRxk



end module DGEMM 


