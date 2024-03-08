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
! vector_routines.f90

!> FORTRAN Subroutine for Matrix Addition on flattened matrices as arrays
!! of type float. Contains C++ wrapper function
!! @param A Addend A, an array representing a Matrix
!! @param B Addend B, an array representing a Matrix
!! @param C Sum C, an array representing the sum of A + B
!! @param mtx_size Assumes same size M x N
SUBROUTINE mtx_add_routine_float_(A, B, C, mtx_size) bind(C)
    USE :: ISO_FORTRAN_ENV
    USE :: ISO_C_BINDING

    INTEGER, INTENT(IN) :: mtx_size
    REAL(KIND=C_FLOAT), DIMENSION(mtx_size, mtx_size), INTENT(IN) :: A, B
    REAL(KIND=C_FLOAT), DIMENSION(mtx_size, mtx_size), INTENT(OUT) :: C

    C = A + B
END SUBROUTINE mtx_add_routine_float_

SUBROUTINE mtx_add_routine_int_(A, B, C, mtx_size) bind(C)
    USE :: ISO_FORTRAN_ENV
    USE :: ISO_C_BINDING

    INTEGER, INTENT(IN) :: mtx_size
    INTEGER(C_INT), DIMENSION(mtx_size, mtx_size), INTENT(IN) :: A, B
    INTEGER(C_INT), DIMENSION(mtx_size, mtx_size), INTENT(OUT) :: C

    C = A + B
END SUBROUTINE mtx_add_routine_int_

!> FORTRAN Subroutine for Matrix Multiplication using Fortran intrinsics.
!! Contains C++ wrapper function
!! @param a Multiplier a, an array representing a Matrix
!! @param b Multiplicand b, an array representing a Matrix
!! @param c Product c, an array representing the sum of a + b
!! @param nrows_a Number of rows
!! @param ncols Number of columns
SUBROUTINE mtx_mult(matrix1, matrix2, result, nrows1, ncols1, ncols2)
    implicit none
    INTEGER, INTENT(IN) :: nrows1, ncols1, ncols2
    REAL, INTENT(IN) :: matrix1(nrows1, ncols1), matrix2(ncols1, ncols2)
    REAL, INTENT(OUT) :: result(nrows1, ncols2)
    INTEGER :: i, j, k

    ! Perform matrix multiplication
    do i = 1, nrows1
        do j = 1, ncols2
            result(i, j) = 0.0
            do k = 1, ncols1
                result(i, j) = result(i, j) + matrix1(i, k)*matrix2(k, j)
            end do
        end do
    end do
END SUBROUTINE mtx_mult

