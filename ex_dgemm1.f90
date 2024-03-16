PROGRAM DGEMM_EX
    EXTERNAL DGEMM

    ! Declarations
    INTEGER, PARAMETER :: M = 3, N = 3, K = 3
    INTEGER :: LDA, LDB, LDC, I, J
    DOUBLE PRECISION :: ALPHA, BETA
    DOUBLE PRECISION :: A(M, K), B(K, N), C(M, N)
    CHARACTER :: TRANSA, TRANSB
    ! MATMUL result mtx
    DOUBLE PRECISION :: MATMUL_MTX_RES(M,N)
    LOGICAL :: EQUAL_MTX
    ! TIMERS
    REAL :: MATMUL_ST, MAMTUL_END, DGEMM_ST, DGEMM_END

    ! Initialize input matrices A and B
    A = RESHAPE([1.0, 2.0, 3.0, &
                 4.0, 5.0, 6.0, &
                 7.0, 8.0, 9.0], [M, K])

    B = RESHAPE([1.0, 0.0, 0.0, &
                 0.0, 1.0, 0.0, &
                 0.0, 0.0, 1.0], [K, N])


    ! PARAMS
    ! 
    ! TRANSA    : N/T/C TRANSPOSE INPUT MATRIX A BEFORE MULT
    ! TRANSB    : N/T/C TRANSPOSE INPUT MATRIX B BEFORE MULT
    ! 
    ! M         : NUMBER OF ROWS FOR MATRIX A & MATRIX C
    ! N         : NUMBER OF COLS FOR MATRIX B & MATRIX C
    ! K         : NUMBER OF COLS FOR MATRIX A & ROWS FOR MATRIX B
    !
    ! ALPHA     : This scalar is multiplied by the product of the matrices A 
    !           and B before adding it to the matrix C. Essentially, it scales 
    !           the contribution of the product of A and B to the resulting 
    !           matrix C.
    !
    ! A         : INPUT MATRIX A
    ! LDA       : LEADING DIMENSION A strides between columns (rows)
    !
    ! B         : INPUT MATRIX B
    ! LDB       : LEADING DIMENSION B strides between rows (columns)
    ! BETA      : This scalar is multiplied by the existing values in matrix C 
    !           before adding the scaled product of matrices A and B. It 
    !           allows for the incorporation of the previous values of C into 
    !           the result
    !
    ! C         : RESULT MATRIX C
    ! LDC       : LEADING DIMENSION C strides between columns (rows)
    

    ! Set parameters
    TRANSA = 'N'
    TRANSB = 'N'
    ALPHA = 1.0
    BETA = 0.0
    LDA = M
    LDB = K
    LDC = M

    CALL CPU_TIME(DGEMM_ST)
    ! Call the DGEMM subroutine for matrix multiplication
    CALL DGEMM(TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC)
    CALL CPU_TIME(DGEMM_END)

    ! Call FORTRAN MATMUL routine for matrix multiplication
    MATMUL_MTX_RES = MATMUL(A, B)

    ! Print the result matrix C
    !PRINT *, "BLAS Matrix C:"
    !DO I = 1, M
    !    DO J = 1, N
    !        PRINT '(F6.2)', C(I, J)

    !    END DO
    !    PRINT *
    !END DO

    !PRINT *, "MATMUL Matrix C:"
    !DO I = 1, M
    !    DO J = 1, N
    !        PRINT '(F6.2)', MATMUL_MTX_RES(I, J)

    !    END DO
    !    PRINT *
    !END DO


    EQUAL_MTX = .TRUE.

    DO I = 1, M
        DO J = 1, N  
            IF (ABS(C(I, J) - MATMUL_MTX_RES(I, J)) > 1.0E-10) THEN

                EQUAL_MTX = .FALSE.

                EXIT
                
            END IF
        END DO
    END DO

    IF (EQUAL_MTX) THEN
        PRINT *, "Matrices C and MATMUL_MTX_RES are equal..."
    ELSE
        PRINT *, "Matrices C and MATMUL_MTX_RES not equal..."
    ENDIF

    PRINT '("DGEMM time = ",f6.3," seconds")', DGEMM_END-DGEMM_ST



END PROGRAM DGEMM_EX

