PROGRAM DGEMM_EX
    EXTERNAL DGEMM

    ! Declarations
    INTEGER, PARAMETER :: M = 3, N = 3, K = 3
    INTEGER :: LDA, LDB, LDC, I, J
    DOUBLE PRECISION :: ALPHA, BETA
    DOUBLE PRECISION :: A(M, K), B(K, N), C(M, N)
    CHARACTER :: TRANSA, TRANSB

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

    ! Call the DGEMM subroutine for matrix multiplication
    CALL DGEMM(TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC)

    ! Print the result matrix C
    PRINT *, "Resultant Matrix C:"
    DO I = 1, M
        DO J = 1, N
            PRINT '(F6.2)', C(I, J)

        END DO
        PRINT *
    END DO

END PROGRAM DGEMM_EX

