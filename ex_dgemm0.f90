PROGRAM DGEMM_EX
    ! the DGEMM routine is ideal, opposed to the naive schoolbook On^3 
    ! implementation, for larger matrices

    EXTERNAL DGEMM

    ! Declarations
    INTEGER, PARAMETER :: M = 1024, N = 1024, K = 1024
    INTEGER :: LDA, LDB, LDC, I, J, KK
    DOUBLE PRECISION :: ALPHA, BETA
    DOUBLE PRECISION :: A(M, K), B(K, N), C(M, N)
    CHARACTER :: TRANSA, TRANSB
    ! MATMUL result mtx
    DOUBLE PRECISION :: MATMUL_MTX_RES(M,N), NAIVE_RES(M,N)
    LOGICAL :: EQUAL_MTX
    ! TIMERS
    REAL :: DGEMM_ST, DGEMM_END, MATMUL_ST, MATMUL_END, NAIVE_ST, NAIVE_END

    ! Initialize input matrices A and B
    CALL RANDOM_SEED()
    CALL RANDOM_NUMBER(A)
    CALL RANDOM_NUMBER(B)


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

    CALL CPU_TIME(MATMUL_ST)
    ! Call FORTRAN MATMUL routine using intrinsics for matrix multiplication
    MATMUL_MTX_RES = MATMUL(A, B)
    CALL CPU_TIME(MATMUL_END)

    ! Naive implementation of matrix multiplication
    NAIVE_RES = 0.0

    CALL CPU_TIME(NAIVE_ST)
    DO I = 1, M
        DO J = 1, N
            DO KK = 1, K
                NAIVE_RES(I, J) = NAIVE_RES(I, J) + A(I, KK) * B(KK, J)
            END DO
        END DO
    END DO

    CALL CPU_TIME(NAIVE_END)


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
            IF (ABS(C(I, J) - NAIVE_RES(I, J)) > 1.0E-10) THEN

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
    PRINT '("MATMUL time = ",f6.3," seconds")', MATMUL_END-MATMUL_ST
    PRINT '("NAIVE time = ",f6.3," seconds")', NAIVE_END-NAIVE_ST


END PROGRAM DGEMM_EX

