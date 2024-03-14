PROGRAM DGEMM_EX
    IMPLICIT NONE

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
            !PRINT '(F6.2, "  ",)', C(I, J)
            PRINT '(F6.2)', C(I, J)

        END DO
        PRINT *
    END DO

END PROGRAM DGEMM_EX

