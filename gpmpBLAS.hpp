/**
 * interface file
 */

#ifndef GPMPBLAS_HPP
#define GPMPBLAS_HPP

#include <stdint.h>

extern "C" {
    void DGEMM_(char        *TRANSA, 
                char        *TRANSB, 
                uint32_t    M, 
                uint32_t    N, 
                uint32_t    K, 
                double      ALPHA, 
                double      *A, 
                uint32_t    LDA, 
                double      *B, 
                uint32_t    LDB, 
                double      BETA, 
                double      *C, 
                uint32_t    LDC);
}

void dgemm(char        *TRANSA, 
           char        *TRANSB, 
           uint32_t    M,  
           uint32_t    N,  
           uint32_t    K,  
           double      ALPHA, 
           double      *A, 
           uint32_t    LDA, 
           double      *B, 
           uint32_t    LDB, 
           double      BETA, 
           double      *C, 
           uint32_t    LDC){

    DGEMM_(TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC);

}

#endif
