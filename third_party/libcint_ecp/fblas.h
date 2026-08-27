/* Minimal shim for PySCF's vhf/fblas.h, trimmed to the single symbol
 * nr_ecp.c/nr_ecp_deriv.c actually call (dgemm_). Engine already links a
 * BLAS providing this symbol at the standard trailing-underscore name. */
#ifndef ENGINE_LIBCINT_ECP_FBLAS_H
#define ENGINE_LIBCINT_ECP_FBLAS_H

#if defined __cplusplus
extern "C" {
#endif

void dgemm_(const char*, const char*,
            const int*, const int*, const int*,
            const double*, const double*, const int*,
            const double*, const int*,
            const double*, double*, const int*);

#if defined __cplusplus
} // end extern "C"
#endif

#endif
