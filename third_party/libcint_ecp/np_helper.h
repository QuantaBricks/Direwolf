/* Minimal shim for PySCF's np_helper/np_helper.h, trimmed to the symbols
 * nr_ecp.c/nr_ecp_deriv.c actually use (NPdset0, MIN, MAX). */
#ifndef ENGINE_LIBCINT_ECP_NP_HELPER_H
#define ENGINE_LIBCINT_ECP_NP_HELPER_H
#include <stddef.h>

#ifndef MIN
#define MIN(X, Y)       ((X) < (Y) ? (X) : (Y))
#endif
#ifndef MAX
#define MAX(X, Y)       ((X) > (Y) ? (X) : (Y))
#endif

static inline void NPdset0(double *p, const size_t n) {
    size_t i;
    for (i = 0; i < n; i++) {
        p[i] = 0;
    }
}

#endif
