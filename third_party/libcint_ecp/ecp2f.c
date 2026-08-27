/*
 * Fortran-callable (trailing-underscore, pass-by-reference) wrappers
 * around PySCF's ECPscalar_cart/ECPscalar_sph, matching the calling
 * convention libcint's own autocode-generated cint1e_*_cart_/cint1e_*_sph_
 * wrappers use (see third_party/libcint's build output, e.g.
 * cint1e_ovlp_cart_ alongside cint1e_ovlp_cart): scalar arguments that
 * the real C function takes by value are taken by pointer here and
 * dereferenced, since Fortran always passes arguments by reference.
 *
 * The ECPOpt* and scratch-cache pointers are passed the same way
 * Engine already passes libcint's CINTOpt* from Fortran (see
 * mod_integrals_core.f90's "0_8" argument to cint1e_ovlp_cart): as a
 * plain integer(8) handle, 0 meaning NULL (let ECPscalar_cart allocate
 * its own scratch internally, same as passing a NULL CINTOpt* asks
 * libcint to run unoptimized).
 */
#include <stdint.h>
#include "cint.h"
#include "nr_ecp.h"

int ECPscalar_cart(double *out, int *dims, int *shls, int *atm, int natm,
                   int *bas, int nbas, double *env, ECPOpt *opt, double *cache);
int ECPscalar_sph(double *out, int *dims, int *shls, int *atm, int natm,
                  int *bas, int nbas, double *env, ECPOpt *opt, double *cache);

/*
 * Gradient kernels (nr_ecp_deriv.c), comp=3 (d/dx,d/dy,d/dz of the bra
 * shell): ECPscalar_ipnuc_cart sums the ECP potential over ALL
 * ECP-bearing atoms (env[AS_NECPBAS] worth of shells); ECPscalar_iprinv_cart
 * restricts to the single ECP atom named by env[AS_RINV_ORIG_ATOM] (set
 * by the caller before each call, same env slot the energy-path kernel
 * leaves at 0). See force.f90's ECP force assembly comment for the full
 * translational-invariance derivation these two feed into.
 */
int ECPscalar_ipnuc_cart(double *out, int *dims, int *shls, int *atm, int natm,
                         int *bas, int nbas, double *env, ECPOpt *opt, double *cache);
int ECPscalar_iprinv_cart(double *out, int *dims, int *shls, int *atm, int natm,
                          int *bas, int nbas, double *env, ECPOpt *opt, double *cache);
int ECPscalar_ipnuc_sph(double *out, int *dims, int *shls, int *atm, int natm,
                        int *bas, int nbas, double *env, ECPOpt *opt, double *cache);
int ECPscalar_iprinv_sph(double *out, int *dims, int *shls, int *atm, int natm,
                         int *bas, int nbas, double *env, ECPOpt *opt, double *cache);

int ecpscalar_cart_(double *out, FINT *dims, FINT *shls,
                     FINT *atm, FINT *natm, FINT *bas, FINT *nbas,
                     double *env, int64_t *opt, int64_t *cache)
{
        return ECPscalar_cart(out, (int *)dims, (int *)shls,
                               (int *)atm, (int)*natm, (int *)bas, (int)*nbas,
                               env, (ECPOpt *)(intptr_t)(*opt),
                               (double *)(intptr_t)(*cache));
}

int ecpscalar_sph_(double *out, FINT *dims, FINT *shls,
                    FINT *atm, FINT *natm, FINT *bas, FINT *nbas,
                    double *env, int64_t *opt, int64_t *cache)
{
        return ECPscalar_sph(out, (int *)dims, (int *)shls,
                              (int *)atm, (int)*natm, (int *)bas, (int)*nbas,
                              env, (ECPOpt *)(intptr_t)(*opt),
                              (double *)(intptr_t)(*cache));
}

int ecpscalar_ipnuc_cart_(double *out, FINT *dims, FINT *shls,
                          FINT *atm, FINT *natm, FINT *bas, FINT *nbas,
                          double *env, int64_t *opt, int64_t *cache)
{
        return ECPscalar_ipnuc_cart(out, (int *)dims, (int *)shls,
                                    (int *)atm, (int)*natm, (int *)bas, (int)*nbas,
                                    env, (ECPOpt *)(intptr_t)(*opt),
                                    (double *)(intptr_t)(*cache));
}

int ecpscalar_iprinv_cart_(double *out, FINT *dims, FINT *shls,
                           FINT *atm, FINT *natm, FINT *bas, FINT *nbas,
                           double *env, int64_t *opt, int64_t *cache)
{
        return ECPscalar_iprinv_cart(out, (int *)dims, (int *)shls,
                                     (int *)atm, (int)*natm, (int *)bas, (int)*nbas,
                                     env, (ECPOpt *)(intptr_t)(*opt),
                                     (double *)(intptr_t)(*cache));
}

int ecpscalar_ipnuc_sph_(double *out, FINT *dims, FINT *shls,
                         FINT *atm, FINT *natm, FINT *bas, FINT *nbas,
                         double *env, int64_t *opt, int64_t *cache)
{
        return ECPscalar_ipnuc_sph(out, (int *)dims, (int *)shls,
                                   (int *)atm, (int)*natm, (int *)bas, (int)*nbas,
                                   env, (ECPOpt *)(intptr_t)(*opt),
                                   (double *)(intptr_t)(*cache));
}

int ecpscalar_iprinv_sph_(double *out, FINT *dims, FINT *shls,
                          FINT *atm, FINT *natm, FINT *bas, FINT *nbas,
                          double *env, int64_t *opt, int64_t *cache)
{
        return ECPscalar_iprinv_sph(out, (int *)dims, (int *)shls,
                                    (int *)atm, (int)*natm, (int *)bas, (int)*nbas,
                                    env, (ECPOpt *)(intptr_t)(*opt),
                                    (double *)(intptr_t)(*cache));
}
