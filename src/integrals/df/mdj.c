/*
 * RI-J three-centre Coulomb work in the Hermite (McMurchie-Davidson) form,
 * the Split-RI-J arrangement: the orbital pair's density (b-pass) or the
 * fitted potential (J-pass) is contracted with the pair's Hermite expansion
 * once per primitive pair, so the per-aux-shell work is only the Hermite
 * Coulomb tensor R_{t+t'} and a dot product - no Cartesian pair components,
 * no cart->sph and no per-triple library call.
 *
 *   (mu nu|K) = sum_ab K_ab sum_t E^{mu nu}_t sum_g c_g lambda (-1)^{|t'|}
 *               E^K_{t'} R_{t+t'}(alpha, P - C)
 *   lambda = 2 pi^{5/2} / (p g sqrt(p+g)),  alpha = p g/(p+g)
 *
 * Basis data are read from the libcint atm/bas/env arrays (orbital shells
 * 0..norb-1, aux shells norb..nbas-1); conventions (common s/p factor,
 * cart->sph matrices) are libcint's, so results equal cint3c2e_{sph,cart}.
 */
#include <math.h>
#include <stdlib.h>
#include <string.h>

extern double *(*c2s_bra_sph[])(double *gsph, int nket, double *gcart, int l);
extern double CINTcommon_fac_sp(int l);
extern void gamma_inc_like(double *f, double t, int m);

#define LMAXO 4          /* orbital shells */
#define LMAXA 6          /* aux shells */
#define LMAXP (2*LMAXO)  /* pair Hermite order */
#define NMAX (LMAXP + LMAXA + 2)   /* +2: second derivatives */
#define DD (NMAX + 1)    /* cube edge of the R table */
#define NCART(l) (((l)+1)*((l)+2)/2)
#define NHERM(l) (((l)+1)*((l)+2)*((l)+3)/6)

/* Boys table: F_n(T) on T = k*BH, Taylor in 6 terms around the nearest node */
#define BH 0.05
#define BTMAX 40.0
#define BNT ((int)(BTMAX/BH) + 2)
#define BNN (NMAX + 7)
static double *btab = NULL;

static int norb, nbas_all, pure;
static const int *g_atm, *g_bas;
static const double *g_env;
/* c2s[l]: nsph x ncart row-major, times the libcint common factor */
static double *c2s[LMAXA + 1];
static int nfun[LMAXA + 1];
/* Hermite index lists: htu[l][k] = (t,u,v) of the k-th Hermite function with t+u+v<=l */
static int hoff[NMAX + 1][NHERM(NMAX)];   /* flat cube offset t*DD*DD+u*DD+v */
static int hsum[NMAX + 1][NHERM(NMAX)];   /* t+u+v */
/* aux primitives */
static int naux, *ax_first, *ax_nprim, *ax_l, *ax_foff, *ax_goff;
static double *ax_exp, *ax_cen, *ax_A;    /* ax_A: per primitive nfun x NHERM(l), sign folded in */
static long ghsize;

static void boys(double t, int m, double *f)
{
   if (t >= BTMAX) {
      double e = 0.0, b = 0.5/t;
      f[0] = 0.886226925452758013649/sqrt(t);
      for (int i = 1; i <= m; i++) f[i] = b*((2*i-1)*f[i-1] - e);
      return;
   }
   int k = (int)(t/BH + 0.5);
   double d = k*BH - t;
   const double *r = btab + (long)k*BNN + m;
   f[m] = r[0] + d*(r[1] + d*(r[2]/2 + d*(r[3]/6 + d*(r[4]/24 + d*(r[5]/120 + d*r[6]/720)))));
   if (m > 0) {
      double e = exp(-t);
      for (int i = m; i > 0; i--) f[i-1] = (2*t*f[i] + e)/(2*i-1);
   }
}

/* recursion tables for the Hermite index list (canonical order, prefix-closed in t+u+v) */
static int rdir[NHERM(NMAX)], rmul[NHERM(NMAX)];
static const int rstr[3] = { DD*DD, DD, 1 };

/* R_{tuv}(alpha, X) for t+u+v <= n into cube R (offsets t*DD*DD+u*DD+v) */
static void rtab(int n, double alpha, double X, double Y, double Z, double T, double *R, double *wk)
{
   double f[NMAX + 1];
   boys(T, n, f);
   if (n == 0) { R[0] = f[0]; return; }
   double m2a = -2.0*alpha, pw[NMAX + 1];
   pw[0] = 1.0;
   for (int i = 1; i <= n; i++) pw[i] = pw[i-1]*m2a;
   const double XYZ[3] = { X, Y, Z };
   /* level j holds R^{(n-j)} for t+u+v <= j; ping-pong between wk and R so the last level lands in R */
   double *cur = (n & 1) ? wk : R;
   cur[0] = f[n]*pw[n];
   for (int j = 1; j <= n; j++) {
      int lev = n - j;
      double *nx = (cur == R) ? wk : R;
      const double *pv = cur;
      nx[0] = f[lev]*pw[lev];
      int nk = NHERM(j);
      for (int k = 1; k < nk; k++) {
         int o = hoff[NMAX][k], st = rstr[rdir[k]], m = rmul[k];
         double val = XYZ[rdir[k]]*pv[o - st];
         if (m) val += m*pv[o - 2*st];
         nx[o] = val;
      }
      cur = nx;
   }
}

/* 1D Hermite coefficients E[i][j][t], i<=li, j<=lj */
#define EI(i,j,t) (((i)*(LMAXO+3) + (j))*(LMAXP+5) + (t))
#define EISZ ((LMAXO+3)*(LMAXO+3)*(LMAXP+5))
static void herm1d(int li, int lj, double p, double PA, double PB, double *E)
{
   double h = 0.5/p;
   E[EI(0,0,0)] = 1.0;
   for (int i = 0; i <= li; i++) {
      for (int j = 0; j <= lj; j++) {
         if (i == 0 && j == 0) continue;
         int pi, pj; double X;
         if (j == 0) { pi = i-1; pj = 0; X = PA; } else { pi = i; pj = j-1; X = PB; }
         int tm = pi + pj;
         for (int t = 0; t <= i + j; t++) {
            double v = 0.0;
            if (t > 0 && t-1 <= tm) v += h*E[EI(pi,pj,t-1)];
            if (t <= tm) v += X*E[EI(pi,pj,t)];
            if (t+1 <= tm) v += (t+1)*E[EI(pi,pj,t+1)];
            E[EI(i,j,t)] = v;
         }
      }
   }
}

static void cart_xyz(int l, int k, int *x, int *y, int *z)
{
   int c = 0;
   for (int lx = l; lx >= 0; lx--)
      for (int ly = l - lx; ly >= 0; ly--) {
         if (c == k) { *x = lx; *y = ly; *z = l - lx - ly; return; }
         c++;
      }
}

static long aoff_of_init(void);
static long v2_init(void);
long mdj_init(const int *atm, const int *bas, const int *nbas, const int *norb_sh,
              const double *env, const int *puream)
{
   g_atm = atm; g_bas = bas; g_env = env;
   nbas_all = *nbas; norb = *norb_sh; pure = *puream;
   if (!btab) {
      btab = malloc(sizeof(double)*(long)BNT*BNN);
      for (int k = 0; k < BNT; k++) gamma_inc_like(btab + (long)k*BNN, k*BH, BNN - 1);
   }
   for (int l = 0; l <= LMAXA; l++) {
      int nc = NCART(l), ns = pure ? 2*l+1 : nc;
      nfun[l] = ns;
      free(c2s[l]);
      c2s[l] = calloc((size_t)ns*nc, sizeof(double));
      double fac = CINTcommon_fac_sp(l);
      if (pure) {
         double *id = calloc((size_t)nc*nc, sizeof(double)), *out = calloc((size_t)nc*ns, sizeof(double));
         for (int k = 0; k < nc; k++) id[k*nc + k] = 1.0;
         double *r = c2s_bra_sph[l](out, nc, id, l);
         for (int k = 0; k < nc; k++)
            for (int s = 0; s < ns; s++) c2s[l][s*nc + k] = fac*r[k*ns + s];
         free(id); free(out);
      } else {
         for (int k = 0; k < nc; k++) c2s[l][k*nc + k] = fac;
      }
   }
   for (int n = 0; n <= NMAX; n++) {
      int k = 0;
      for (int s = 0; s <= n; s++)
         for (int t = s; t >= 0; t--)
            for (int u = s - t; u >= 0; u--) {
               int v = s - t - u;
               hoff[n][k] = t*DD*DD + u*DD + v; hsum[n][k] = s; k++;
            }
   }
   for (int k = 1; k < NHERM(NMAX); k++) {
      int o = hoff[NMAX][k], t = o/(DD*DD), u = (o/DD)%DD, v = o%DD;
      if (t > 0) { rdir[k] = 0; rmul[k] = t - 1; }
      else if (u > 0) { rdir[k] = 1; rmul[k] = u - 1; }
      else { rdir[k] = 2; rmul[k] = v - 1; }
   }
   /* aux primitives */
   naux = nbas_all - norb;
   free(ax_first); free(ax_nprim); free(ax_l); free(ax_foff); free(ax_goff);
   ax_first = malloc(sizeof(int)*naux); ax_nprim = malloc(sizeof(int)*naux);
   ax_l = malloc(sizeof(int)*naux); ax_foff = malloc(sizeof(int)*naux); ax_goff = malloc(sizeof(int)*naux);
   int np = 0, fo = 0; long go = 0;
   for (int k = 0; k < naux; k++) {
      const int *b = bas + 8*(norb + k);
      if (b[3] != 1 || b[1] > LMAXA) return -1;
      ax_first[k] = np; ax_nprim[k] = b[2]; ax_l[k] = b[1];
      ax_foff[k] = fo; ax_goff[k] = (int)go;
      np += b[2]; fo += nfun[b[1]]; go += (long)b[2]*NHERM(b[1]);
   }
   for (int i = 0; i < norb; i++) {
      const int *b = bas + 8*i;
      if (b[3] != 1 || b[1] > LMAXO) return -1;
   }
   ghsize = go;
  free(ax_exp); free(ax_cen); free(ax_A);
   ax_exp = malloc(sizeof(double)*np); ax_cen = malloc(sizeof(double)*3*np);
   long asz = 0;
   for (int k = 0; k < naux; k++) asz += (long)ax_nprim[k]*nfun[ax_l[k]]*NHERM(ax_l[k]);
   ax_A = calloc(asz, sizeof(double));
   long ao = 0;
   for (int k = 0; k < naux; k++) {
      const int *b = bas + 8*(norb + k);
      int l = b[1], nc = NCART(l), nh = NHERM(l), ns = nfun[l];
      const double *C = env + atm[6*b[0] + 1];
      for (int g = 0; g < b[2]; g++) {
         int q = ax_first[k] + g;
         double gam = env[b[5] + g], cg = env[b[6] + g];
         ax_exp[q] = gam; ax_cen[3*q] = C[0]; ax_cen[3*q+1] = C[1]; ax_cen[3*q+2] = C[2];
         /* one-centre 1D coefficients e[i][t] */
         double e[LMAXA+1][LMAXA+1];
         memset(e, 0, sizeof(e));
         e[0][0] = 1.0;
         for (int i = 1; i <= l; i++)
            for (int t = 0; t <= i; t++)
               e[i][t] = (t > 0 ? 0.5/gam*e[i-1][t-1] : 0.0) + (t+1 <= i-1 ? (t+1)*e[i-1][t+1] : 0.0);
         double *A = ax_A + ao;   /* ns x nh, row-major */
         for (int c = 0; c < nc; c++) {
            int x, y, z; cart_xyz(l, c, &x, &y, &z);
            for (int h = 0; h < nh; h++) {
               int o = hoff[l][h], t = o/(DD*DD), u = (o/DD)%DD, v = o%DD;
               if (t > x || u > y || v > z) continue;
               double val = cg*e[x][t]*e[y][u]*e[z][v]*(((t+u+v) & 1) ? -1.0 : 1.0);
               if (val == 0.0) continue;
               for (int s = 0; s < ns; s++) A[s*nh + h] += c2s[l][s*nc + c]*val;
            }
         }
         ao += (long)ns*nh;
      }
   }
   (void)aoff_of_init();
   return v2_init();
}

static long aoff_of(int k, int g);
static long aoff_of_init(void) { return naux > 0 ? aoff_of(0, 0) : 0; }
static long aoff_of(int k, int g)
{
   /* offset of primitive g of aux shell k into ax_A */
   static long *tab = NULL; static int tn = -1;
   if (tn != naux) {
      free(tab); tab = malloc(sizeof(long)*(naux + 1));
      long o = 0;
      for (int i = 0; i < naux; i++) { tab[i] = o; o += (long)ax_nprim[i]*nfun[ax_l[i]]*NHERM(ax_l[i]); }
      tab[naux] = o; tn = naux;
   }
   return tab[k] + (long)g*nfun[ax_l[k]]*NHERM(ax_l[k]);
}

/* shell pair primitive data */
typedef struct {
   int li, lj, ni, nj, np;
   double A[3], B[3];
   const double *ai, *aj, *ci, *cj;
} pair_t;

static void pair_setup(int ish, int jsh, pair_t *P)
{
   const int *bi = g_bas + 8*ish, *bj = g_bas + 8*jsh;
   P->li = bi[1]; P->lj = bj[1];
   const double *A = g_env + g_atm[6*bi[0] + 1], *B = g_env + g_atm[6*bj[0] + 1];
   memcpy(P->A, A, 3*sizeof(double)); memcpy(P->B, B, 3*sizeof(double));
   P->ni = bi[2]; P->nj = bj[2];
   P->ai = g_env + bi[5]; P->aj = g_env + bj[5];
   P->ci = g_env + bi[6]; P->cj = g_env + bj[6];
}

/* Cartesian block from a function block: Xc(ci,cj) = c2s_i^T X c2s_j; X is di x dj column-major */
static void to_cart(int li, int lj, const double *X, int ldx, double *Xc)
{
   int nci = NCART(li), ncj = NCART(lj), di = nfun[li], dj = nfun[lj];
   double tmp[NCART(LMAXO)*(2*LMAXO+1)];  /* nci x dj */
   for (int b = 0; b < dj; b++)
      for (int a = 0; a < nci; a++) {
         double s = 0.0;
         for (int m = 0; m < di; m++) s += c2s[li][m*nci + a]*X[b*ldx + m];
         tmp[b*nci + a] = s;
      }
   for (int bb = 0; bb < ncj; bb++)
      for (int a = 0; a < nci; a++) {
         double s = 0.0;
         for (int b = 0; b < dj; b++) s += tmp[b*nci + a]*c2s[lj][b*ncj + bb];
         Xc[bb*nci + a] = s;
      }
}

#define SCREEN_PRIM 1.0e-15

/* b-pass for one shell pair: G[aux hermite slots] += contributions of the density block Dsph (di x dj, ld di) */
void mdj_bpass_pair_ref(const int *ish_, const int *jsh_, const double *D, const double *schw,
                    const double *auxb, const double *cut, double *G)
{
   pair_t P; pair_setup(*ish_, *jsh_, &P);
   int li = P.li, lj = P.lj, L = li + lj, nci = NCART(li), ncj = NCART(lj), nhL = NHERM(L);
   double Dc[NCART(LMAXO)*NCART(LMAXO)], Ex[EISZ],
          Ey[EISZ], Ez[EISZ], Dh[NHERM(LMAXP)];
   double R[DD*DD*DD], wk[DD*DD*DD];
   to_cart(li, lj, D, nfun[li], Dc);
   double dmax = 0.0;
   for (int k = 0; k < nci*ncj; k++) dmax = fmax(dmax, fabs(Dc[k]));
   if (dmax == 0.0) return;
   double AB2 = (P.A[0]-P.B[0])*(P.A[0]-P.B[0]) + (P.A[1]-P.B[1])*(P.A[1]-P.B[1]) + (P.A[2]-P.B[2])*(P.A[2]-P.B[2]);
   int ix[NCART(LMAXO)], iy[NCART(LMAXO)], iz[NCART(LMAXO)], jx[NCART(LMAXO)], jy[NCART(LMAXO)], jz[NCART(LMAXO)];
   for (int a = 0; a < nci; a++) cart_xyz(li, a, ix+a, iy+a, iz+a);
   for (int b = 0; b < ncj; b++) cart_xyz(lj, b, jx+b, jy+b, jz+b);
   for (int pa = 0; pa < P.ni; pa++)
   for (int pb = 0; pb < P.nj; pb++) {
      double a = P.ai[pa], b = P.aj[pb], p = a + b;
      double Kab = exp(-a*b/p*AB2)*P.ci[pa]*P.cj[pb];
      if (fabs(Kab)*dmax*pow(M_PI/p, 1.5) < SCREEN_PRIM) continue;
      double Pc[3];
      for (int d = 0; d < 3; d++) Pc[d] = (a*P.A[d] + b*P.B[d])/p;
      herm1d(li, lj, p, Pc[0]-P.A[0], Pc[0]-P.B[0], Ex);
      herm1d(li, lj, p, Pc[1]-P.A[1], Pc[1]-P.B[1], Ey);
      herm1d(li, lj, p, Pc[2]-P.A[2], Pc[2]-P.B[2], Ez);
      for (int h = 0; h < nhL; h++) {
         int o = hoff[L][h], t = o/(DD*DD), u = (o/DD)%DD, v = o%DD;
         double s = 0.0;
         for (int bb = 0; bb < ncj; bb++)
            for (int aa = 0; aa < nci; aa++) {
               if (t > ix[aa]+jx[bb] || u > iy[aa]+jy[bb] || v > iz[aa]+jz[bb]) continue;
               s += Dc[bb*nci + aa]*Ex[EI(ix[aa],jx[bb],t)]*Ey[EI(iy[aa],jy[bb],u)]*Ez[EI(iz[aa],jz[bb],v)];
            }
         Dh[h] = Kab*s;
      }
      for (int k = 0; k < naux; k++) {
         if (*schw*auxb[k] < *cut) continue;
         int lk = ax_l[k], nhk = NHERM(lk), n = L + lk;
         for (int g = 0; g < ax_nprim[k]; g++) {
            int q = ax_first[k] + g;
            double gam = ax_exp[q], pg = p + gam, al = p*gam/pg;
            double X = Pc[0] - ax_cen[3*q], Y = Pc[1] - ax_cen[3*q+1], Z = Pc[2] - ax_cen[3*q+2];
            rtab(n, al, X, Y, Z, al*(X*X + Y*Y + Z*Z), R, wk);
            double lam = 34.986836655249725693/(p*gam*sqrt(pg));
            double *Gq = G + ax_goff[k] + (long)g*nhk;
            for (int h2 = 0; h2 < nhk; h2++) {
               const double *Rs = R + hoff[lk][h2];
               double s = 0.0;
               for (int h = 0; h < nhL; h++) s += Dh[h]*Rs[hoff[L][h]];
               Gq[h2] += lam*s;
            }
         }
      }
   }
}

/* bvec(aux functions) += A^T G per aux primitive */
void mdj_bfinal_ref(const double *G, double *bvec)
{
   for (int k = 0; k < naux; k++) {
      int l = ax_l[k], nh = NHERM(l), ns = nfun[l];
      for (int g = 0; g < ax_nprim[k]; g++) {
         const double *A = ax_A + aoff_of(k, g), *Gq = G + ax_goff[k] + (long)g*nh;
         for (int s = 0; s < ns; s++) {
            double v = 0.0;
            for (int h = 0; h < nh; h++) v += A[s*nh + h]*Gq[h];
            bvec[ax_foff[k] + s] += v;
         }
      }
   }
}

/* H[aux hermite slots] = A c (per aux primitive) */
void mdj_jprep_ref(const double *c, double *H)
{
   for (int k = 0; k < naux; k++) {
      int l = ax_l[k], nh = NHERM(l), ns = nfun[l];
      for (int g = 0; g < ax_nprim[k]; g++) {
         const double *A = ax_A + aoff_of(k, g);
         double *Hq = H + ax_goff[k] + (long)g*nh;
         for (int h = 0; h < nh; h++) {
            double v = 0.0;
            for (int s = 0; s < ns; s++) v += A[s*nh + h]*c[ax_foff[k] + s];
            Hq[h] = v;
         }
      }
   }
}

/* J-pass for one shell pair: Jsph (di x dj, ld di) = sum_K (mu nu|K) c_K, aux shells kept when schw*cb[k] >= cut */
void mdj_jpass_pair_ref(const int *ish_, const int *jsh_, const double *H, const double *schw,
                    const double *auxb, const double *cut0, const double *cb, const double *cut, double *Jsph)
{
   pair_t P; pair_setup(*ish_, *jsh_, &P);
   int li = P.li, lj = P.lj, L = li + lj, nci = NCART(li), ncj = NCART(lj), nhL = NHERM(L);
   int di = nfun[li], dj = nfun[lj];
   double Jc[NCART(LMAXO)*NCART(LMAXO)], Ex[EISZ],
          Ey[EISZ], Ez[EISZ], V[NHERM(LMAXP)];
   double R[DD*DD*DD], wk[DD*DD*DD];
   memset(Jc, 0, sizeof(double)*nci*ncj);
   double AB2 = (P.A[0]-P.B[0])*(P.A[0]-P.B[0]) + (P.A[1]-P.B[1])*(P.A[1]-P.B[1]) + (P.A[2]-P.B[2])*(P.A[2]-P.B[2]);
   int ix[NCART(LMAXO)], iy[NCART(LMAXO)], iz[NCART(LMAXO)], jx[NCART(LMAXO)], jy[NCART(LMAXO)], jz[NCART(LMAXO)];
   for (int a = 0; a < nci; a++) cart_xyz(li, a, ix+a, iy+a, iz+a);
   for (int b = 0; b < ncj; b++) cart_xyz(lj, b, jx+b, jy+b, jz+b);
   for (int pa = 0; pa < P.ni; pa++)
   for (int pb = 0; pb < P.nj; pb++) {
      double a = P.ai[pa], b = P.aj[pb], p = a + b;
      double Kab = exp(-a*b/p*AB2)*P.ci[pa]*P.cj[pb];
      if (fabs(Kab)*pow(M_PI/p, 1.5) < SCREEN_PRIM) continue;
      double Pc[3];
      for (int d = 0; d < 3; d++) Pc[d] = (a*P.A[d] + b*P.B[d])/p;
      memset(V, 0, sizeof(double)*nhL);
      for (int k = 0; k < naux; k++) {
         if (*schw*cb[k] < *cut || *schw*auxb[k] < *cut0) continue;
         int lk = ax_l[k], nhk = NHERM(lk), n = L + lk;
         for (int g = 0; g < ax_nprim[k]; g++) {
            int q = ax_first[k] + g;
            double gam = ax_exp[q], pg = p + gam, al = p*gam/pg;
            double X = Pc[0] - ax_cen[3*q], Y = Pc[1] - ax_cen[3*q+1], Z = Pc[2] - ax_cen[3*q+2];
            rtab(n, al, X, Y, Z, al*(X*X + Y*Y + Z*Z), R, wk);
            double lam = 34.986836655249725693/(p*gam*sqrt(pg));
            const double *Hq = H + ax_goff[k] + (long)g*nhk;
            for (int h2 = 0; h2 < nhk; h2++) {
               double w = lam*Hq[h2];
               if (w == 0.0) continue;
               const double *Rs = R + hoff[lk][h2];
               for (int h = 0; h < nhL; h++) V[h] += w*Rs[hoff[L][h]];
            }
         }
      }
      herm1d(li, lj, p, Pc[0]-P.A[0], Pc[0]-P.B[0], Ex);
      herm1d(li, lj, p, Pc[1]-P.A[1], Pc[1]-P.B[1], Ey);
      herm1d(li, lj, p, Pc[2]-P.A[2], Pc[2]-P.B[2], Ez);
      for (int bb = 0; bb < ncj; bb++)
         for (int aa = 0; aa < nci; aa++) {
            double s = 0.0;
            for (int h = 0; h < nhL; h++) {
               int o = hoff[L][h], t = o/(DD*DD), u = (o/DD)%DD, v = o%DD;
               if (t > ix[aa]+jx[bb] || u > iy[aa]+jy[bb] || v > iz[aa]+jz[bb]) continue;
               s += V[h]*Ex[EI(ix[aa],jx[bb],t)]*Ey[EI(iy[aa],jy[bb],u)]*Ez[EI(iz[aa],jz[bb],v)];
            }
            Jc[bb*nci + aa] += Kab*s;
         }
   }
   /* Jsph = c2s_i Jc c2s_j^T */
   double tmp[(2*LMAXO+1)*NCART(LMAXO)];
   for (int bb = 0; bb < ncj; bb++)
      for (int m = 0; m < di; m++) {
         double s = 0.0;
         for (int aa = 0; aa < nci; aa++) s += c2s[li][m*nci + aa]*Jc[bb*nci + aa];
         tmp[bb*di + m] = s;
      }
   for (int n = 0; n < dj; n++)
      for (int m = 0; m < di; m++) {
         double s = 0.0;
         for (int bb = 0; bb < ncj; bb++) s += tmp[bb*di + m]*c2s[lj][n*ncj + bb];
         Jsph[n*di + m] = s;
      }
}

/* ------------------------------------------------------------------------
 * Vectorised form: aux primitives grouped by l, NB primitives per SIMD
 * block, compact Hermite index lists. G/H layout: per l group, Hermite
 * component major, primitive minor (gbase[l] + h*nqp[l] + q).
 * ---------------------------------------------------------------------- */
#define NB 8
static int lmaxa_used, nq[LMAXA + 1], nqp[LMAXA + 1];
static long gbase[LMAXA + 1], gsize2;
static double *q_exp[LMAXA + 1], *q_cx[LMAXA + 1], *q_cy[LMAXA + 1], *q_cz[LMAXA + 1];
static int *q_shell[LMAXA + 1], *q_prim[LMAXA + 1], *q_atom[LMAXA + 1];
static int cid[DD][DD][DD];
static int cdir[NHERM(NMAX)], cp1[NHERM(NMAX)], cp2[NHERM(NMAX)], cmul[NHERM(NMAX)];
static int *sumtab[LMAXP + 3][LMAXA + 1];

static long v2_init(void)
{
   for (int k = 0; k < NHERM(NMAX); k++) {
      int o = hoff[NMAX][k], t = o/(DD*DD), u = (o/DD)%DD, v = o%DD;
      cid[t][u][v] = k;
   }
   for (int k = 1; k < NHERM(NMAX); k++) {
      int o = hoff[NMAX][k], t = o/(DD*DD), u = (o/DD)%DD, v = o%DD;
      if (t > 0) { cdir[k] = 0; cmul[k] = t-1; cp1[k] = cid[t-1][u][v]; cp2[k] = t > 1 ? cid[t-2][u][v] : 0; }
      else if (u > 0) { cdir[k] = 1; cmul[k] = u-1; cp1[k] = cid[t][u-1][v]; cp2[k] = u > 1 ? cid[t][u-2][v] : 0; }
      else { cdir[k] = 2; cmul[k] = v-1; cp1[k] = cid[t][u][v-1]; cp2[k] = v > 1 ? cid[t][u][v-2] : 0; }
   }
   for (int L = 0; L <= LMAXP + 2; L++)
      for (int lk = 0; lk <= LMAXA; lk++) {
         free(sumtab[L][lk]);
         int nh = NHERM(L), nk = NHERM(lk);
         sumtab[L][lk] = malloc(sizeof(int)*nh*nk);
         for (int h = 0; h < nh; h++)
            for (int h2 = 0; h2 < nk; h2++) {
               int o1 = hoff[NMAX][h], o2 = hoff[NMAX][h2];
               int t = o1/(DD*DD) + o2/(DD*DD), u = (o1/DD)%DD + (o2/DD)%DD, v = o1%DD + o2%DD;
               sumtab[L][lk][h*nk + h2] = cid[t][u][v];
            }
      }
   lmaxa_used = 0;
   for (int l = 0; l <= LMAXA; l++) nq[l] = 0;
   for (int k = 0; k < naux; k++) { nq[ax_l[k]] += ax_nprim[k]; if (ax_l[k] > lmaxa_used) lmaxa_used = ax_l[k]; }
   gsize2 = 0;
   for (int l = 0; l <= LMAXA; l++) {
      nqp[l] = (nq[l] + NB - 1)/NB*NB;
      gbase[l] = gsize2;
      gsize2 += (long)NHERM(l)*nqp[l];
      free(q_exp[l]); free(q_cx[l]); free(q_cy[l]); free(q_cz[l]); free(q_shell[l]); free(q_prim[l]); free(q_atom[l]);
      q_exp[l] = malloc(sizeof(double)*(nqp[l] + 1)); q_cx[l] = malloc(sizeof(double)*(nqp[l] + 1));
      q_cy[l] = malloc(sizeof(double)*(nqp[l] + 1)); q_cz[l] = malloc(sizeof(double)*(nqp[l] + 1));
      q_shell[l] = malloc(sizeof(int)*(nqp[l] + 1)); q_prim[l] = malloc(sizeof(int)*(nqp[l] + 1));
      q_atom[l] = malloc(sizeof(int)*(nqp[l] + 1));
      int q = 0;
      for (int k = 0; k < naux; k++) {
         if (ax_l[k] != l) continue;
         for (int g = 0; g < ax_nprim[k]; g++) {
            int r = ax_first[k] + g;
            q_exp[l][q] = ax_exp[r]; q_cx[l][q] = ax_cen[3*r]; q_cy[l][q] = ax_cen[3*r+1]; q_cz[l][q] = ax_cen[3*r+2];
            q_shell[l][q] = k; q_prim[l][q] = g; q_atom[l][q] = g_bas[8*(norb + k)]; q++;
         }
      }
      for (; q < nqp[l]; q++) {   /* padding: far away, weight zero */
         q_exp[l][q] = 1.0; q_cx[l][q] = 1.0e3; q_cy[l][q] = 1.0e3; q_cz[l][q] = 1.0e3;
         q_shell[l][q] = -1; q_prim[l][q] = 0; q_atom[l][q] = -1;
      }
   }
   return gsize2;
}

/* F_0..F_n for NB lanes, F[m*NB + v] */
static inline void boys_lanes(int n, const double *T, double *F)
{
   for (int v = 0; v < NB; v++) {
      double t = T[v];
      if (t >= BTMAX) {
         double b = 0.5/t;
         F[v] = 0.886226925452758013649/sqrt(t);
         for (int i = 1; i <= n; i++) F[i*NB + v] = b*(2*i-1)*F[(i-1)*NB + v];
      } else {
         int k = (int)(t/BH + 0.5);
         double d = k*BH - t;
         const double *r = btab + (long)k*BNN + n;
         double fn = r[0] + d*(r[1] + d*(r[2]/2 + d*(r[3]/6 + d*(r[4]/24 + d*(r[5]/120 + d*r[6]/720)))));
         F[n*NB + v] = fn;
         if (n > 0) {
            double e = exp(-t);
            for (int i = n; i > 0; i--) F[(i-1)*NB + v] = (2*t*F[i*NB + v] + e)/(2*i-1);
         }
      }
   }
}

/* R (compact index list, NB lanes) up to order n; returns the buffer holding level 0 */
static inline double *rtab_lanes(int n, const double *al, const double *X, const double *Y, const double *Z,
                                 const double *F, double *B0, double *B1)
{
   double pw[NB];
   for (int v = 0; v < NB; v++) pw[v] = 1.0;
   if (n == 0) {
      for (int v = 0; v < NB; v++) B0[v] = F[v];
      return B0;
   }
   double m2a[NB], pwl[(NMAX + 1)*NB];
   for (int v = 0; v < NB; v++) { m2a[v] = -2.0*al[v]; pwl[v] = 1.0; }
   for (int i = 1; i <= n; i++)
      for (int v = 0; v < NB; v++) pwl[i*NB + v] = pwl[(i-1)*NB + v]*m2a[v];
   double *cur = (n & 1) ? B1 : B0;
   for (int v = 0; v < NB; v++) cur[v] = F[n*NB + v]*pwl[n*NB + v];
   const double *XYZ[3] = { X, Y, Z };
   for (int j = 1; j <= n; j++) {
      int lev = n - j;
      double *nx = (cur == B0) ? B1 : B0;
      const double *pv = cur;
      for (int v = 0; v < NB; v++) nx[v] = F[lev*NB + v]*pwl[lev*NB + v];
      int nk = NHERM(j);
      for (int k = 1; k < nk; k++) {
         const double *xd = XYZ[cdir[k]], *a1 = pv + cp1[k]*NB;
         double *o = nx + k*NB;
         if (cmul[k]) {
            const double *a2 = pv + cp2[k]*NB; double m = cmul[k];
#pragma omp simd
            for (int v = 0; v < NB; v++) o[v] = xd[v]*a1[v] + m*a2[v];
         } else {
#pragma omp simd
            for (int v = 0; v < NB; v++) o[v] = xd[v]*a1[v];
         }
      }
      cur = nx;
   }
   (void)pw;
   return cur;
}

#define LAMC 34.986836655249725693   /* 2 pi^{5/2} */

static void pair_hermite_setup(const pair_t *P, int pa, int pb, double AB2, double *p_, double *Kab_, double *Pc)
{
   double a = P->ai[pa], b = P->aj[pb], p = a + b;
   *p_ = p;
   *Kab_ = exp(-a*b/p*AB2)*P->ci[pa]*P->cj[pb];
   for (int d = 0; d < 3; d++) Pc[d] = (a*P->A[d] + b*P->B[d])/p;
}

/* lane block setup: T, lambda (masked), al, X, Y, Z */
static inline void lanes_setup(int lk, int q0, double p, const double *Pc, double schw,
                               const double *bnd1, double cut1, const double *bnd2, double cut2,
                               double *al, double *X, double *Y, double *Z, double *T, double *lam)
{
   for (int v = 0; v < NB; v++) {
      int q = q0 + v, k = q_shell[lk][q];
      double g = q_exp[lk][q], pg = p + g;
      al[v] = p*g/pg;
      X[v] = Pc[0] - q_cx[lk][q]; Y[v] = Pc[1] - q_cy[lk][q]; Z[v] = Pc[2] - q_cz[lk][q];
      T[v] = al[v]*(X[v]*X[v] + Y[v]*Y[v] + Z[v]*Z[v]);
      int keep = k >= 0 && schw*bnd1[k >= 0 ? k : 0] >= cut1 && (!bnd2 || schw*bnd2[k >= 0 ? k : 0] >= cut2);
      lam[v] = keep ? LAMC/(p*g*sqrt(pg)) : 0.0;
   }
}

void mdj_bpass_pair(const int *ish_, const int *jsh_, const double *D, const double *schw,
                    const double *auxb, const double *cut, double *G)
{
   pair_t P; pair_setup(*ish_, *jsh_, &P);
   int li = P.li, lj = P.lj, L = li + lj, nci = NCART(li), ncj = NCART(lj), nhL = NHERM(L);
   double Dc[NCART(LMAXO)*NCART(LMAXO)], Ex[EISZ],
          Ey[EISZ], Ez[EISZ], Dh[NHERM(LMAXP)];
   double B0[NHERM(NMAX)*NB], B1[NHERM(NMAX)*NB], F[(NMAX + 1)*NB];
   double al[NB], X[NB], Y[NB], Z[NB], T[NB], lam[NB], acc[NB];
   to_cart(li, lj, D, nfun[li], Dc);
   double dmax = 0.0;
   for (int k = 0; k < nci*ncj; k++) dmax = fmax(dmax, fabs(Dc[k]));
   if (dmax == 0.0) return;
   double AB2 = (P.A[0]-P.B[0])*(P.A[0]-P.B[0]) + (P.A[1]-P.B[1])*(P.A[1]-P.B[1]) + (P.A[2]-P.B[2])*(P.A[2]-P.B[2]);
   int ix[NCART(LMAXO)], iy[NCART(LMAXO)], iz[NCART(LMAXO)], jx[NCART(LMAXO)], jy[NCART(LMAXO)], jz[NCART(LMAXO)];
   for (int a = 0; a < nci; a++) cart_xyz(li, a, ix+a, iy+a, iz+a);
   for (int b = 0; b < ncj; b++) cart_xyz(lj, b, jx+b, jy+b, jz+b);
   for (int pa = 0; pa < P.ni; pa++)
   for (int pb = 0; pb < P.nj; pb++) {
      double p, Kab, Pc[3];
      pair_hermite_setup(&P, pa, pb, AB2, &p, &Kab, Pc);
      if (fabs(Kab)*dmax*pow(M_PI/p, 1.5) < SCREEN_PRIM) continue;
      herm1d(li, lj, p, Pc[0]-P.A[0], Pc[0]-P.B[0], Ex);
      herm1d(li, lj, p, Pc[1]-P.A[1], Pc[1]-P.B[1], Ey);
      herm1d(li, lj, p, Pc[2]-P.A[2], Pc[2]-P.B[2], Ez);
      for (int h = 0; h < nhL; h++) {
         int o = hoff[L][h], t = o/(DD*DD), u = (o/DD)%DD, v = o%DD;
         double s = 0.0;
         for (int bb = 0; bb < ncj; bb++)
            for (int aa = 0; aa < nci; aa++) {
               if (t > ix[aa]+jx[bb] || u > iy[aa]+jy[bb] || v > iz[aa]+jz[bb]) continue;
               s += Dc[bb*nci + aa]*Ex[EI(ix[aa],jx[bb],t)]*Ey[EI(iy[aa],jy[bb],u)]*Ez[EI(iz[aa],jz[bb],v)];
            }
         Dh[h] = Kab*s;
      }
      for (int lk = 0; lk <= lmaxa_used; lk++) {
         if (nq[lk] == 0) continue;
         int n = L + lk, nhk = NHERM(lk);
         const int *S = sumtab[L][lk];
         for (int q0 = 0; q0 < nqp[lk]; q0 += NB) {
            lanes_setup(lk, q0, p, Pc, *schw, auxb, *cut, NULL, 0.0, al, X, Y, Z, T, lam);
            boys_lanes(n, T, F);
            double *R = rtab_lanes(n, al, X, Y, Z, F, B0, B1);
            for (int h2 = 0; h2 < nhk; h2++) {
               for (int v = 0; v < NB; v++) acc[v] = 0.0;
               for (int h = 0; h < nhL; h++) {
                  const double *r = R + S[h*nhk + h2]*NB; double d = Dh[h];
#pragma omp simd
                  for (int v = 0; v < NB; v++) acc[v] += d*r[v];
               }
               double *Gp = G + gbase[lk] + (long)h2*nqp[lk] + q0;
#pragma omp simd
               for (int v = 0; v < NB; v++) Gp[v] += lam[v]*acc[v];
            }
         }
      }
   }
}

void mdj_bfinal(const double *G, double *bvec)
{
   for (int l = 0; l <= lmaxa_used; l++) {
      int nh = NHERM(l), ns = nfun[l];
      for (int q = 0; q < nq[l]; q++) {
         int k = q_shell[l][q], g = q_prim[l][q];
         const double *A = ax_A + aoff_of(k, g);
         for (int s = 0; s < ns; s++) {
            double v = 0.0;
            for (int h = 0; h < nh; h++) v += A[s*nh + h]*G[gbase[l] + (long)h*nqp[l] + q];
            bvec[ax_foff[k] + s] += v;
         }
      }
   }
}

void mdj_jprep(const double *c, double *H)
{
   memset(H, 0, sizeof(double)*gsize2);
   for (int l = 0; l <= lmaxa_used; l++) {
      int nh = NHERM(l), ns = nfun[l];
      for (int q = 0; q < nq[l]; q++) {
         int k = q_shell[l][q], g = q_prim[l][q];
         const double *A = ax_A + aoff_of(k, g);
         for (int h = 0; h < nh; h++) {
            double v = 0.0;
            for (int s = 0; s < ns; s++) v += A[s*nh + h]*c[ax_foff[k] + s];
            H[gbase[l] + (long)h*nqp[l] + q] = v;
         }
      }
   }
}

void mdj_jpass_pair(const int *ish_, const int *jsh_, const double *H, const double *schw,
                    const double *auxb, const double *cut0, const double *cb, const double *cut, double *Jsph)
{
   pair_t P; pair_setup(*ish_, *jsh_, &P);
   int li = P.li, lj = P.lj, L = li + lj, nci = NCART(li), ncj = NCART(lj), nhL = NHERM(L);
   int di = nfun[li], dj = nfun[lj];
   double Jc[NCART(LMAXO)*NCART(LMAXO)], Ex[EISZ],
          Ey[EISZ], Ez[EISZ], V[NHERM(LMAXP)];
   double B0[NHERM(NMAX)*NB], B1[NHERM(NMAX)*NB], F[(NMAX + 1)*NB], W[NHERM(LMAXA)*NB];
   double al[NB], X[NB], Y[NB], Z[NB], T[NB], lam[NB], acc[NB];
   memset(Jc, 0, sizeof(double)*nci*ncj);
   double AB2 = (P.A[0]-P.B[0])*(P.A[0]-P.B[0]) + (P.A[1]-P.B[1])*(P.A[1]-P.B[1]) + (P.A[2]-P.B[2])*(P.A[2]-P.B[2]);
   int ix[NCART(LMAXO)], iy[NCART(LMAXO)], iz[NCART(LMAXO)], jx[NCART(LMAXO)], jy[NCART(LMAXO)], jz[NCART(LMAXO)];
   for (int a = 0; a < nci; a++) cart_xyz(li, a, ix+a, iy+a, iz+a);
   for (int b = 0; b < ncj; b++) cart_xyz(lj, b, jx+b, jy+b, jz+b);
   for (int pa = 0; pa < P.ni; pa++)
   for (int pb = 0; pb < P.nj; pb++) {
      double p, Kab, Pc[3];
      pair_hermite_setup(&P, pa, pb, AB2, &p, &Kab, Pc);
      if (fabs(Kab)*pow(M_PI/p, 1.5) < SCREEN_PRIM) continue;
      memset(V, 0, sizeof(double)*nhL);
      for (int lk = 0; lk <= lmaxa_used; lk++) {
         if (nq[lk] == 0) continue;
         int n = L + lk, nhk = NHERM(lk);
         const int *S = sumtab[L][lk];
         for (int q0 = 0; q0 < nqp[lk]; q0 += NB) {
            lanes_setup(lk, q0, p, Pc, *schw, cb, *cut, auxb, *cut0, al, X, Y, Z, T, lam);
            int any = 0;
            for (int v = 0; v < NB; v++) any |= (lam[v] != 0.0);
            if (!any) continue;
            boys_lanes(n, T, F);
            double *R = rtab_lanes(n, al, X, Y, Z, F, B0, B1);
            for (int h2 = 0; h2 < nhk; h2++) {
               const double *Hp = H + gbase[lk] + (long)h2*nqp[lk] + q0;
               for (int v = 0; v < NB; v++) W[h2*NB + v] = lam[v]*Hp[v];
            }
            for (int h = 0; h < nhL; h++) {
               for (int v = 0; v < NB; v++) acc[v] = 0.0;
               for (int h2 = 0; h2 < nhk; h2++) {
                  const double *r = R + S[h*nhk + h2]*NB, *w = W + h2*NB;
#pragma omp simd
                  for (int v = 0; v < NB; v++) acc[v] += w[v]*r[v];
               }
               double s = 0.0;
               for (int v = 0; v < NB; v++) s += acc[v];
               V[h] += s;
            }
         }
      }
      herm1d(li, lj, p, Pc[0]-P.A[0], Pc[0]-P.B[0], Ex);
      herm1d(li, lj, p, Pc[1]-P.A[1], Pc[1]-P.B[1], Ey);
      herm1d(li, lj, p, Pc[2]-P.A[2], Pc[2]-P.B[2], Ez);
      for (int bb = 0; bb < ncj; bb++)
         for (int aa = 0; aa < nci; aa++) {
            double s = 0.0;
            for (int h = 0; h < nhL; h++) {
               int o = hoff[L][h], t = o/(DD*DD), u = (o/DD)%DD, v = o%DD;
               if (t > ix[aa]+jx[bb] || u > iy[aa]+jy[bb] || v > iz[aa]+jz[bb]) continue;
               s += V[h]*Ex[EI(ix[aa],jx[bb],t)]*Ey[EI(iy[aa],jy[bb],u)]*Ez[EI(iz[aa],jz[bb],v)];
            }
            Jc[bb*nci + aa] += Kab*s;
         }
   }
   double tmp[(2*LMAXO+1)*NCART(LMAXO)];
   for (int bb = 0; bb < ncj; bb++)
      for (int m = 0; m < di; m++) {
         double s = 0.0;
         for (int aa = 0; aa < nci; aa++) s += c2s[li][m*nci + aa]*Jc[bb*nci + aa];
         tmp[bb*di + m] = s;
      }
   for (int n2 = 0; n2 < dj; n2++)
      for (int m = 0; m < di; m++) {
         double s = 0.0;
         for (int bb = 0; bb < ncj; bb++) s += tmp[bb*di + m]*c2s[lj][n2*ncj + bb];
         Jsph[n2*di + m] = s;
      }
}

/* ------------------------------------------------------------------------
 * Several densities at once (TDDFT / CPSCF batches): one R table per
 * (primitive pair, aux block) serves all nb densities.
 *   b-pass: D (di*dj x nb, leading dimension ldD), G (gsize2 x nb, ldG)
 *   J-pass: H (gsize2 x nb, ldH) -> J (di*dj x nb, ldJ); cb is the aux
 *   shell bound of c maximised over the batch
 * ---------------------------------------------------------------------- */
void mdj_bpass_multi(const int *ish_, const int *jsh_, const double *D, const int *ldD, const int *nb_,
                     const double *schw, const double *auxb, const double *cut, double *G, const long *ldG)
{
   pair_t P; pair_setup(*ish_, *jsh_, &P);
   int nb = *nb_;
   int li = P.li, lj = P.lj, L = li + lj, nci = NCART(li), ncj = NCART(lj), nhL = NHERM(L);
   double Ex[EISZ], Ey[EISZ], Ez[EISZ];
   double B0[NHERM(NMAX)*NB], B1[NHERM(NMAX)*NB], F[(NMAX + 1)*NB];
   double al[NB], X[NB], Y[NB], Z[NB], T[NB], lam[NB];
   double *Dc = malloc(sizeof(double)*nci*ncj*nb), *Dh = malloc(sizeof(double)*nhL*nb),
          *acc = malloc(sizeof(double)*NB*nb);
   double dmax = 0.0;
   for (int b = 0; b < nb; b++) {
      to_cart(li, lj, D + (long)b*(*ldD), nfun[li], Dc + (long)b*nci*ncj);
      for (int k = 0; k < nci*ncj; k++) dmax = fmax(dmax, fabs(Dc[(long)b*nci*ncj + k]));
   }
   if (dmax == 0.0) { free(Dc); free(Dh); free(acc); return; }
   double AB2 = (P.A[0]-P.B[0])*(P.A[0]-P.B[0]) + (P.A[1]-P.B[1])*(P.A[1]-P.B[1]) + (P.A[2]-P.B[2])*(P.A[2]-P.B[2]);
   int ix[NCART(LMAXO)], iy[NCART(LMAXO)], iz[NCART(LMAXO)], jx[NCART(LMAXO)], jy[NCART(LMAXO)], jz[NCART(LMAXO)];
   for (int a = 0; a < nci; a++) cart_xyz(li, a, ix+a, iy+a, iz+a);
   for (int b = 0; b < ncj; b++) cart_xyz(lj, b, jx+b, jy+b, jz+b);
   for (int pa = 0; pa < P.ni; pa++)
   for (int pb = 0; pb < P.nj; pb++) {
      double p, Kab, Pc[3];
      pair_hermite_setup(&P, pa, pb, AB2, &p, &Kab, Pc);
      if (fabs(Kab)*dmax*pow(M_PI/p, 1.5) < SCREEN_PRIM) continue;
      herm1d(li, lj, p, Pc[0]-P.A[0], Pc[0]-P.B[0], Ex);
      herm1d(li, lj, p, Pc[1]-P.A[1], Pc[1]-P.B[1], Ey);
      herm1d(li, lj, p, Pc[2]-P.A[2], Pc[2]-P.B[2], Ez);
      for (int h = 0; h < nhL; h++) {
         int o = hoff[L][h], t = o/(DD*DD), u = (o/DD)%DD, v = o%DD;
         for (int b = 0; b < nb; b++) Dh[(long)b*nhL + h] = 0.0;
         for (int bb = 0; bb < ncj; bb++)
            for (int aa = 0; aa < nci; aa++) {
               if (t > ix[aa]+jx[bb] || u > iy[aa]+jy[bb] || v > iz[aa]+jz[bb]) continue;
               double e = Kab*Ex[EI(ix[aa],jx[bb],t)]*Ey[EI(iy[aa],jy[bb],u)]*Ez[EI(iz[aa],jz[bb],v)];
               for (int b = 0; b < nb; b++) Dh[(long)b*nhL + h] += e*Dc[(long)b*nci*ncj + bb*nci + aa];
            }
      }
      for (int lk = 0; lk <= lmaxa_used; lk++) {
         if (nq[lk] == 0) continue;
         int n = L + lk, nhk = NHERM(lk);
         const int *S = sumtab[L][lk];
         for (int q0 = 0; q0 < nqp[lk]; q0 += NB) {
            lanes_setup(lk, q0, p, Pc, *schw, auxb, *cut, NULL, 0.0, al, X, Y, Z, T, lam);
            boys_lanes(n, T, F);
            double *R = rtab_lanes(n, al, X, Y, Z, F, B0, B1);
            for (int h2 = 0; h2 < nhk; h2++) {
               memset(acc, 0, sizeof(double)*NB*nb);
               for (int h = 0; h < nhL; h++) {
                  const double *r = R + S[h*nhk + h2]*NB;
                  for (int b = 0; b < nb; b++) {
                     double d = Dh[(long)b*nhL + h], *ab = acc + (long)b*NB;
#pragma omp simd
                     for (int v = 0; v < NB; v++) ab[v] += d*r[v];
                  }
               }
               for (int b = 0; b < nb; b++) {
                  double *Gp = G + (long)b*(*ldG) + gbase[lk] + (long)h2*nqp[lk] + q0, *ab = acc + (long)b*NB;
#pragma omp simd
                  for (int v = 0; v < NB; v++) Gp[v] += lam[v]*ab[v];
               }
            }
         }
      }
   }
   free(Dc); free(Dh); free(acc);
}

void mdj_jpass_multi(const int *ish_, const int *jsh_, const double *H, const long *ldH, const int *nb_,
                     const double *schw, const double *auxb, const double *cut0, const double *cb,
                     const double *cut, double *Jsph, const int *ldJ)
{
   pair_t P; pair_setup(*ish_, *jsh_, &P);
   int nb = *nb_;
   int li = P.li, lj = P.lj, L = li + lj, nci = NCART(li), ncj = NCART(lj), nhL = NHERM(L);
   int di = nfun[li], dj = nfun[lj];
   double Ex[EISZ], Ey[EISZ], Ez[EISZ];
   double B0[NHERM(NMAX)*NB], B1[NHERM(NMAX)*NB], F[(NMAX + 1)*NB];
   double al[NB], X[NB], Y[NB], Z[NB], T[NB], lam[NB];
   double *Jc = calloc((size_t)nci*ncj*nb, sizeof(double)), *V = malloc(sizeof(double)*nhL*nb),
          *W = malloc(sizeof(double)*NHERM(LMAXA)*NB*nb);
   double AB2 = (P.A[0]-P.B[0])*(P.A[0]-P.B[0]) + (P.A[1]-P.B[1])*(P.A[1]-P.B[1]) + (P.A[2]-P.B[2])*(P.A[2]-P.B[2]);
   int ix[NCART(LMAXO)], iy[NCART(LMAXO)], iz[NCART(LMAXO)], jx[NCART(LMAXO)], jy[NCART(LMAXO)], jz[NCART(LMAXO)];
   for (int a = 0; a < nci; a++) cart_xyz(li, a, ix+a, iy+a, iz+a);
   for (int b = 0; b < ncj; b++) cart_xyz(lj, b, jx+b, jy+b, jz+b);
   for (int pa = 0; pa < P.ni; pa++)
   for (int pb = 0; pb < P.nj; pb++) {
      double p, Kab, Pc[3];
      pair_hermite_setup(&P, pa, pb, AB2, &p, &Kab, Pc);
      if (fabs(Kab)*pow(M_PI/p, 1.5) < SCREEN_PRIM) continue;
      memset(V, 0, sizeof(double)*nhL*nb);
      for (int lk = 0; lk <= lmaxa_used; lk++) {
         if (nq[lk] == 0) continue;
         int n = L + lk, nhk = NHERM(lk);
         const int *S = sumtab[L][lk];
         for (int q0 = 0; q0 < nqp[lk]; q0 += NB) {
            lanes_setup(lk, q0, p, Pc, *schw, cb, *cut, auxb, *cut0, al, X, Y, Z, T, lam);
            int any = 0;
            for (int v = 0; v < NB; v++) any |= (lam[v] != 0.0);
            if (!any) continue;
            boys_lanes(n, T, F);
            double *R = rtab_lanes(n, al, X, Y, Z, F, B0, B1);
            for (int b = 0; b < nb; b++)
               for (int h2 = 0; h2 < nhk; h2++) {
                  const double *Hp = H + (long)b*(*ldH) + gbase[lk] + (long)h2*nqp[lk] + q0;
                  double *w = W + ((long)b*nhk + h2)*NB;
                  for (int v = 0; v < NB; v++) w[v] = lam[v]*Hp[v];
               }
            for (int h = 0; h < nhL; h++) {
               for (int b = 0; b < nb; b++) {
                  double acc[NB];
                  for (int v = 0; v < NB; v++) acc[v] = 0.0;
                  for (int h2 = 0; h2 < nhk; h2++) {
                     const double *r = R + S[h*nhk + h2]*NB, *w = W + ((long)b*nhk + h2)*NB;
#pragma omp simd
                     for (int v = 0; v < NB; v++) acc[v] += w[v]*r[v];
                  }
                  double s = 0.0;
                  for (int v = 0; v < NB; v++) s += acc[v];
                  V[(long)b*nhL + h] += s;
               }
            }
         }
      }
      herm1d(li, lj, p, Pc[0]-P.A[0], Pc[0]-P.B[0], Ex);
      herm1d(li, lj, p, Pc[1]-P.A[1], Pc[1]-P.B[1], Ey);
      herm1d(li, lj, p, Pc[2]-P.A[2], Pc[2]-P.B[2], Ez);
      for (int bb = 0; bb < ncj; bb++)
         for (int aa = 0; aa < nci; aa++)
            for (int h = 0; h < nhL; h++) {
               int o = hoff[L][h], t = o/(DD*DD), u = (o/DD)%DD, v = o%DD;
               if (t > ix[aa]+jx[bb] || u > iy[aa]+jy[bb] || v > iz[aa]+jz[bb]) continue;
               double e = Kab*Ex[EI(ix[aa],jx[bb],t)]*Ey[EI(iy[aa],jy[bb],u)]*Ez[EI(iz[aa],jz[bb],v)];
               for (int b = 0; b < nb; b++) Jc[(long)b*nci*ncj + bb*nci + aa] += e*V[(long)b*nhL + h];
            }
   }
   double tmp[(2*LMAXO+1)*NCART(LMAXO)];
   for (int b = 0; b < nb; b++) {
      const double *Jb = Jc + (long)b*nci*ncj;
      double *Js = Jsph + (long)b*(*ldJ);
      for (int bb = 0; bb < ncj; bb++)
         for (int m = 0; m < di; m++) {
            double s = 0.0;
            for (int aa = 0; aa < nci; aa++) s += c2s[li][m*nci + aa]*Jb[bb*nci + aa];
            tmp[bb*di + m] = s;
         }
      for (int n2 = 0; n2 < dj; n2++)
         for (int m = 0; m < di; m++) {
            double s = 0.0;
            for (int bb = 0; bb < ncj; bb++) s += tmp[bb*di + m]*c2s[lj][n2*ncj + bb];
            Js[n2*di + m] = s;
         }
   }
   free(Jc); free(V); free(W);
}

/* sum_t E^x_{ix jx tx} E^y E^z V_{t} over the Hermite functions of the product */
static inline double ev3(const int *i3, const int *j3, const double *Ex, const double *Ey, const double *Ez, const double *V)
{
   double s = 0.0;
   for (int tx = 0; tx <= i3[0] + j3[0]; tx++) {
      double ex = Ex[EI(i3[0], j3[0], tx)];
      if (ex == 0.0) continue;
      for (int ty = 0; ty <= i3[1] + j3[1]; ty++) {
         double ey = ex*Ey[EI(i3[1], j3[1], ty)];
         if (ey == 0.0) continue;
         for (int tz = 0; tz <= i3[2] + j3[2]; tz++)
            s += ey*Ez[EI(i3[2], j3[2], tz)]*V[cid[tx][ty][tz]];
      }
   }
   return s;
}

/* ------------------------------------------------------------------------
 * Gradient of E = sum_{mu nu in the pair} D_{mu nu} sum_K c_K (mu nu|K)
 * (D as for the b-pass, H = A c as from mdj_jprep):
 *   orbital centres: d/dA_x G_i = 2a G_{i+1} - i G_{i-1}, contracted with the
 *     pair's Hermite potential V_t = sum_aux lambda H R_{t+t'} up to L+1
 *   aux centres: dE/dC_x = -lambda sum_{t,t'} Dh_t H_t' R_{t+t'+e_x}, i.e.
 *     the lane-resolved potential at the shifted Hermite index
 * gA, gB (3 each) are set; gC (3 per atom) is added into.
 * ---------------------------------------------------------------------- */
void mdj_grad_pair(const int *ish_, const int *jsh_, const double *D, const double *H,
                   const double *schw, const double *auxb, const double *cut0,
                   double *gA, double *gB, double *gC)
{
   pair_t P; pair_setup(*ish_, *jsh_, &P);
   int li = P.li, lj = P.lj, L = li + lj, L1 = L + 1, nci = NCART(li), ncj = NCART(lj);
   int nhL = NHERM(L), nhL1 = NHERM(L1);
   double Dc[NCART(LMAXO)*NCART(LMAXO)], Ex[EISZ], Ey[EISZ], Ez[EISZ], Dh[NHERM(LMAXP)], V[NHERM(LMAXP + 1)];
   double B0[NHERM(NMAX)*NB], B1[NHERM(NMAX)*NB], F[(NMAX + 1)*NB];
   double Vl[NHERM(LMAXP + 1)*NB], W[NHERM(LMAXA)*NB];
   double al[NB], X[NB], Y[NB], Z[NB], T[NB], lam[NB];
   int shx[3][NHERM(LMAXP)];
   for (int d = 0; d < 3; d++) gA[d] = gB[d] = 0.0;
   to_cart(li, lj, D, nfun[li], Dc);
   double dmax = 0.0;
   for (int k = 0; k < nci*ncj; k++) dmax = fmax(dmax, fabs(Dc[k]));
   if (dmax == 0.0) return;
   for (int h = 0; h < nhL; h++) {
      int o = hoff[NMAX][h], t = o/(DD*DD), u = (o/DD)%DD, v = o%DD;
      shx[0][h] = cid[t+1][u][v]; shx[1][h] = cid[t][u+1][v]; shx[2][h] = cid[t][u][v+1];
   }
   double AB2 = (P.A[0]-P.B[0])*(P.A[0]-P.B[0]) + (P.A[1]-P.B[1])*(P.A[1]-P.B[1]) + (P.A[2]-P.B[2])*(P.A[2]-P.B[2]);
   int ic[NCART(LMAXO)][3], jc[NCART(LMAXO)][3];
   for (int a = 0; a < nci; a++) cart_xyz(li, a, &ic[a][0], &ic[a][1], &ic[a][2]);
   for (int b = 0; b < ncj; b++) cart_xyz(lj, b, &jc[b][0], &jc[b][1], &jc[b][2]);
   for (int pa = 0; pa < P.ni; pa++)
   for (int pb = 0; pb < P.nj; pb++) {
      double p, Kab, Pc[3];
      pair_hermite_setup(&P, pa, pb, AB2, &p, &Kab, Pc);
      if (fabs(Kab)*dmax*pow(M_PI/p, 1.5) < SCREEN_PRIM) continue;
      double a = P.ai[pa], b = P.aj[pb];
      herm1d(li + 1, lj + 1, p, Pc[0]-P.A[0], Pc[0]-P.B[0], Ex);
      herm1d(li + 1, lj + 1, p, Pc[1]-P.A[1], Pc[1]-P.B[1], Ey);
      herm1d(li + 1, lj + 1, p, Pc[2]-P.A[2], Pc[2]-P.B[2], Ez);
      for (int h = 0; h < nhL; h++) {
         int o = hoff[NMAX][h], t = o/(DD*DD), u = (o/DD)%DD, v = o%DD;
         double s = 0.0;
         for (int bb = 0; bb < ncj; bb++)
            for (int aa = 0; aa < nci; aa++) {
               if (t > ic[aa][0]+jc[bb][0] || u > ic[aa][1]+jc[bb][1] || v > ic[aa][2]+jc[bb][2]) continue;
               s += Dc[bb*nci + aa]*Ex[EI(ic[aa][0],jc[bb][0],t)]*Ey[EI(ic[aa][1],jc[bb][1],u)]*Ez[EI(ic[aa][2],jc[bb][2],v)];
            }
         Dh[h] = Kab*s;
      }
      memset(V, 0, sizeof(double)*nhL1);
      for (int lk = 0; lk <= lmaxa_used; lk++) {
         if (nq[lk] == 0) continue;
         int n = L1 + lk, nhk = NHERM(lk);
         const int *S = sumtab[L1][lk];
         for (int q0 = 0; q0 < nqp[lk]; q0 += NB) {
            lanes_setup(lk, q0, p, Pc, *schw, auxb, *cut0, NULL, 0.0, al, X, Y, Z, T, lam);
            int any = 0;
            for (int v = 0; v < NB; v++) any |= (lam[v] != 0.0);
            if (!any) continue;
            boys_lanes(n, T, F);
            double *R = rtab_lanes(n, al, X, Y, Z, F, B0, B1);
            for (int h2 = 0; h2 < nhk; h2++) {
               const double *Hp = H + gbase[lk] + (long)h2*nqp[lk] + q0;
               for (int v = 0; v < NB; v++) W[h2*NB + v] = lam[v]*Hp[v];
            }
            for (int h = 0; h < nhL1; h++) {
               double *vl = Vl + h*NB;
               for (int v = 0; v < NB; v++) vl[v] = 0.0;
               for (int h2 = 0; h2 < nhk; h2++) {
                  const double *r = R + S[h*nhk + h2]*NB, *w = W + h2*NB;
#pragma omp simd
                  for (int v = 0; v < NB; v++) vl[v] += w[v]*r[v];
               }
               double s = 0.0;
               for (int v = 0; v < NB; v++) s += vl[v];
               V[h] += s;
            }
            for (int d = 0; d < 3; d++) {
               double u[NB];
               for (int v = 0; v < NB; v++) u[v] = 0.0;
               for (int h = 0; h < nhL; h++) {
                  const double *vl = Vl + shx[d][h]*NB; double dh = Dh[h];
#pragma omp simd
                  for (int v = 0; v < NB; v++) u[v] += dh*vl[v];
               }
               for (int v = 0; v < NB; v++) {
                  int at = q_atom[lk][q0 + v];
                  if (at >= 0) gC[3*at + d] -= u[v];
               }
            }
         }
      }
      /* orbital centres */
      for (int bb = 0; bb < ncj; bb++)
         for (int aa = 0; aa < nci; aa++) {
            double dcv = Dc[bb*nci + aa]*Kab;
            if (dcv == 0.0) continue;
            for (int d = 0; d < 3; d++) {
               int i3[3] = { ic[aa][0], ic[aa][1], ic[aa][2] }, j3[3] = { jc[bb][0], jc[bb][1], jc[bb][2] };
               double tA, tB;
               i3[d]++;  tA = 2.0*a*ev3(i3, j3, Ex, Ey, Ez, V); i3[d]--;
               if (i3[d] > 0) { i3[d]--; tA -= (i3[d] + 1)*ev3(i3, j3, Ex, Ey, Ez, V); i3[d]++; }
               j3[d]++;  tB = 2.0*b*ev3(i3, j3, Ex, Ey, Ez, V); j3[d]--;
               if (j3[d] > 0) { j3[d]--; tB -= (j3[d] + 1)*ev3(i3, j3, Ex, Ey, Ez, V); j3[d]++; }
               gA[d] += dcv*tA; gB[d] += dcv*tB;
            }
         }
   }
}

/* ------------------------------------------------------------------------
 * Hessian kernels. Everything is built from the same Hermite potential
 * V^K_t = sum_{aux primitives on atom K} lambda H R_{t+t'} of a primitive
 * pair; the pair's side is contracted with derivative operators on the two
 * orbital centres (D_A,x f(i) = 2a f(i+e_x) - i_x f(i-e_x), and B alike),
 * the aux centre follows from translational invariance of every triple
 * integral, d_C = -(d_A + d_B), separately for each aux atom K.
 * ---------------------------------------------------------------------- */
typedef struct { double c; int di[3], dj[3]; } term_t;
typedef struct { int n, cen[2], dir[2]; } op_t;

static int apply_D(term_t *T, int n, int cen, int d, const int *i3, const int *j3, double a, double b)
{
   term_t out[16]; int m = 0;
   double ex = cen == 0 ? a : b;
   for (int k = 0; k < n; k++) {
      term_t p = T[k];
      p.c = T[k].c*2.0*ex;
      if (cen == 0) p.di[d]++; else p.dj[d]++;
      out[m++] = p;
      int cur = cen == 0 ? i3[d] + T[k].di[d] : j3[d] + T[k].dj[d];
      if (cur > 0) {
         term_t q = T[k];
         q.c = -T[k].c*cur;
         if (cen == 0) q.di[d]--; else q.dj[d]--;
         out[m++] = q;
      }
   }
   memcpy(T, out, sizeof(term_t)*m);
   return m;
}

static int op_terms(term_t *T, const op_t *op, const int *i3, const int *j3, double a, double b)
{
   int n = 1;
   memset(T, 0, sizeof(term_t));
   T[0].c = 1.0;
   for (int s = 0; s < op->n; s++) n = apply_D(T, n, op->cen[s], op->dir[s], i3, j3, a, b);
   return n;
}

static int htu[NHERM(NMAX)][3];
static int htu_ready = 0;
static void htu_init(void)
{
   for (int h = 0; h < NHERM(NMAX); h++) {
      int o = hoff[NMAX][h];
      htu[h][0] = o/(DD*DD); htu[h][1] = (o/DD)%DD; htu[h][2] = o%DD;
   }
   htu_ready = 1;
}

/* row[h] += w * sum_terms c * Ex Ey Ez at the hermite indices htu[h] */
static void add_terms(double *row, int nh, const term_t *T, int n, const int *i3, const int *j3, double w,
                      const double *Ex, const double *Ey, const double *Ez)
{
   for (int k = 0; k < n; k++) {
      int ix = i3[0] + T[k].di[0], iy = i3[1] + T[k].di[1], iz = i3[2] + T[k].di[2];
      int jx = j3[0] + T[k].dj[0], jy = j3[1] + T[k].dj[1], jz = j3[2] + T[k].dj[2];
      if (ix < 0 || iy < 0 || iz < 0 || jx < 0 || jy < 0 || jz < 0) continue;
      double wc = w*T[k].c;
      int sx = ix + jx, sy = iy + jy, sz = iz + jz;
      for (int h = 0; h < nh; h++) {
         int t = htu[h][0], u = htu[h][1], v = htu[h][2];
         if (t > sx || u > sy || v > sz) continue;
         row[h] += wc*Ex[EI(ix,jx,t)]*Ey[EI(iy,jy,u)]*Ez[EI(iz,jz,v)];
      }
   }
}

/* per-thread potentials by aux atom */
static __thread double *tl_V = NULL;
static __thread size_t tl_Vcap = 0;
static __thread char *tl_touched = NULL;
static __thread int *tl_list = NULL;
static __thread int tl_nat = 0;

typedef struct { double *V; int nh, nlist, nat; char *touched; int *list; } vat_t;

static void vat_prepare(vat_t *vt, int nat, int nh)
{
   size_t need = (size_t)nat*nh;
   if (need > tl_Vcap) {
      free(tl_V);
      tl_V = calloc(need, sizeof(double));
      tl_Vcap = need;
      if (!tl_V) abort();
      /* a fresh buffer is all zero: forget the old touch flags */
      if (tl_touched) memset(tl_touched, 0, tl_nat);
   }
   if (nat > tl_nat) {
      free(tl_touched); free(tl_list);
      tl_touched = calloc(nat, 1); tl_list = malloc(sizeof(int)*nat);
      tl_nat = nat;
   }
   vt->V = tl_V; vt->nh = nh; vt->nlist = 0; vt->nat = nat; vt->touched = tl_touched; vt->list = tl_list;
}

static void vat_clear(vat_t *vt)
{
   for (int k = 0; k < vt->nlist; k++) {
      int a = vt->list[k];
      memset(vt->V + (size_t)a*vt->nh, 0, sizeof(double)*vt->nh);
      vt->touched[a] = 0;
   }
   vt->nlist = 0;
}

/* V^K_t (t of order <= Lt) for a primitive pair of exponent p at Pc, all aux atoms */
static void vat_accumulate(vat_t *vt, int Lt, const double *H, double p, const double *Pc, double schw,
                           const double *auxb, double cut0, const double *cb, double cut)
{
   double B0[NHERM(NMAX)*NB], B1[NHERM(NMAX)*NB], F[(NMAX + 1)*NB], W[NHERM(LMAXA)*NB];
   double al[NB], X[NB], Y[NB], Z[NB], T[NB], lam[NB];
   int nh = NHERM(Lt);
   double Vl[NHERM(LMAXP + 2)*NB];
   for (int lk = 0; lk <= lmaxa_used; lk++) {
      if (nq[lk] == 0) continue;
      int n = Lt + lk, nhk = NHERM(lk);
      const int *S = sumtab[Lt][lk];
      for (int q0 = 0; q0 < nqp[lk]; q0 += NB) {
         lanes_setup(lk, q0, p, Pc, schw, cb ? cb : auxb, cb ? cut : cut0, cb ? auxb : NULL, cut0, al, X, Y, Z, T, lam);
         int any = 0;
         for (int v = 0; v < NB; v++) any |= (lam[v] != 0.0);
         if (!any) continue;
         boys_lanes(n, T, F);
         double *R = rtab_lanes(n, al, X, Y, Z, F, B0, B1);
         for (int h2 = 0; h2 < nhk; h2++) {
            const double *Hp = H + gbase[lk] + (long)h2*nqp[lk] + q0;
            for (int v = 0; v < NB; v++) W[h2*NB + v] = lam[v]*Hp[v];
         }
         for (int h = 0; h < nh; h++) {
            double acc[NB];
            for (int v = 0; v < NB; v++) acc[v] = 0.0;
            for (int h2 = 0; h2 < nhk; h2++) {
               const double *r = R + S[h*nhk + h2]*NB, *w = W + h2*NB;
#pragma omp simd
               for (int v = 0; v < NB; v++) acc[v] += w[v]*r[v];
            }
            for (int v = 0; v < NB; v++) Vl[h*NB + v] = acc[v];
         }
         for (int v = 0; v < NB; v++) {
            if (lam[v] == 0.0) continue;
            int at = q_atom[lk][q0 + v];
            if (at < 0) continue;
            if (!vt->touched[at]) { vt->touched[at] = 1; vt->list[vt->nlist++] = at; }
            double *Va = vt->V + (size_t)at*nh;
            for (int h = 0; h < nh; h++) Va[h] += Vl[h*NB + v];
         }
      }
   }
}

/* ---- pass 1: b and its first derivatives, per aux hermite slot ----
 * G0 += b; GI(:,x) += -d_A b, GJ(:,x) += -d_B b, GK(:,x) += d_A b + d_B b
 * (libcint's ip convention: bx of the aux atom is -d_C b = d_A b + d_B b). */
void mdj_b1_pair(const int *ish_, const int *jsh_, const double *D, const double *schw,
                 const double *auxb, const double *cut, const long *ldG, double *G0, double *GI, double *GJ, double *GK)
{
   if (!htu_ready) htu_init();
   pair_t P; pair_setup(*ish_, *jsh_, &P);
   int li = P.li, lj = P.lj, L = li + lj, nci = NCART(li), ncj = NCART(lj);
   int nh1 = NHERM(L + 1);
   double Dc[NCART(LMAXO)*NCART(LMAXO)], Ex[EISZ], Ey[EISZ], Ez[EISZ];
   double Dh[7][NHERM(LMAXP + 1)];
   double B0[NHERM(NMAX)*NB], B1[NHERM(NMAX)*NB], F[(NMAX + 1)*NB];
   double al[NB], X[NB], Y[NB], Z[NB], T[NB], lam[NB];
   static const op_t ops[7] = {
      {0,{0,0},{0,0}}, {1,{0,0},{0,0}}, {1,{0,0},{1,0}}, {1,{0,0},{2,0}},
      {1,{1,0},{0,0}}, {1,{1,0},{1,0}}, {1,{1,0},{2,0}} };
   to_cart(li, lj, D, nfun[li], Dc);
   double dmax = 0.0;
   for (int k = 0; k < nci*ncj; k++) dmax = fmax(dmax, fabs(Dc[k]));
   if (dmax == 0.0) return;
   double AB2 = (P.A[0]-P.B[0])*(P.A[0]-P.B[0]) + (P.A[1]-P.B[1])*(P.A[1]-P.B[1]) + (P.A[2]-P.B[2])*(P.A[2]-P.B[2]);
   int ic[NCART(LMAXO)][3], jc[NCART(LMAXO)][3];
   for (int a = 0; a < nci; a++) cart_xyz(li, a, &ic[a][0], &ic[a][1], &ic[a][2]);
   for (int b = 0; b < ncj; b++) cart_xyz(lj, b, &jc[b][0], &jc[b][1], &jc[b][2]);
   for (int pa = 0; pa < P.ni; pa++)
   for (int pb = 0; pb < P.nj; pb++) {
      double p, Kab, Pc[3];
      pair_hermite_setup(&P, pa, pb, AB2, &p, &Kab, Pc);
      if (fabs(Kab)*dmax*pow(M_PI/p, 1.5) < SCREEN_PRIM) continue;
      double a = P.ai[pa], b = P.aj[pb];
      herm1d(li + 1, lj + 1, p, Pc[0]-P.A[0], Pc[0]-P.B[0], Ex);
      herm1d(li + 1, lj + 1, p, Pc[1]-P.A[1], Pc[1]-P.B[1], Ey);
      herm1d(li + 1, lj + 1, p, Pc[2]-P.A[2], Pc[2]-P.B[2], Ez);
      memset(Dh, 0, sizeof(Dh));
      for (int bb = 0; bb < ncj; bb++)
         for (int aa = 0; aa < nci; aa++) {
            double w = Dc[bb*nci + aa]*Kab;
            if (w == 0.0) continue;
            for (int c = 0; c < 7; c++) {
               term_t T4[4];
               int nt = op_terms(T4, &ops[c], ic[aa], jc[bb], a, b);
               add_terms(Dh[c], nh1, T4, nt, ic[aa], jc[bb], w, Ex, Ey, Ez);
            }
         }
      for (int lk = 0; lk <= lmaxa_used; lk++) {
         if (nq[lk] == 0) continue;
         int n = L + 1 + lk, nhk = NHERM(lk);
         const int *S = sumtab[L + 1][lk];
         for (int q0 = 0; q0 < nqp[lk]; q0 += NB) {
            lanes_setup(lk, q0, p, Pc, *schw, auxb, *cut, NULL, 0.0, al, X, Y, Z, T, lam);
            int any = 0;
            for (int v = 0; v < NB; v++) any |= (lam[v] != 0.0);
            if (!any) continue;
            boys_lanes(n, T, F);
            double *R = rtab_lanes(n, al, X, Y, Z, F, B0, B1);
            for (int h2 = 0; h2 < nhk; h2++) {
               double acc[7][NB];
               for (int c = 0; c < 7; c++) for (int v = 0; v < NB; v++) acc[c][v] = 0.0;
               for (int h = 0; h < nh1; h++) {
                  const double *r = R + S[h*nhk + h2]*NB;
                  for (int c = 0; c < 7; c++) {
                     double d = Dh[c][h];
                     if (d == 0.0) continue;
#pragma omp simd
                     for (int v = 0; v < NB; v++) acc[c][v] += d*r[v];
                  }
               }
               long slot = gbase[lk] + (long)h2*nqp[lk] + q0;
               for (int v = 0; v < NB; v++) {
                  double l = lam[v];
                  G0[slot + v] += l*acc[0][v];
                  for (int x = 0; x < 3; x++) {
                     double aA = l*acc[1 + x][v], aB = l*acc[4 + x][v];
                     GI[x*(*ldG) + slot + v] -= aA;
                     GJ[x*(*ldG) + slot + v] -= aB;
                     GK[x*(*ldG) + slot + v] += aA + aB;
                  }
               }
            }
         }
      }
   }
}

/* ---- pass 3: second derivatives of F = sum D (mn|P) c_P, accumulated into
 * d2(3,nat,3,nat) (column-major, atoms 0-based): the blocks AA, AB, BB of the
 * orbital centres and, per aux atom K, the aux blocks by translational invariance. */
void mdj_hess_pair(const int *ish_, const int *jsh_, const double *D, const double *H, const double *schw,
                   const double *auxb, const double *cut, const int *nat_, const int *atomI_, const int *atomJ_,
                   double *d2)
{
   if (!htu_ready) htu_init();
   pair_t P; pair_setup(*ish_, *jsh_, &P);
   int nat = *nat_, aI = *atomI_, aJ = *atomJ_;
   int li = P.li, lj = P.lj, L = li + lj, nci = NCART(li), ncj = NCART(lj);
   int nh2 = NHERM(L + 2);
   double Dc[NCART(LMAXO)*NCART(LMAXO)], Ex[EISZ], Ey[EISZ], Ez[EISZ];
   double Wm[21][NHERM(LMAXP + 2)];
   static op_t ops[21];
   static int ops_ready = 0;
   if (!ops_ready) {
      int c = 0;
      for (int x = 0; x < 3; x++) for (int y = x; y < 3; y++) { ops[c].n = 2; ops[c].cen[0] = 0; ops[c].cen[1] = 0; ops[c].dir[0] = x; ops[c].dir[1] = y; c++; }
      for (int x = 0; x < 3; x++) for (int y = 0; y < 3; y++) { ops[c].n = 2; ops[c].cen[0] = 0; ops[c].cen[1] = 1; ops[c].dir[0] = x; ops[c].dir[1] = y; c++; }
      for (int x = 0; x < 3; x++) for (int y = x; y < 3; y++) { ops[c].n = 2; ops[c].cen[0] = 1; ops[c].cen[1] = 1; ops[c].dir[0] = x; ops[c].dir[1] = y; c++; }
      ops_ready = 1;
   }
   to_cart(li, lj, D, nfun[li], Dc);
   double dmax = 0.0;
   for (int k = 0; k < nci*ncj; k++) dmax = fmax(dmax, fabs(Dc[k]));
   if (dmax == 0.0) return;
   double AB2 = (P.A[0]-P.B[0])*(P.A[0]-P.B[0]) + (P.A[1]-P.B[1])*(P.A[1]-P.B[1]) + (P.A[2]-P.B[2])*(P.A[2]-P.B[2]);
   int ic[NCART(LMAXO)][3], jc[NCART(LMAXO)][3];
   for (int a = 0; a < nci; a++) cart_xyz(li, a, &ic[a][0], &ic[a][1], &ic[a][2]);
   for (int b = 0; b < ncj; b++) cart_xyz(lj, b, &jc[b][0], &jc[b][1], &jc[b][2]);
   vat_t vt; vat_prepare(&vt, nat, nh2);
   for (int pa = 0; pa < P.ni; pa++)
   for (int pb = 0; pb < P.nj; pb++) {
      double p, Kab, Pc[3];
      pair_hermite_setup(&P, pa, pb, AB2, &p, &Kab, Pc);
      if (fabs(Kab)*dmax*pow(M_PI/p, 1.5) < SCREEN_PRIM) continue;
      double a = P.ai[pa], b = P.aj[pb];
      herm1d(li + 2, lj + 2, p, Pc[0]-P.A[0], Pc[0]-P.B[0], Ex);
      herm1d(li + 2, lj + 2, p, Pc[1]-P.A[1], Pc[1]-P.B[1], Ey);
      herm1d(li + 2, lj + 2, p, Pc[2]-P.A[2], Pc[2]-P.B[2], Ez);
      memset(Wm, 0, sizeof(Wm));
      for (int bb = 0; bb < ncj; bb++)
         for (int aa = 0; aa < nci; aa++) {
            double w = Dc[bb*nci + aa]*Kab;
            if (w == 0.0) continue;
            for (int c = 0; c < 21; c++) {
               term_t T4[4];
               int nt = op_terms(T4, &ops[c], ic[aa], jc[bb], a, b);
               add_terms(Wm[c], nh2, T4, nt, ic[aa], jc[bb], w, Ex, Ey, Ez);
            }
         }
      vat_accumulate(&vt, L + 2, H, p, Pc, *schw, auxb, *cut, NULL, 0.0);
      for (int k = 0; k < vt.nlist; k++) {
         int K = vt.list[k];
         const double *Vk = vt.V + (size_t)K*nh2;
         double val[21];
         for (int c = 0; c < 21; c++) {
            double s = 0.0;
            for (int h = 0; h < nh2; h++) s += Wm[c][h]*Vk[h];
            val[c] = s;
         }
         double AA[3][3], AB[3][3], BB[3][3];
         int c = 0;
         for (int x = 0; x < 3; x++) for (int y = x; y < 3; y++) { AA[x][y] = AA[y][x] = val[c++]; }
         for (int x = 0; x < 3; x++) for (int y = 0; y < 3; y++) { AB[x][y] = val[c++]; }
         for (int x = 0; x < 3; x++) for (int y = x; y < 3; y++) { BB[x][y] = BB[y][x] = val[c++]; }
#define D2(x,a,y,b) d2[(x) + 3*(a) + 3*(long)nat*(y) + 9*(long)nat*(long)(b)]
         for (int x = 0; x < 3; x++)
            for (int y = 0; y < 3; y++) {
               D2(x,aI,y,aI) += AA[x][y];
               D2(x,aJ,y,aJ) += BB[x][y];
               D2(x,aI,y,aJ) += AB[x][y];
               D2(x,aJ,y,aI) += AB[y][x];
               double ik = -(AA[x][y] + AB[x][y]);          /* (I,K): x on A, y on C */
               D2(x,aI,y,K) += ik;  D2(y,K,x,aI) += ik;
               double jk = -(AB[y][x] + BB[x][y]);          /* (J,K): x on B, y on C */
               D2(x,aJ,y,K) += jk;  D2(y,K,x,aJ) += jk;
               D2(x,K,y,K) += AA[x][y] + AB[x][y] + AB[y][x] + BB[x][y];
            }
#undef D2
      }
      vat_clear(&vt);
   }
}

/* ---- Jx: first derivatives of the fitted Coulomb matrix block, F_mn = sum_P (mn|P) c_P ----
 * SA = d_A F, SB = d_B F (spherical di x dj x 3, no AO norms), and for each aux atom K
 * touched SK = (d_A + d_B) F^K = -d_C F^K, K listed in Klist (0-based), SK(:,:,:,k). */
static __thread double *tl_E6 = NULL;
static __thread size_t tl_E6cap = 0;
static __thread double *tl_Kacc = NULL;
static __thread size_t tl_Kcap = 0;
static __thread char *tl_kt = NULL;
static __thread int *tl_klist = NULL;
static __thread int tl_knat = 0;

static void cart2sph_blk(int li, int lj, const double *Jc, double *Js)
{
   int nci = NCART(li), ncj = NCART(lj), di = nfun[li], dj = nfun[lj];
   double tmp[(2*LMAXO+1)*NCART(LMAXO)];
   for (int bb = 0; bb < ncj; bb++)
      for (int m = 0; m < di; m++) {
         double s = 0.0;
         for (int aa = 0; aa < nci; aa++) s += c2s[li][m*nci + aa]*Jc[bb*nci + aa];
         tmp[bb*di + m] = s;
      }
   for (int n2 = 0; n2 < dj; n2++)
      for (int m = 0; m < di; m++) {
         double s = 0.0;
         for (int bb = 0; bb < ncj; bb++) s += tmp[bb*di + m]*c2s[lj][n2*ncj + bb];
         Js[n2*di + m] = s;
      }
}

void mdj_jx_pair(const int *ish_, const int *jsh_, const double *H, const double *schw, const double *auxb,
                 const double *cut, const int *nat_, double *SA, double *SB, int *nK, int *Klist, double *SK)
{
   if (!htu_ready) htu_init();
   pair_t P; pair_setup(*ish_, *jsh_, &P);
   int nat = *nat_;
   int li = P.li, lj = P.lj, L = li + lj, nci = NCART(li), ncj = NCART(lj), nab = nci*ncj, di = nfun[li], dj = nfun[lj];
   int nh1 = NHERM(L + 1);
   double Ex[EISZ], Ey[EISZ], Ez[EISZ];
   static const op_t ops[6] = { {1,{0,0},{0,0}}, {1,{0,0},{1,0}}, {1,{0,0},{2,0}},
                                {1,{1,0},{0,0}}, {1,{1,0},{1,0}}, {1,{1,0},{2,0}} };
   size_t need = (size_t)9*nab*nh1;     /* A_x, B_x, (A+B)_x */
   if (need > tl_E6cap) { free(tl_E6); tl_E6 = malloc(need*sizeof(double)); tl_E6cap = need; if (!tl_E6) abort(); }
   if ((size_t)nat*3*nab > tl_Kcap) {
      free(tl_Kacc); tl_Kcap = (size_t)nat*3*225;
      if ((size_t)nat*3*nab > tl_Kcap) tl_Kcap = (size_t)nat*3*nab;
      tl_Kacc = calloc(tl_Kcap, sizeof(double)); if (!tl_Kacc) abort();
      if (tl_kt) memset(tl_kt, 0, tl_knat);
   }
   if (nat > tl_knat) { free(tl_kt); free(tl_klist); tl_kt = calloc(nat, 1); tl_klist = malloc(sizeof(int)*nat); tl_knat = nat; }
   double *EA = tl_E6, *EB = tl_E6 + (size_t)3*nab*nh1, *ES = tl_E6 + (size_t)6*nab*nh1;
   double Acc[3][NCART(LMAXO)*NCART(LMAXO)], Bcc[3][NCART(LMAXO)*NCART(LMAXO)];
   memset(Acc, 0, sizeof(Acc)); memset(Bcc, 0, sizeof(Bcc));
   int nk = 0;
   double AB2 = (P.A[0]-P.B[0])*(P.A[0]-P.B[0]) + (P.A[1]-P.B[1])*(P.A[1]-P.B[1]) + (P.A[2]-P.B[2])*(P.A[2]-P.B[2]);
   int ic[NCART(LMAXO)][3], jc[NCART(LMAXO)][3];
   for (int a = 0; a < nci; a++) cart_xyz(li, a, &ic[a][0], &ic[a][1], &ic[a][2]);
   for (int b = 0; b < ncj; b++) cart_xyz(lj, b, &jc[b][0], &jc[b][1], &jc[b][2]);
   vat_t vt; vat_prepare(&vt, nat, nh1);
   double Vtot[NHERM(LMAXP + 1)];
   for (int pa = 0; pa < P.ni; pa++)
   for (int pb = 0; pb < P.nj; pb++) {
      double p, Kab, Pc[3];
      pair_hermite_setup(&P, pa, pb, AB2, &p, &Kab, Pc);
      if (fabs(Kab)*pow(M_PI/p, 1.5) < SCREEN_PRIM) continue;
      double a = P.ai[pa], b = P.aj[pb];
      herm1d(li + 1, lj + 1, p, Pc[0]-P.A[0], Pc[0]-P.B[0], Ex);
      herm1d(li + 1, lj + 1, p, Pc[1]-P.A[1], Pc[1]-P.B[1], Ey);
      herm1d(li + 1, lj + 1, p, Pc[2]-P.A[2], Pc[2]-P.B[2], Ez);
      memset(tl_E6, 0, need*sizeof(double));
      for (int bb = 0; bb < ncj; bb++)
         for (int aa = 0; aa < nci; aa++) {
            int ab = bb*nci + aa;
            for (int c = 0; c < 6; c++) {
               term_t T4[4];
               int nt = op_terms(T4, &ops[c], ic[aa], jc[bb], a, b);
               double *row = (c < 3 ? EA + ((size_t)c*nab + ab)*nh1 : EB + ((size_t)(c-3)*nab + ab)*nh1);
               add_terms(row, nh1, T4, nt, ic[aa], jc[bb], Kab, Ex, Ey, Ez);
            }
         }
      for (size_t k = 0; k < (size_t)3*nab*nh1; k++) ES[k] = EA[k] + EB[k];
      vat_accumulate(&vt, L + 1, H, p, Pc, *schw, auxb, *cut, NULL, 0.0);
      for (int h = 0; h < nh1; h++) Vtot[h] = 0.0;
      for (int k = 0; k < vt.nlist; k++) {
         int K = vt.list[k];
         const double *Vk = vt.V + (size_t)K*nh1;
         for (int h = 0; h < nh1; h++) Vtot[h] += Vk[h];
         if (!tl_kt[K]) { tl_kt[K] = 1; tl_klist[nk++] = K; memset(tl_Kacc + (size_t)K*3*nab, 0, sizeof(double)*3*nab); }
         for (int x = 0; x < 3; x++)
            for (int ab = 0; ab < nab; ab++) {
               const double *e = ES + ((size_t)x*nab + ab)*nh1;
               double s = 0.0;
               for (int h = 0; h < nh1; h++) s += e[h]*Vk[h];
               tl_Kacc[((size_t)K*3 + x)*nab + ab] += s;
            }
      }
      for (int x = 0; x < 3; x++)
         for (int ab = 0; ab < nab; ab++) {
            const double *ea = EA + ((size_t)x*nab + ab)*nh1, *eb = EB + ((size_t)x*nab + ab)*nh1;
            double sa = 0.0, sb = 0.0;
            for (int h = 0; h < nh1; h++) { sa += ea[h]*Vtot[h]; sb += eb[h]*Vtot[h]; }
            Acc[x][ab] += sa; Bcc[x][ab] += sb;
         }
      vat_clear(&vt);
   }
   for (int x = 0; x < 3; x++) {
      cart2sph_blk(li, lj, Acc[x], SA + (size_t)x*di*dj);
      cart2sph_blk(li, lj, Bcc[x], SB + (size_t)x*di*dj);
   }
   for (int k = 0; k < nk; k++) {
      int K = tl_klist[k];
      Klist[k] = K;
      for (int x = 0; x < 3; x++)
         cart2sph_blk(li, lj, tl_Kacc + ((size_t)K*3 + x)*nab, SK + ((size_t)k*3 + x)*di*dj);
      tl_kt[K] = 0;
   }
   *nK = nk;
}

/* ------------------------------------------------------------------------
 * Raw 3-centre block of one shell pair: B(K + nauxf*(m + di*n)) = (mu_m nu_n | K), aux function fastest, zero for the aux
 * shells cut by schw*auxb[k] < cut; spherical (or Cartesian) orbital functions, the same unnormalised convention as the
 * passes above (the caller multiplies by NorVEC / NorVECAux). The J-pass with unit aux coefficients: per primitive pair the
 * Hermite potential of every aux function, V(h, f) = sum_g lam A_g(f, h2) R(h + h2), contracted with the pair's Hermite
 * expansion E(ab, h). wk: nci*ncj*nauxf doubles of scratch (Cartesian accumulator).
 * ---------------------------------------------------------------------- */
/* Bc(ab, f) += sum_h E(ab, h) V(h, f), f over the functions of one aux shell */
static inline void b3c_flush(int nab, int nhL, int ns, const double *Em, const double *V, double *Bc)
{
   for (int f = 0; f < ns; f++) {
      double *Bk = Bc + (long)f*nab;
      const double *Vf = V + f*nhL;
      for (int h = 0; h < nhL; h++) {
         double vh = Vf[h];
         if (vh == 0.0) continue;
         const double *Eh = Em + h*nab;
#pragma omp simd
         for (int ab = 0; ab < nab; ab++) Bk[ab] += Eh[ab]*vh;
      }
   }
}

long mdj_nauxf(void) { return naux > 0 ? (long)ax_foff[naux-1] + nfun[ax_l[naux-1]] : 0; }

void mdj_b3c_pair(const int *ish_, const int *jsh_, const double *schw, const double *auxb, const double *cut,
                  double *B, double *wk)
{
   pair_t P; pair_setup(*ish_, *jsh_, &P);
   const long nf = mdj_nauxf();
   int li = P.li, lj = P.lj, L = li + lj, nci = NCART(li), ncj = NCART(lj), nhL = NHERM(L), nab = nci*ncj;
   int di = nfun[li], dj = nfun[lj];
   double Ex[EISZ], Ey[EISZ], Ez[EISZ];
   double B0[NHERM(NMAX)*NB], B1[NHERM(NMAX)*NB], F[(NMAX + 1)*NB];
   double al[NB], X[NB], Y[NB], Z[NB], T[NB], lam[NB];
   double Em[NCART(LMAXO)*NCART(LMAXO)*NHERM(LMAXP)];        /* E(ab, h), ab fastest */
   double V[NHERM(LMAXP)*(2*LMAXA+1+16)];                     /* V(h, f) of one aux shell, h fastest */
   double *Bc = wk;                                           /* Bc(ab, K) */
   for (int k = 0; k < naux; k++)
      if (*schw*auxb[k] >= *cut) memset(Bc + (long)ax_foff[k]*nab, 0, sizeof(double)*nab*nfun[ax_l[k]]);
   double AB2 = (P.A[0]-P.B[0])*(P.A[0]-P.B[0]) + (P.A[1]-P.B[1])*(P.A[1]-P.B[1]) + (P.A[2]-P.B[2])*(P.A[2]-P.B[2]);
   int ix[NCART(LMAXO)], iy[NCART(LMAXO)], iz[NCART(LMAXO)], jx[NCART(LMAXO)], jy[NCART(LMAXO)], jz[NCART(LMAXO)];
   for (int a = 0; a < nci; a++) cart_xyz(li, a, ix+a, iy+a, iz+a);
   for (int b = 0; b < ncj; b++) cart_xyz(lj, b, jx+b, jy+b, jz+b);
   for (int pa = 0; pa < P.ni; pa++)
   for (int pb = 0; pb < P.nj; pb++) {
      double a = P.ai[pa], b = P.aj[pb], p = a + b;
      double Kab = exp(-a*b/p*AB2)*P.ci[pa]*P.cj[pb];
      if (fabs(Kab)*pow(M_PI/p, 1.5) < SCREEN_PRIM) continue;
      double Pc[3];
      for (int d = 0; d < 3; d++) Pc[d] = (a*P.A[d] + b*P.B[d])/p;
      herm1d(li, lj, p, Pc[0]-P.A[0], Pc[0]-P.B[0], Ex);
      herm1d(li, lj, p, Pc[1]-P.A[1], Pc[1]-P.B[1], Ey);
      herm1d(li, lj, p, Pc[2]-P.A[2], Pc[2]-P.B[2], Ez);
      for (int h = 0; h < nhL; h++) {
         int o = hoff[L][h], t = o/(DD*DD), u = (o/DD)%DD, v = o%DD;
         for (int bb = 0; bb < ncj; bb++)
            for (int aa = 0; aa < nci; aa++) {
               double e = 0.0;
               if (t <= ix[aa]+jx[bb] && u <= iy[aa]+jy[bb] && v <= iz[aa]+jz[bb])
                  e = Kab*Ex[EI(ix[aa],jx[bb],t)]*Ey[EI(iy[aa],jy[bb],u)]*Ez[EI(iz[aa],jz[bb],v)];
               Em[h*nab + bb*nci + aa] = e;
            }
      }
      /* aux primitives in lane blocks of one l (ordered by shell, then primitive): V(h, f) of the current aux shell
       * accumulates over its primitives and goes into Bc when the shell changes */
      for (int lk = 0; lk <= lmaxa_used; lk++) {
         if (nq[lk] == 0) continue;
         int n = L + lk, nhk = NHERM(lk), ns = nfun[lk], cur = -1;
         const int *S = sumtab[L][lk];
         for (int q0 = 0; q0 < nqp[lk]; q0 += NB) {
            lanes_setup(lk, q0, p, Pc, *schw, auxb, *cut, NULL, 0.0, al, X, Y, Z, T, lam);
            int any = 0;
            for (int v = 0; v < NB; v++) any |= (lam[v] != 0.0);
            if (!any) continue;
            boys_lanes(n, T, F);
            double *R = rtab_lanes(n, al, X, Y, Z, F, B0, B1);
            for (int v = 0; v < NB; v++) {
               if (lam[v] == 0.0) continue;
               int k = q_shell[lk][q0 + v];
               if (k != cur) {
                  if (cur >= 0) b3c_flush(nab, nhL, ns, Em, V, Bc + (long)ax_foff[cur]*nab);
                  cur = k;
                  memset(V, 0, sizeof(double)*nhL*ns);
               }
               const double *A = ax_A + aoff_of(k, q_prim[lk][q0 + v]);    /* A(f, h2) row-major */
               for (int h2 = 0; h2 < nhk; h2++) {
                  double Rg[NHERM(LMAXP)];
                  for (int h = 0; h < nhL; h++) Rg[h] = lam[v]*R[S[h*nhk + h2]*NB + v];
                  for (int f = 0; f < ns; f++) {
                     double w = A[f*nhk + h2];
                     if (w == 0.0) continue;
                     double *Vf = V + f*nhL;
#pragma omp simd
                     for (int h = 0; h < nhL; h++) Vf[h] += w*Rg[h];
                  }
               }
            }
         }
         if (cur >= 0) b3c_flush(nab, nhL, ns, Em, V, Bc + (long)ax_foff[cur]*nab);
      }
   }
   /* B(K, :) = c2s_i Bc(:, K) c2s_j^T for the surviving aux shells, 0 for the others */
   double tmp[(2*LMAXO+1)*NCART(LMAXO)];
   for (int k = 0; k < naux; k++) {
      if (*schw*auxb[k] < *cut) {
         for (int f = 0; f < nfun[ax_l[k]]; f++)
            for (int e = 0; e < di*dj; e++) B[ax_foff[k] + f + nf*e] = 0.0;
         continue;
      }
      for (int f = 0; f < nfun[ax_l[k]]; f++) {
         long K = ax_foff[k] + f;
         const double *Jc = Bc + K*nab;
         for (int bb = 0; bb < ncj; bb++)
            for (int m = 0; m < di; m++) {
               double s = 0.0;
               for (int aa = 0; aa < nci; aa++) s += c2s[li][m*nci + aa]*Jc[bb*nci + aa];
               tmp[bb*di + m] = s;
            }
         for (int nn = 0; nn < dj; nn++)
            for (int m = 0; m < di; m++) {
               double s = 0.0;
               for (int bb = 0; bb < ncj; bb++) s += tmp[bb*di + m]*c2s[lj][nn*ncj + bb];
               B[K + nf*(long)(nn*di + m)] = s;
            }
      }
   }
}

/* ------------------------------------------------------------------------
 * Read access to the aux primitive tables of the vectorised layout, for a device back end: per l the
 * primitive count nq, the padded count nqp, the offset gbase of the l group in G / H (G[gbase + h*nqp + q]), the
 * exponents, centres and aux shell of every primitive (shell -1: padding). Returns the highest aux l in use, -1 when
 * the tables are not built.
 * ---------------------------------------------------------------------- */
int mdj_aux_tables(int l, int *nq_, int *nqp_, long *gbase_, const double **ex, const double **cx, const double **cy,
                   const double **cz, const int **sh)
{
   if (naux <= 0 || l < 0 || l > LMAXA) return -1;
   *nq_ = nq[l]; *nqp_ = nqp[l]; *gbase_ = gbase[l];
   *ex = q_exp[l]; *cx = q_cx[l]; *cy = q_cy[l]; *cz = q_cz[l]; *sh = q_shell[l];
   return lmaxa_used;
}
long mdj_gsize(void) { return gsize2; }
