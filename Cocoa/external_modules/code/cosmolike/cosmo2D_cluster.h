#ifndef __COSMOLIKE_COSMO2D_CLUSTERH
#define __COSMOLIKE_COSMO2D_CLUSTERH
#ifdef __cplusplus
extern "C" {
#endif

// ----------------------------------------------------------------------------
// Naming convention: (same as cosmo2D.h)
// ----------------------------------------------------------------------------
// c = cluster position ("c" as in "cluster")
// g = galaxy positions ("g" as in "galaxy")
// k = kappa CMB ("k" as in "kappa")
// s = kappa from source galaxies ("s" as in "shear")

// ----------------------------------------------------------------------------
// Threading
// ----------------------------------------------------------------------------
// Thread loops in lambda_obs is not allowed. Most functions update static arrays 
// when varying lambda_obs in cluster_utils. With lambda_obs fixed, loops on 
// redshift bins can be threaded using the standard loop unrolling technique

// -------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------
// Correlation Functions (real Space) - Full Sky - bin average
// -------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------
// nt = theta bin

// nl = lambda_obs bin, ni = cluster redshift bin, nj = source redshift bin
double w_gammat_cluster_tomo(int nt, int nl, int ni, int nj, int limber);

// nl{1,2} = lambda_obs bins, n{i,j} = cluster redshift bins
double w_cc_tomo(int nt, int nl1, int nl2, int ni, int nj, int limber);

// nl = lambda_obs bin, ni = cluster redshift bin, nj = galaxy redshift bin
double w_cg_tomo(int nt, int nl, int ni, int nj, int limber);

// -------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------
// Correlation Functions (real space) - flat sky
// -------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------

// nl = lambda_obs bin, ni = cluster redshift bin, nj = source redshift bin
double w_gammat_cluster_tomo_flatsky(double theta, int nl, int ni, int nj, int limber);

// nl{1,2} = lambda_obs bins, n{i,j} = cluster redshift bins
double w_cc_tomo_flatsky(double theta, int nl1, int nl2, int ni, int nj, int limber);

// nl = lambda_obs bin, ni = cluster redshift bin, nj = galaxy redshift bin
double w_cg_tomo_flatsky(double theta, int nl, int ni, int nj, int limber);

// -------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------
// Non Limber (Angular Power Spectrum)
// -------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------
void f_chi_for_Psi_cluster_cl(double* chi, int Nchi, double* fchi, int ni, int nl, double zmin,
double zmax);

void f_chi_for_Psi_cluster_cl_RSD(double* chi, int Nchi, double* fchi, int ni, int nl, double zmin, 
double zmax);

void f_chi_for_Psi_cluster_cl_Mag(double* chi, int Nchi, double* fchi, int ni, int nl, double zmax);

// nl = lambda_obs bin, ni = cluster redshift bin, nj = galaxy redshift bin
void C_cg_tomo(int L, int nl, int ni, int nj, double* Cl, double dev, double tol);

// nl{1,2} = lambda_obs bins, n{i,j} = cluster redshift bins
void C_cc_tomo(int L, int nl1, int nl2, int ni, int nj, double* Cl, double dev, double tol);

// -------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------
// Limber Approximation (Angular Power Spectrum)
// -------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------

// nl = lambda_obs bin, ni = cluster redshift bin, nj = source redshift bin
double C_cs_tomo_limber_nointerp(double l, int nl, int ni, int nj, int use_linear_ps);
double C_cs_tomo_limber(double l, int nl, int ni, int nj);

// nl{1,2} = lambda_obs bins, n{i,j} = cluster redshift bins
double C_cc_tomo_limber_nointerp(double l, int nl1, int nl2, int ni, int nj, int use_linear_ps);
double C_cc_tomo_limber(double l, int nl1, int nl2, int ni, int nj);

// nl = lambda_obs bin, ni = cluster redshift bin, nj = galaxy redshift bin
double C_cg_tomo_limber_nointerp(double l, int nl, int ni, int nj, int use_linear_ps);
double C_cg_tomo_limber(double l, int nl, int ni, int nj);

// -------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------
// cluster number counts
// -------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------

// nl = lambda_obs bin, ni = cluster redshift bin
double projected_average_number_counts(int nl, int nz);

#ifdef __cplusplus
}
#endif
#endif // HEADER GUARD