#ifndef __COSMOLIKE_COSMO2D_H
#define __COSMOLIKE_COSMO2D_H
#ifdef __cplusplus
extern "C" {
#endif

// ----------------------------------------------------------------------------
// Naming convention:
// ----------------------------------------------------------------------------
// c = cluster position ("c" as in "cluster")
// g = galaxy positions ("g" as in "galaxy")
// k = kappa CMB ("k" as in "kappa")
// s = kappa from source galaxies ("s" as in "shear")

// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// Correlation Functions (real Space) - Full Sky - bin average
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------

// ss in real space has a special name
double xi_pm_tomo(const int pm, const int nt, const int ni, const int nj, const int limber);

// gs in real space has a special name
double w_gammat_tomo(const int nt, const int ni, const int nj, const int limber);

double w_gg_tomo(const int nt, const int ni, const int nj, const int limber);

double w_gk_tomo(const int nt, const int ni, const int limber);

double w_ks_tomo(const int nt, const int ni, const int limber);

// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// Correlation Functions (real space) - flat sky
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------

// ss in real space has a special name
double xi_pm_tomo_flatsky(const int pm, double theta, const int ni, const int nj,
  const int limber);

// gs in real space has a special name
double w_gammat_tomo_flatsky(double theta, const int ni, const int nj, const int limber);

// WARNING: C_gg beyond linear bias for cross-tomography bins not yet supported
double w_gg_tomo_flatsky(double theta, const int ni, const int nj, const int limber);

double w_gk_tomo_flatsky(double theta, const int ni, const int limber);

double w_ks_tomo_flatsky(double theta, const int ni, const int limber);

// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// Non Limber (Angular Power Spectrum)
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------

void f_chi_for_Psi_sh_IA(double* chi, int Nchi, double* f_chi, const int nj, double zmin, double zmax);

void f_chi_for_Psi_sh(double* chi, int Nchi, double* f_chi, const int nj, double zmax);

void f_chi_for_Psi_cl(double* chi, int Nchi, double* f_chi, const int ni, double zmin, double zmax);

void f_chi_for_Psi_cl_RSD(double* chi, int Nchi, double* f_chi, const int ni, double zmin, double zmax);

void f_chi_for_Psi_cl_Mag(double* chi, int Nchi, double* f_chi, const int ni, double zmax);

void C_gl_tomo(int L, int nl, int ns, double* Cl, double dev, double tolerance);

void C_cl_tomo(int L, const int ni, const int nj, double* Cl, double dev, double tolerance);

// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// Limber Approximation (Angular Power Spectrum)
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------

double C_ss_tomo_TATT_EE_limber_nointerp(const double l, const int ni, const int nj,
const int init_static_vars_only); // IA=(5||6)

double C_ss_tomo_TATT_EE_limber(const double l, const int ni, const int nj); // IA=(5||6)

double C_ss_tomo_TATT_BB_limber_nointerp(const double l, const int ni, const int nj,
const int init_static_vars_only); // IA=(5||6)

double C_ss_tomo_TATT_BB_limber(const double l, const int ni, const int nj); // IA=(5||6)


double C_ss_tomo_limber_nointerp(const double l, const int ni, const int nj, const int use_linear_ps,
const int init_static_vars_only);

double C_ss_tomo_limber(const double l, const int ni, const int nj); // IA=(0||3||4)


// works with IA=(0||3||4||5||6)
double C_gs_tomo_limber_nointerp(const double l, const int ni, const int nj, const int use_linear_ps,
const int init_static_vars_only);

double C_gs_tomo_limber(const double l, const int ni, const int nj);

// ----------------------------------------------------------------------------
// All functions below can run w/ like.IA=0 || like.IA=3 || like.IA=4
// ----------------------------------------------------------------------------

double C_gg_tomo_limber_nointerp(const double l, const int ni, const int nj, const int use_linear_ps,
const int init_static_vars_only);

double C_gg_tomo_limber(const double l, const int ni, const int nj);


double C_gk_tomo_limber_nointerp(const double l, int nl, const int use_linear_ps,
const int init_static_vars_only);

double C_gk_tomo_limber(const double l, const int ni);


double C_ks_tomo_limber_nointerp(const double l, int ns, const int use_linear_ps,
const int init_static_vars_only);

double C_ks_tomo_limber(const double l, const int ni);


double C_kk_limber_nointerp(const double l, const int use_linear_ps, const int init_static_vars_only);

double C_kk_limber(const double l);

#ifdef __cplusplus
}
#endif
#endif // HEADER GUARD