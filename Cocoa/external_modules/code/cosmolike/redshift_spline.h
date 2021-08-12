#ifndef __COSMOLIKE_REDSHIFT_SPLINE_H
#define __COSMOLIKE_REDSHIFT_SPLINE_H
#ifdef __cplusplus
extern "C" {
#endif

// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// integration boundary routines
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------

double amin_source(const int ni);

double amax_source(const int ni);

double amin_lens(const int ni);

double amax_lens(const int ni);

// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// redshift overlap tests, allowed tomography combinations
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------

int test_kmax(double l, int zl); // test whether the (l,zl) bin is in the linear clustering regime

int test_zoverlap(int zl, int zs); // test whether source bin zs is behind lens bin zl

int test_zoverlap_c(int zc, int zs); // test whether source bin zs is behind lens bin zl (clusters)


int ZL(const int Nbin); // find zlens bin of tomo combination (galaxy-galaxy lensing)

int ZS(const int Nbin); // find zsource bin of tomo combination (galaxy-galaxy lensing)

int N_ggl(const int zl, const int zs); // find tomo bin number tomography combination


int ZCL(const int Nbin); // find zlens bin of tomo combination (cluster-galaxy lensing)

int ZCS(const int Nbin); // find zsource bin of tomo combination (cluster-galaxy lensing)

int N_cgl(const int zl, const int zs); // find tomo bin number tomography combination


int Z1(const int Nbin); // find z1 bin of tomo combination (cosmic shear)

int Z2(const int Nbin); // find z2 bin of tomo combination (cosmic shear)

// find tomo bin number tomography combination (cosmic shear)
int N_shear(const int z1, const int z2);


int ZCL1(const int Nbin); // find z1 bin of tomo combination (galaxy clustering)

int ZCL2(const int Nbin); // find z2 bin of tomo combination (galaxy clustering)

// find tomo bin number tomography combination (galaxy clustering)
int N_CL(const int z1, const int z2);


int ZCCL1(const int Nbin); // find z1 bin of tomo combination (cluster clustering)

int ZCCL2(const int Nbin); // find z2 bin of tomo combination (cluster clustering)

// find tomo bin number tomography combination (cluster clustering)
int N_CCL(const int z1, const int z2);


int ZCGCL1(const int ni); // find z1 bin of tomo combination (cluster-galaxy cross clustering)

int ZCGCL2(const int nj); // find z2 bin of tomo combination (cluster-galaxy cross clustering)

// find tomo bin number tomography combination (cluster-galaxy cross clustering)
int N_CGCL(const int ni, const int nj); // ni = Cluster Nbin, nj = Galaxy Nbin

// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// shear routines (for redshift distributions, including photo-zs * (optional))
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------

// returns n(ztrue | j), works only with binned distributions; 
double zdistr_photoz(double z, int j); // j =-1 -> no tomography; j>= 0 -> tomography bin j

double sigma_zphot_shear(double z, int nz); // compute photo-z scatter (sigma(z))

double bias_zphot_shear(double z, int nz); // compute photo-z bias

// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// clustering routines (for redshift distributions, including photo-zs * (optional))
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------

// returns n(ztrue | j), works only with binned distributions;
double pf_photoz(double zz, int j); // j =-1 -> no tomography; j>= 0 -> tomography bin j

double sigma_zphot_clustering(double z, int nz); // compute photo-z scatter (sigma(z))

double bias_zphot_clustering(double z, int nz); // compute photo-z bias

// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// cluster routines (for redshift distributions, including photo-zs * (optional))
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------

//return pf(z,j) based on redshift file with structure z[i] nz[0][i]..nz[tomo.clustering_Nbin-1][i]
double pf_cluster_histo_n(double z, void* params); 

double pz_cluster(double zz, int nz, int nN); // int_zmin^zmax dz_obs p(z_obs|z_true)

// simplfied selection function, disregards evolution of N-M relation+mass function within z bin
double zdistr_cluster(int nz,  double z, double chi, double hoverh0) ;

// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// routines calculating the number of source and lens galaxies per bin
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------

// returns n_gal for shear tomography bin j, works only with binned
double nsource(int j); // distributions; j =-1 -> no tomography; j>= 0 -> tomography bin j

// returns n_gal for clustering tomography bin j, works only with binned
double nlens(int j); // distributions; j =-1 -> no tomography; j>= 0 -> tomography bin j

double zmean(int j); // mean true redshift of (clustering/lens) galaxies within redshift bin j

double zmean_source(int j); // mean true redshift of source galaxies in tomography bin j

// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// lens efficiencies
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------

double g_cmb(const double a); // lens efficiency for CMB lensing

// lens efficiency of source galaxies in tomography bin ni
double g_tomo(const double a, const int ni); 

// lens efficiency of source galaxies in tomography bin zbin
double g2_tomo(const double a, const int ni) ;

// lens efficiency of *lens* galaxies in tomography bin zbin
// used in magnification calculations
double g_lens(const double a, const int ni);

// lens efficiency of lens cluster in tomo bin nz, lambda bin nl
// used in magnification calculations
double g_lens_cluster(const double a, const int nz, const int nl);

// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// routines for association of a pair redshift bin numbers and power spectrum tomography bin
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------

double ggl_efficiency(int zl, int zs);


#ifdef __cplusplus
}
#endif
#endif // HEADER GUARD