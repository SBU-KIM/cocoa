#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_spline2d.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_integration.h>

#include "basics.h"
#include "bias.h"
#include "cosmo3D.h"
#include "cluster_util.h"
#include "halo.h"
#include "recompute.h"
#include "radial_weights.h"
#include "redshift_spline.h"
#include "structs.h"

#include "log.c/src/log.h"

static int INTERPOLATE_SURVEY_AREA = 1;
static int ASSUME_CONSTANT_LAMBD_HALO_EXCLUSION = 1;

static int GSL_WORKSPACE_SIZE = 250;
static double M_PIVOT = 5E14; //Msun/h
static double log_M_min = 12.0;
static double log_M_max = 15.9; 
static double LOG10_E = 0.4342944819;

static int has_b2_galaxies()
{
  int res = 0;
  for(int i=0; i<tomo.clustering_Nbin; i++) 
  {
    if(gbias.b2[i])
    {
      res = 1;
    }
  }
  return res 
}

static int binned_P_lambda_obs_given_M_size_z_table = 10;
static int binned_P_lambda_obs_given_M_size_M_table = 50;
static double binned_P_lambda_obs_given_M_zmin_table = 0.20;
static double binned_P_lambda_obs_given_M_zmax_table = 0.80;
static int binned_p_cm_size_a_table = 30;

// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// BUZZARD
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------

double buzzard_P_lambda_obs_given_M_nointerp(double obs_lambda, void* params)
{      
  double* array = (double*) params;
  
  const double mass = array[0];
  const double cluster_z = array[1];
  const double lnlm0 = nuisance.cluster_MOR[0]; 
  const double Alm = nuisance.cluster_MOR[1];
  const double sigma_lnlm_intrinsic = nuisance.cluster_MOR[2];
  const double Blm = nuisance.cluster_MOR[3];
  const double Mpiv = M_PIVOT; //Msun/h
  const double lnlm = lnlm0 + (Alm)*log(mass/Mpiv) + Blm*log((1+cluster_z)/1.45);
  const double sigma_total = (lnlm>0) ? 
    sqrt(sigma_lnlm_intrinsic*sigma_lnlm_intrinsic+(exp(lnlm)-1.)/exp(2*lnlm)) : 
    sqrt(sigma_lnlm_intrinsic*sigma_lnlm_intrinsic);
  const double x = 1.0/2.0*(log(obs_lambda)-lnlm)*(log(obs_lambda)-lnlm)/pow(sigma_total,2.0);
  
  return exp(-x)/M_SQRTPI/M_SQRT2/sigma_total/obs_lambda;
}

// \int_(bin_lambda_obs_min)^(bin_lambda_obs_max) \dlambda_obs P(\lambda_obs|M)
// (see for example https://arxiv.org/pdf/1810.09456.pdf - eq3 qnd 6) 
double buzzard_binned_P_lambda_obs_given_M_nointerp(int nl, double M, double z)
{
  double params[2] = {M, z};

  const double bin_lambda_obs_min = Cluster.N_min[nl];
  const double bin_lambda_obs_max = Cluster.N_max[nl];

  return int_gsl_integrate_medium_precision(buzzard_P_lambda_obs_given_M_nointerp,
    (void*) params, bin_lambda_obs_min, bin_lambda_obs_max, NULL, 1000);
}

// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// SDSS
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------

// Cocoa: we try to avoid reading of files in the cosmolike_core code 
// Cocoa: (reading is done in the C++/python interface)
void setup_SDSS_P_true_lambda_given_mass(int* io_nintrinsic_sigma, double** io_intrinsic_sigma, 
int* io_natsrgm, double** io_atsrgm, double** io_alpha, double** io_sigma, int io)
{
  static int nintrinsic_sigma;
  static int natsrgm;  
  static double* intrinsic_sigma = NULL;
  static double* atsrgm = NULL;
  static double* alpha = NULL;
  static double* sigma = NULL;

  if (io == 1) // IO == 1 IMPLES THAT IO_XXX WILL COPIED TO LOCAL XXX
  {
    if (intrinsic_sigma != NULL)
    {
      free(intrinsic_sigma);
    }
    else if (atsrgm != NULL)
    {
      free(atsrgm);
    }
    else if (alpha != NULL)
    {
      free(alpha);
    }
    else if (sigma != NULL)
    {
      free(sigma);
    }

    nintrinsic_sigma = (*io_nintrinsic_sigma);
    natsrgm = (*io_natsrgm);

    if (!(nintrinsic_sigma > 5) || !(natsrgm > 5))
    {
      log_fatal("array to small for 2D interpolation");
      exit(1);
    }

    intrinsic_sigma = (double*) malloc(nintrinsic_sigma*sizeof(double));
    atsrgm = (double*) malloc(natsrgm*sizeof(double));
    alpha = (double*) malloc(natsrgm*nintrinsic_sigma*sizeof(double));
    sigma = (double*) malloc(natsrgm*nintrinsic_sigma*sizeof(double));
    if(intrinsic_sigma == NULL || atsrgm == NULL || alpha == NULL || sigma == NULL) 
    {
      log_fatal("fail allocation");
      exit(1);
    }
    
    for (int i=0; i<nintrinsic_sigma; i++)
    {
      intrinsic_sigma[i] = (*io_intrinsic_sigma)[i];
      for (int j=0; j<natsrgm; j++)
      {
        if (i == 0) 
        {
          atsrgm[j] = (*io_atsrgm)[j];
        }
        alpha[i*nintrinsic_sigma + j] = (*io_alpha)[i*nintrinsic_sigma + j];
        sigma[i*nintrinsic_sigma + j] = (*io_sigma)[i*nintrinsic_sigma + j];
      }
    }
  }
  else  // IO != 1 IMPLES THAT LOCAL XXX WILL BE COPIED TO IO_XXX
  {
    if (intrinsic_sigma == NULL || atsrgm == NULL || alpha == NULL || sigma == NULL ||
        io_intrinsic_sigma == NULL || io_atsrgm == NULL || io_alpha == NULL || io_sigma == NULL)
    {
      log_fatal("array/pointers not allocated\n");
      exit(1);
    }
    if((*io_intrinsic_sigma) != NULL)
    {
      free((*io_intrinsic_sigma));
      (*io_intrinsic_sigma) = NULL;
    }
    else if((*io_atsrgm) != NULL)
    {
      free((*io_atsrgm));
      (*io_atsrgm) = NULL;
    }
    else if((*io_alpha) != NULL)
    {
      free((*io_alpha));
      (*io_alpha) = NULL;
    }
    else if((*io_sigma) != NULL)
    {
      free((*io_sigma));
      (*io_sigma) = NULL;
    }

    (*io_intrinsic_sigma) = (double*) malloc(nintrinsic_sigma*sizeof(double));
    (*io_atsrgm) = (double*) malloc(natsrgm*sizeof(double));
    (*io_alpha) = (double*) malloc(nintrinsic_sigma*natsrgm*sizeof(double));
    (*io_sigma) = (double*) malloc(nintrinsic_sigma*natsrgm*sizeof(double));
    if(*io_intrinsic_sigma == NULL || *io_atsrgm == NULL || *io_alpha == NULL || *io_sigma == NULL) 
    {
      log_fatal("fail allocation");
      exit(1);
    }

    for (int i=0; i<nintrinsic_sigma; i++)
    {
      (*io_intrinsic_sigma)[i] = intrinsic_sigma[i];
      for (int j=0; j<natsrgm; j++)
      {
        if (i == 0) 
        {
          (*io_atsrgm)[j] = atsrgm[j];
        }
        (*io_alpha)[i*nintrinsic_sigma + j] = alpha[i*nintrinsic_sigma + j];
        (*io_sigma)[i*nintrinsic_sigma + j] = sigma[i*nintrinsic_sigma + j];
      }
    }
  }
}

// SKEW-NORMAL APPROXIMATION eq B1 of https://arxiv.org/pdf/1810.09456.pdf
double SDSS_P_true_lambda_given_mass(const double true_lambda, const double mass, const double z)
{
  static int first = 0;
  static gsl_spline2d* falpha = NULL; // skewness of the skew-normal distribution
  static gsl_spline2d* fsigma = NULL; // variance of the skew-normal distribution

  if(first == 0) 
  {
    first = 1;
    
    int nintrinsic_sigma;
    int natsrgm;
    double** intrinsic_sigma;
    double** atsrgm;
    double** tmp_alpha;
    double** tmp_sigma;
    double* alpha;
    double* sigma;

    intrinsic_sigma = (double**) malloc(1*sizeof(double*));
    atsrgm = (double**) malloc(1*sizeof(double*));
    tmp_alpha = (double**) malloc(1*sizeof(double*));
    tmp_sigma = (double**) malloc(1*sizeof(double*));
    if (intrinsic_sigma == NULL || atsrgm == NULL || tmp_alpha == NULL || tmp_sigma == NULL)
    {
      log_fatal("fail allocation");
      exit(1);
    }
    else 
    {
      (*intrinsic_sigma) = NULL;
      (*atsrgm) = NULL;
      (*tmp_alpha) = NULL;
      (*tmp_sigma) = NULL; 
    }
   
    setup_SDSS_P_true_lambda_given_mass(&nintrinsic_sigma, intrinsic_sigma, &natsrgm, atsrgm, 
      tmp_alpha, tmp_sigma, 0);

    const gsl_interp2d_type* T = gsl_interp2d_bilinear;
    falpha = gsl_spline2d_alloc(T, nintrinsic_sigma, natsrgm);
    fsigma = gsl_spline2d_alloc(T, nintrinsic_sigma, natsrgm);
    if (falpha == NULL || fsigma == NULL)
    {
      log_fatal("fail allocation");
      exit(1);
    }

    // we don't want to guess the appropriate GSL z array ordering in z = f(x, y)
    alpha = (double*) malloc(nintrinsic_sigma*natsrgm*sizeof(double));
    sigma = (double*) malloc(nintrinsic_sigma*natsrgm*sizeof(double));
    if (alpha == NULL || sigma == NULL)
    {
      log_fatal("fail allocation");
      exit(1);
    }
    
    for (int i=0; i<nintrinsic_sigma; i++) 
    {
      for (int j=0; j<natsrgm; j++) 
      {
        int status = 0;
        status = gsl_spline2d_set(falpha, alpha, i, j, (*tmp_alpha)[i*nintrinsic_sigma+j]);
        if (status) 
        {
          log_fatal(gsl_strerror(status));
          exit(1);
        }
        status = gsl_spline2d_set(fsigma, sigma, i, j, (*tmp_sigma)[i*nintrinsic_sigma+j]);
        if (status) 
        {
          log_fatal(gsl_strerror(status));
          exit(1);
        }
      }
    }

    int status = 0;
    status = gsl_spline2d_init(falpha, intrinsic_sigma[0], atsrgm[0], alpha, 
      nintrinsic_sigma, natsrgm);
    if (status) 
    {
      log_fatal(gsl_strerror(status));
      exit(1);
    }
    status = gsl_spline2d_init(fsigma, intrinsic_sigma[0], atsrgm[0], sigma, 
      nintrinsic_sigma, natsrgm);
    if (status) 
    {
      log_fatal(gsl_strerror(status));
      exit(1);
    }

    free(intrinsic_sigma[0]);
    free(atsrgm[0]);
    free(tmp_alpha[0]);
    free(tmp_sigma[0]);
    free(intrinsic_sigma);
    free(atsrgm);
    free(tmp_alpha);
    free(tmp_sigma);
    free(alpha); // GSL SPLINE 2D copies the array
    free(sigma); // GSL SPLINE 2D copies the array
  }

  const double mass_min = pow(10.0, nuisance.cluster_MOR[0]);
  const double mass_M1 = pow(10.0, nuisance.cluster_MOR[1]);
  const double intrinsic_alpha = nuisance.cluster_MOR[2];
  const double intrinsic_sigma = nuisance.cluster_MOR[3]; //intrisic scatter, mass-richness relation
  
  // atsrgm = average true satellite richness given mass
  const double tmp = (mass - mass_min)/(mass_M1 - mass_min);
  double atsrgm = pow(tmp, intrinsic_alpha)*pow(((1+z)/1.45), nuisance.cluster_MOR[4]); 
  if (atsrgm > 160) 
  {
    atsrgm = 160;
  }
  else if (atsrgm < 1) 
  {
    atsrgm = 1;
  }
  
  int status = 0;
  double alpha = 0.0;
  double sigma = 0.0;
  status = gsl_spline2d_eval_e(falpha, intrinsic_sigma, atsrgm, NULL, NULL, &alpha);
  if (status) 
  {
    log_fatal(gsl_strerror(status));
    exit(1);
  }
  status = gsl_spline2d_eval_e(fsigma, intrinsic_sigma, atsrgm, NULL, NULL, &sigma);
  if (status) 
  {
    log_fatal(gsl_strerror(status));
    exit(1);
  }

  const double x = 1.0/(M_SQRT2*abs(sigma));
  const double result = exp(-(true_lambda-atsrgm)*(true_lambda-atsrgm)*x*x)*x/M_SQRTPI;
  return result*gsl_sf_erfc(-alpha*(true_lambda-atsrgm)*x);
}

// Cocoa: we try to avoid reading of files in the cosmolike_core code 
// Cocoa: (reading is done in the C++/python interface)
void setup_SDSS_P_lambda_obs_given_true_lambda(int* io_nz, double** io_z, int* io_nlambda, 
double** io_lambda, double** io_tau, double** io_mu, double** io_sigma, double** io_fmask, 
double** io_fprj, int io)
{
  static int nz;
  static int nlambda; 
  static double* z = NULL; 
  static double* lambda = NULL;
  static double* tau = NULL;
  static double* mu = NULL;
  static double* sigma = NULL;
  static double* fmask = NULL;
  static double* fprj = NULL;

  if (io == 1) // IO == 1 IMPLES THAT IO_XXX WILL COPIED TO LOCAL XXX
  {
    if (z != NULL)
    {
      free(z);
    }
    else if (lambda != NULL)
    {
      free(lambda);
    }
    else if (tau != NULL)
    {
      free(tau);
    }
    else if (mu != NULL)
    {
      free(mu);
    }
    else if (sigma != NULL)
    {
      free(sigma);
    }
    else if (fmask != NULL)
    {
      free(fmask);
    }
    else if (fprj != NULL)
    {
      free(fprj);
    }

    nz = (*io_nz);
    nlambda = (*io_nlambda);
    if (!(nz > 5) || !(nlambda > 5))
    {
      log_fatal("array to small for 2D interpolation");
      exit(1);
    }

    z = (double*) malloc(nz*sizeof(double));
    lambda = (double*) malloc(nlambda*sizeof(double));
    tau = (double*) malloc(nz*nlambda*sizeof(double));
    mu = (double*) malloc(nz*nlambda*sizeof(double));
    sigma = (double*) malloc(nz*nlambda*sizeof(double));
    fmask = (double*) malloc(nz*nlambda*sizeof(double));
    fprj = (double*) malloc(nz*nlambda*sizeof(double));
    if(z == NULL || lambda == NULL || tau == NULL || mu == NULL || sigma == NULL || 
       fmask == NULL || fprj == NULL) 
    {
      log_fatal("fail allocation");
      exit(1);
    }

    for (int i=0; i<nz; i++)
    {
      z[i] = (*io_z)[i];
      for (int j=0; j<nlambda; j++)
      {
        if (i == 0) 
        {
          lambda[j] = (*io_lambda)[j];
        }
        tau[i*nz + j] = (*io_tau)[i*nz + j];
        mu[i*nz + j] = (*io_mu)[i*nz + j];
        sigma[i*nz + j] = (*io_sigma)[i*nz + j];
        fmask[i*nz + j] = (*io_fmask)[i*nz + j];
        fprj[i*nz + j] = (*io_fprj)[i*nz + j];
      }
    }
  }
  else
  {
    // IO != 1 IMPLES THAT LOCAL H(Z) WILL BE COPIED TO IO_chi(Z)
    if (z == NULL || lambda == NULL || tau == NULL || mu == NULL || sigma == NULL || fmask == NULL 
        || fprj == NULL || io_lambda == NULL || io_tau == NULL || io_mu == NULL || io_sigma == NULL
        || io_fmask == NULL || io_fprj == NULL)
    {
      log_fatal("array/pointer not allocated");
      exit(1);
    }
    if((*io_z) != NULL)
    {
      free((*io_z));
      (*io_z) = NULL;
    }
    else if((*io_lambda) != NULL)
    {
      free((*io_lambda));
      (*io_lambda) = NULL;
    }
    else if((*io_tau) != NULL)
    {
      free((*io_tau));
      (*io_tau) = NULL;
    }
    else if((*io_mu) != NULL)
    {
      free((*io_mu));
      (*io_mu) = NULL;
    }
    else if((*io_sigma) != NULL)
    {
      free((*io_sigma));
      (*io_sigma) = NULL;
    }
    else if((*io_fmask) != NULL)
    {
      free((*io_fmask));
      (*io_fmask) = NULL;
    }   
    else if((*io_fprj) != NULL)
    {
      free((*io_fprj));
      (*io_fprj) = NULL;
    }
   
    (*io_z) = (double*) malloc(nz*sizeof(double));
    (*io_lambda) = (double*) malloc(nlambda*sizeof(double));
    (*io_tau) = (double*) malloc(nz*nlambda*sizeof(double));
    (*io_mu) = (double*) malloc(nz*nlambda*sizeof(double));
    (*io_sigma) = (double*) malloc(nz*nlambda*sizeof(double));
    (*io_fmask) = (double*) malloc(nz*nlambda*sizeof(double));
    (*io_fprj) = (double*) malloc(nz*nlambda*sizeof(double));
    if((*io_z) == NULL || (*io_lambda) == NULL || (*io_tau) == NULL || (*io_mu) == NULL ||
       (*io_sigma) == NULL || (*io_fmask) == NULL || (*io_fprj) == NULL) 
    {
      log_fatal("fail allocation");
      exit(1);
    }

    for (int i=0; i<nz; i++)
    {
      (*io_z)[i] = z[i];
      for (int j=0; j<nlambda; j++)
      {
        if (i == 0) 
        {
          (*io_lambda)[j] = lambda[j];
        }
        (*io_tau)[i*nz + j] = tau[i*nz + j];
        (*io_mu)[i*nz + j]  = mu[i*nz + j];
        (*io_sigma)[i*nz + j] = sigma[i*nz + j];
        (*io_fmask)[i*nz + j] = fmask[i*nz + j];
        (*io_fprj)[i*nz + j] = fprj[i*nz + j];
      }
    }
  }
}

double SDSS_P_lambda_obs_given_true_lambda(double observed_lambda, double true_lambda, double zz) 
{
  static int first = 0;
  static gsl_spline2d* ftau;
  static gsl_spline2d* fmu;
  static gsl_spline2d* fsigma;
  static gsl_spline2d* ffmask;
  static gsl_spline2d* ffprj;

  if(first == 0) 
  {
    first = 1;

    int nz;
    int nlambda;
    double** tmp_tau;
    double** tmp_mu;
    double** tmp_sigma;
    double** tmp_fmask;
    double** tmp_fprj;
    double** z;
    double** lambda;
    double* tau;
    double* mu;
    double* sigma;
    double* fmask;
    double* fprj;

    z = (double**) malloc(1*sizeof(double*));
    lambda = (double**) malloc(1*sizeof(double*));
    tmp_tau = (double**) malloc(1*sizeof(double*));
    tmp_mu = (double**) malloc(1*sizeof(double*));
    tmp_sigma = (double**) malloc(1*sizeof(double*));
    tmp_fmask = (double**) malloc(1*sizeof(double*));
    tmp_fprj = (double**) malloc(1*sizeof(double*));
    if (z == NULL || lambda == NULL || tmp_tau == NULL || tmp_mu == NULL || tmp_sigma == NULL
        || tmp_fmask == NULL || tmp_fprj == NULL)
    {
      log_fatal("fail allocation");
      exit(1);
    }
    else 
    {
      (*z) = NULL;
      (*lambda) = NULL;
      (*tmp_tau) = NULL;
      (*tmp_mu) = NULL; 
      (*tmp_sigma) = NULL; 
      (*tmp_fmask) = NULL; 
      (*tmp_fprj) = NULL; 
    }
   
    setup_SDSS_P_lambda_obs_given_true_lambda(&nz, z, &nlambda, lambda, tmp_tau, tmp_mu, tmp_sigma, 
      tmp_fmask, tmp_fprj, 0);

    const gsl_interp2d_type *T = gsl_interp2d_bilinear;
    ftau = gsl_spline2d_alloc(T, nz, nlambda);
    fmu = gsl_spline2d_alloc(T, nz, nlambda);
    fsigma = gsl_spline2d_alloc(T, nz, nlambda);
    ffmask = gsl_spline2d_alloc(T, nz, nlambda);    
    ffprj = gsl_spline2d_alloc(T, nz, nlambda);
    if (ftau == NULL || fmu == NULL || fsigma == NULL || ffmask == NULL || ffprj == NULL)
    {
      log_fatal("fail allocation");
      exit(1);
    }

    // we don't want to guess the appropriate GSL z array ordering in z = f(x, y)
    tau = (double*) malloc(nz*nlambda*sizeof(double));
    mu = (double*) malloc(nz*nlambda*sizeof(double));
    sigma = (double*) malloc(nz*nlambda*sizeof(double));
    fmask = (double*) malloc(nz*nlambda*sizeof(double));
    fprj = (double*) malloc(nz*nlambda*sizeof(double));
    if (tau == NULL || mu == NULL || sigma == NULL || fmask == NULL || fprj == NULL)
    {
      log_fatal("fail allocation");
      exit(1);
    }

    #pragma omp parallel for collapse(2)
    for (int i=0; i<nz; i++) 
    {
      for (int j=0; j<nlambda; j++) 
      {
        int status = 0;
        status = gsl_spline2d_set(ftau, tau, i, j, (*tmp_tau)[i*nz + j]);
        if (status) 
        {
          log_fatal(gsl_strerror(status));
          exit(1);
        }
        status = gsl_spline2d_set(fmu, mu, i, j, (*tmp_mu)[i*nz + j]);
        if (status) 
        {
          log_fatal(gsl_strerror(status));
          exit(1);
        }
        status = gsl_spline2d_set(fsigma, sigma, i, j, (*tmp_sigma)[i*nz + j]);
        if (status) 
        {
          log_fatal(gsl_strerror(status));
          exit(1);
        }
        status = gsl_spline2d_set(ffmask, fmask, i, j, (*tmp_fmask)[i*nz + j]);
        if (status) 
        {
          log_fatal(gsl_strerror(status));
          exit(1);
        }
        status = gsl_spline2d_set(ffprj, fprj, i, j, (*tmp_fprj)[i*nz + j]);
        if (status) 
        {
          log_fatal(gsl_strerror(status));
          exit(1);
        }
      }
    }

    int status = 0;
    status = gsl_spline2d_init(ftau, (*z), (*lambda), tau, nz, nlambda);
    if (status) 
    {
      log_fatal(gsl_strerror(status));
      exit(1);
    }
    status = gsl_spline2d_init(fmu, (*z), (*lambda), mu, nz, nlambda);
    if (status) 
    {
      log_fatal(gsl_strerror(status));
      exit(1);
    }
    status = gsl_spline2d_init(fsigma, (*z), (*lambda), sigma, nz, nlambda);
    if (status) 
    {
      log_fatal(gsl_strerror(status));
      exit(1);
    }
    status = gsl_spline2d_init(ffmask, (*z), (*lambda), fmask, nz, nlambda);
    if (status) 
    {
      log_fatal(gsl_strerror(status));
      exit(1);
    }
    status = gsl_spline2d_init(ffprj, (*z), (*lambda), fprj, nz, nlambda);
    if (status) 
    {
      log_fatal(gsl_strerror(status));
      exit(1);
    }

    free(z[0]);
    free(lambda[0]);
    free(tmp_tau[0]);
    free(tmp_mu[0]);
    free(tmp_sigma[0]);
    free(tmp_fmask[0]);
    free(tmp_fprj[0]);
    free(z);
    free(lambda);
    free(tmp_tau);
    free(tmp_mu);
    free(tmp_sigma);
    free(tmp_fmask);
    free(tmp_fprj);
    free(tau);   // GSL SPLINE 2D copies the array
    free(mu);    // GSL SPLINE 2D copies the array
    free(sigma); // GSL SPLINE 2D copies the array
    free(fmask); // GSL SPLINE 2D copies the array
    free(fprj);  // GSL SPLINE 2D copies the array
  }
  

  double tau, mu, sigma, fmask, fprj;
  {
    int status = 0;
    status = gsl_spline2d_eval_e(ftau, zz, true_lambda, NULL, NULL, &tau);
    if (status) 
    {
      log_fatal(gsl_strerror(status));
      exit(1);
    }  
    status = gsl_spline2d_eval_e(fmu, zz, true_lambda, NULL, NULL, &mu);
    if (status) 
    {
      log_fatal(gsl_strerror(status));
      exit(1);
    }  
    status = gsl_spline2d_eval_e(fsigma, zz, true_lambda, NULL, NULL, &sigma);
    if (status) 
    {
      log_fatal(gsl_strerror(status));
      exit(1);
    }
    status = gsl_spline2d_eval_e(ffmask, zz, true_lambda, NULL, NULL, &fmask);
    if (status) 
    {
      log_fatal(gsl_strerror(status));
      exit(1);
    }
    status = gsl_spline2d_eval_e(ffprj, zz, true_lambda, NULL, NULL, &fprj);
    if (status) 
    {
      log_fatal(gsl_strerror(status));
      exit(1);
    }
  }

  const double x = 1.0/(M_SQRT2*abs(sigma));
  const double y = 1.0/true_lambda;
  const double j = exp(0.5*tau*(2*mu+tau*sigma*sigma-2*observed_lambda));

  double r0 = exp(-1.*(observed_lambda-mu)*(observed_lambda-mu)*x*x);
  r0 *= (1-fmask)*(1-fprj)*x/M_SQRTPI;

  double r1 = 0.5*((1-fmask)*fprj*tau + fmask*fprj*y)*j;
  r1 *= gsl_sf_erfc((mu+tau*sigma*sigma-observed_lambda)*x);

  double r2 = 0.5*fmask*y;
  r2 *= gsl_sf_erfc((mu-observed_lambda-true_lambda)*x) - gsl_sf_erfc((mu-observed_lambda)*x);

  double r3 = 0.5*fmask*fprj*y*exp(-tau*true_lambda)*j;
  r3 *= gsl_sf_erfc((mu+tau*sigma*sigma-observed_lambda-true_lambda)*x);
  
  return r0 + r1 + r2 - r3;
}

double SDSS_P_lambda_obs_lambda_true_given_M(double true_lambda, const double obs_lambda, 
const double M, const double z)
{
  const double r1 = SDSS_P_lambda_obs_given_true_lambda(obs_lambda, true_lambda, z);
  const double r2 = SDSS_P_true_lambda_given_mass(true_lambda, M, z);
  return r1*r2;
}

double SDSS_P_lambda_obs_lambda_true_given_M_wrapper(double true_lambda, void* params)
{
  return SDSS_P_lambda_obs_lambda_true_given_M(true_lambda, obs_lambda, M, z);
}

double SDSS_P_lambda_obs_given_M(double obs_lambda, const double M, const double z)
{  
  double params_in[3] = {M, z, obs_lambda}; 
  const double true_lambda_min =  3.;
  const double true_lambda_max =  160.;
  return int_gsl_integrate_medium_precision(SDSS_P_lambda_obs_lambda_true_given_M_wrapper, 
      (void*) params_in, true_lambda_min, true_lambda_max, NULL, GSL_WORKSPACE_SIZE);
}

double SDSS_P_lambda_obs_given_M_wrapper(double obs_lambda, void* params)
{
  double* ar = (double*) params; 
  const double M = ar[0];
  const double z = ar[1];
  return SDSS_P_lambda_obs_given_M(obs_lambda, M, z);
}

double SDSS_binned_P_lambda_obs_given_M(int nl, double M, double z)
{
  double params[2] = {M, z};
  const int nl_min = Cluster.N_min[nl];
  const int nl_max = Cluster.N_max[nl];
  return int_gsl_integrate_medium_precision(SDSS_P_lambda_obs_given_M_wrapper, (void*) params, 
    (double) nl_min, (double) nl_max, NULL, GSL_WORKSPACE_SIZE);
}

// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// INTERFACE
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------

// \int_(bin_lambda_obs_min)^(bin_lambda_obs_max) \dlambda_obs P(\lambda_obs|M)
// (see for example https://arxiv.org/pdf/1810.09456.pdf - eq 3 and eq 6) 
double binned_P_lambda_obs_given_M_noiterp(int nl, double M, double z)
{
  if (strcmp(Cluster.model, "SDSS") == 0) 
  {
    return SDSS_binned_P_lambda_obs_given_M(nl, M, z);
    /*
      const int nsizex = 20;
      const int nsizey = 20;
      const double true_lambda_min = 3.;
      const double true_lambda_max = 160.;
      const double obs_lambda_min = Cluster.N_min[nl];
      const double obs_lambda_max = Cluster.N_max[nl];
      const int dobs_lambda = (obs_lambda_max-obs_lambda_min)/((double) nsizex - 1);
      const int dtrue_lambda = (true_lambda_max-true_lambda_min)/((double) nsizey - 1);

      for(int i=0; i<nsizex; i++) // integrate on lambda_obs inside the bin nl
      {
        for(int j=0; j<nsizey; j++) // integrate on lambda_true
        {
          const double obs_lambda = obs_lambda_min + i*dobs_lambda;
          const double true_lambda = true_lambda_min + i*dtrue_lambda;
          const double wx = (i == 0 || nsizex - 1) ? 0.5*dobs_lambda : dobs_lambda;
          const double wy = (j == 0 || nsizey - 1) ? 0.5*dtrue_lambda : dtrue_lambda;
          sum[i*nsizey + j] = wx*wy*SDSS_P_lambda_obs_lambda_true_given_M(true_lambda, obs_lambda, M, z);
        }
      }
      double sum = 0.0;
      for(int i=0; i<nsizex; i++) // integrate on lambda_obs inside the bin nl
      {
        for(int j=0; j<nsizey; j++) // integrate on lambda_true
        {
          sum += sum[i*nsizey + j];
        }
      }
      return sum;
    */
  }
  else if (strcmp(Cluster.model, "BUZZARD") == 0)
  {
    return buzzard_binned_P_lambda_obs_given_M(nl, M, z);
  }
  else
  {
    log_fatal("Cluster.model not implemented");
    exit(1);
  }
}

double binned_P_lambda_obs_given_M(int nl, double M, double z) 
{
  static cosmopara C;
  static nuisancepara N;
  static double*** table = 0;

  const int nsizez = binned_P_lambda_obs_given_M_size_z_table;
  const int nsizem = binned_P_lambda_obs_given_M_size_M_table;
  const int nsizel = Cluster.N200_Nbin;
  const double logmmin = log_M_min;
  const double logmmax = log_M_max; 
  const double zmin = binned_P_lambda_obs_given_M_zmin_table; 
  const double zmax = binned_P_lambda_obs_given_M_zmax_table;
  const double dm = (logmmax - logmmin)/((double) nm - 1.0);
  const double dz = (zmax-zmin)/((double) nsizez - 1.0);
  const double mmin = pow(10.0, logmmin);
  const double mmax = pow(10.0, logmmax);
  
  if (table == 0)
  {
    table = (double***) malloc(sizeof(double**)*nsizel);
    for(int i=0; i<nz; i++)
    {
      table[i] = (double**) malloc(sizeof(double*)*nsizez);
      for(int j=0; j<nz; i++)
      {
        table[i][j] = (double*) malloc(sizeof(double)*nsizem);
      }
    }
  }
  if (recompute_clusters(C, N))
  { 
    if (strcmp(Cluster.model, "SDSS") == 0) // init static vars
    {
      SDSS_P_lambda_obs_given_true_lambda(50, 50, z);
      SDSS_P_true_lambda_given_mass(50, M, z);
    }
    #pragma omp parallel for collapse(3)
    for (int i=0; i<nsizel; i++) 
    {
      for (int j=0; j<nsizez; j++) 
      {
        for (int k=0; k<nsizem; k++) 
        {
          table[i][j][k] = 
            binned_P_lambda_obs_given_M_noiterp(i, pow(10, logmmin + k*dm), zmin + j*dz);
        }
      }
    }
    update_cosmopara(&C);
    update_nuisance(&N);
  }
  if (z > zmin && z < zmax && M > mmin && M < mmax)
  {
    const double logM = log10(M);
    return interpol2d(table[nl], nz, zmin, zmax, dz, z, nm, logmmin, logmmax, dm, logM, 1, 1);
  }
  else
  {
    return 0.0;
  }
}

// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// CLUSTER BIAS
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------

double B1_x_BSF(double M, double a)
{ // cluster bias including selection bias
  double B1;
  if (Cluster.bias_model == 0) 
  {
    B1 = B1(mass,a);
  } 
  else if (Cluster.bias_model == 1) 
  {
    B1 = tinker_emu_B1_tab(M, a);
  }
  else 
  {
    log_fatal("Cluster bias model %i not implemented", Cluster.bias_model);
    exit(0); 
  }
   
  double BS;
  if (Cluster.N_SF == 1) 
  {
    BS = nuisance.cluster_selection[0]; 
  }
  else if (Cluster.N_SF == 2)
  {
    BS = nuisance.cluster_selection[0]*pow((M/M_PIVOT), nuisance.cluster_selection[1]); 
  }
  else
  {
    log_fatal("Cluster selection bias model %i not implemented", Cluster.bias_model);
    exit(0); 
  }

  return B1 * BS;
}

double B2_x_BSF(double mass, double a)
{ // cluster bias (b2) including selection bias

  double B1;
  if (Cluster.bias_model == 0) 
  {
    B1 = B1(mass,a);
  } 
  else if (Cluster.bias_model == 1) 
  {
    B1 = tinker_emu_B1_tab(M, a);
  }
  else 
  {
    log_fatal("Cluster bias model %i not implemented", Cluster.bias_model);
    exit(0); 
  }

  const double B2 = b2_from_b1(B1);
  
  double BS;
  if (Cluster.N_SF == 1) 
  {
    BS = nuisance.cluster_selection[0]; 
  }
  else if (Cluster.N_SF == 2)
  {
    BS = nuisance.cluster_selection[0]*pow((M/M_PIVOT), nuisance.cluster_selection[1]); 
  }
  else
  {
    log_fatal("Cluster selection bias model %i not implemented", Cluster.bias_model);
    exit(0); 
  }
  
  return B2 * BS;
}

double B1M1_x_BS(double mass, double a)
{
  double B1M1;
  if (Cluster.bias_model == 0) 
  {
    B1M1 = B1(mass,a) -1;
  } 
  else if (Cluster.bias_model == 1) 
  {
    B1M1 = tinker_emu_B1_tab(M, a) -1;
  }
  else 
  {
    log_fatal("Cluster bias model %i not implemented", Cluster.bias_model);
    exit(0); 
  }

  double BS;
  if (Cluster.N_SF == 1) 
  {
    BS = nuisance.cluster_selection[0]; 
  }
  else if (Cluster.N_SF == 2)
  {
    BS = nuisance.cluster_selection[0]*pow((M/M_PIVOT), nuisance.cluster_selection[1]); 
  }
  else
  {
    log_fatal("Cluster selection bias model %i not implemented", Cluster.bias_model);
    exit(0); 
  }

  return B1M1 * BS;
}

double dndlnM_times_binned_P_lambda_obs_given_M(double lnM, void* params)
{
  double *ar = (double*) params; 
  const int nl = (int) ar[0];
  const double z = ar[1];
  const double M = exp(lnM);
  const double a = 1./(1. + z);

  double mfunc; 
  if (Cluster.hmf_model == 0)
  {
    mfunc = massfunc(M, a);
  }
  else if (Cluster.hmf_model == 1)
  {
    mfunc = tinker_emu_massfunc_tab(M, a);
  }
  else
  {
    log_fatal("massfunc model %i not implemented", Cluster.hmf_model);
    exit(1); 
  }
  return mfunc*M*binned_P_lambda_obs_given_M(nl, M , z);
}

double int_for_weighted_B1(double lnM, void* params)
{
  double *ar = (double*) params; //nl, z
  const double z = ar[1];
  const double M = exp(lnM);  
  const double a = 1./(1. + z);
  return B1_x_BSF(M, a)*dndlnM_times_binned_P_lambda_obs_given_M(lnM, params); 
}

double weighted_B1_nointerp(int nl, double z, int init_static_vars_only)
{  
  const double ln_M_min = 
    (Cluster.N_min[nlambda] > 1E10) ? log(Cluster.N_min[nlambda]) : log_M_min/LOG10_E;
  const double ln_M_max = 
    (Cluster.N_min[nlambda] > 1E10) ? log(Cluster.N_max[nlambda]) : log_M_max/LOG10_E;
 
  double param[2] = {(double) nl, z};

  if(int init_static_vars_only == 1)
  {
    const double r1 = int_for_weighted_B1(ln_M_min, (void*) param);
    const double r2 = dndlnM_times_binned_P_lambda_obs_given_M(M_min, (void*) param);
    return 0.0;
  }
  else
  {
    const double r1 = int_gsl_integrate_low_precision(int_for_weighted_B1, 
      (void*) param, ln_M_min, ln_M_max, NULL, GSL_WORKSPACE_SIZE);
  
    const double r2 = int_gsl_integrate_low_precision(dndlnM_times_binned_P_lambda_obs_given_M, 
      (void*) param, ln_M_min, ln_M_max, NULL, GSL_WORKSPACE_SIZE);
  }
  return (r2 <= 0) ? 0.0 : r1/r2; 
}

double weighted_B1(int nl, double z)
{
  static cosmopara C;
  static nuisancepara N;
  static double** table;

  const int nasize = Ntable.N_a;
  const int nlsize = Cluster.N200_Nbin;
  const double zmin = fmax(tomo.cluster_zmin[0]-0.05, 0.01); 
  const double zmax = tomo.cluster_zmax[tomo.cluster_Nbin-1] + 0.05;
  const double amin = 1./(1.+zmax); 
  const double amax = 1./(1.+zmin);
  const double da = (amax-amin)/(nasize-1.);

  if (table == 0) 
  {
    table = (double**) malloc(sizeof(double*)*nlsize);
    for(int i=0; i<nlsize; i++)
    {
      table[i] = (double*) malloc(sizeof(double)*nasize);
    }
  }
  if (recompute_clusters(C, N))
  {    
    weighted_B1_nointerp(0, 1.0/amin - 1.0, 1); // init static vars only
    #pragma omp parallel for collapse(2)
    for (int n=0; n<nlsize; n++)
    {
      for (int i=0; i<nasize; i++)
      {
        const double a = amin + i*da;
        table[n][i] = weighted_B1_nointerp(n, 1./a - 1.0, 0);
      }
    }
    update_cosmopara(&C);
    update_nuisance(&N);
  }
  const double a = 1./(z + 1.);
  return interpol(table[nl], nasize, amin, amax, da, a, 0., 0.);
}

double int_for_weighted_B2(double lnM, void* params)
{
  double *ar = (double*) params; //nl, z
  const double z = ar[1];
  const double M = exp(lnM);  
  const double a = 1./(1. + z);
  return B2_x_BSF(M, a)*dndlnM_times_binned_P_lambda_obs_given_M(lnM, params); 
}

double weighted_B2_nointerp(int nl, double z, int init_static_vars_only)
{  
  const double ln_M_min = 
    (Cluster.N_min[nlambda] > 1E10) ? log(Cluster.N_min[nlambda]) : log_M_min/LOG10_E; 
  const double ln_M_max = 
    (Cluster.N_min[nlambda] > 1E10) ? log(Cluster.N_max[nlambda]) : log_M_max/LOG10_E; 
 
   double param[2] = {(double) nl, z};

  if(int init_static_vars_only == 1)
  {
    const double r1 = int_for_weighted_B2(ln_M_min, (void*) param);
    const double r2 = dndlnM_times_binned_P_lambda_obs_given_M(ln_M_min, (void*) param);
    return 0.0;
  }
  else
  {
    const double r1 = int_gsl_integrate_low_precision(int_for_weighted_B2, 
      (void*) param, ln_M_min, ln_M_max, NULL, GSL_WORKSPACE_SIZE);

    const double r2 = int_gsl_integrate_low_precision(dndlnM_times_binned_P_lambda_obs_given_M, 
      (void*) param, ln_M_min, ln_M_max, NULL, GSL_WORKSPACE_SIZE);

    return (r2 <= 0) ? 0.0 : r1/r2; 
  }
}

double weighted_B2(int nl, double z)
{
  static cosmopara C;
  static nuisancepara N;
  static double** table;

  const int nasize = Ntable.N_a;
  const int nlsize = Cluster.N200_Nbin;
  const double zmin = fmax(tomo.cluster_zmin[0]-0.05, 0.01); 
  const double zmax = tomo.cluster_zmax[tomo.cluster_Nbin-1] + 0.05;
  const double amin = 1./(1.+zmax); 
  const double amax = 1./(1.+zmin);
  const double da = (amax-amin)/(nasize-1.);

  if (table == 0) 
  {
    table = (double**) malloc(sizeof(double*)*nlsize);
    for(int i=0; i<nlsize; i++)
    {
      table[i] = (double*) malloc(sizeof(double)*nasize);
    }
  }
  if (recompute_clusters(C, N))
  {    
    weighted_B2_nointerp(0, 1.0/amin - 1., 1); // init static vars only
    #pragma omp parallel for collapse(2)
    for (int n=0; n<nlsize; n++)
    {
      for (int i=0; i<nasize; i++)
      {
        const double a = amin + i*da;
        table[n][i] = weighted_B2_nointerp(n, 1./a - 1.0, 0);
      }
    }
    update_cosmopara(&C);
    update_nuisance(&N);
  }
  const double a = 1./(z + 1.);
  return interpol(table[nl], nasize, amin, amax, da, a, 0., 0.);
}

double int_for_weighted_B1M1(double lnM, void* params)
{
  double *ar = (double*) params;
  const double z = ar[1];
  const double M = exp(lnM);  
  const double a = 1.0/(1.0 + z);
  return B1M1_x_BS(M, a)*dndlnM_times_binned_P_lambda_obs_given_M(lnM, params); 
}

double weighted_B1M1_nointerp(int nl, double z, int init_static_vars_only)
{  
  const double ln_M_min = 
    (Cluster.N_min[nlambda] > 1E10) ? log(Cluster.N_min[nlambda]) : log_M_min/LOG10_E; 
  const double ln_M_max = 
    (Cluster.N_min[nlambda] > 1E10) ? log(Cluster.N_max[nlambda]) : log_M_max/LOG10_E;

  double params[2] = {(double) nl, z};

  if (init_static_vars_only == 1)
  {
    const double r1 = int_for_weighted_B1M1(ln_M_min, (void*) param);
    const double r2 = dndlnM_times_binned_P_lambda_obs_given_M(ln_M_min, (void*) params);
    return 0.0;
  }
  else
  {
    const double r1 = int_gsl_integrate_low_precision(int_for_weighted_B1M1, 
      (void*) param, ln_M_min, ln_M_max, NULL, GSL_WORKSPACE_SIZE);

    const double r2 = int_gsl_integrate_low_precision(dndlnM_times_binned_P_lambda_obs_given_M, 
      (void*) params, ln_M_min, ln_M_max, NULL, GSL_WORKSPACE_SIZE);
    
    return (r2 <= 0) ? 0.0 : r1/r2; 
  }
}

double weighted_B1M1(int nl, double z)
{
  static cosmopara C;
  static nuisancepara N;
  static double** table;

  const int nasize = Ntable.N_a;
  const int nlsize = Cluster.N200_Nbin;
  
  const double zmin = fmax(tomo.cluster_zmin[0]-0.05, 0.01); 
  const double zmax = tomo.cluster_zmax[tomo.cluster_Nbin-1] + 0.05;
  const double amin = 1./(1.+zmax); 
  const double amax = 1./(1.+zmin);
  const double da = (amax-amin)/(nasize-1.);

  if (table == 0) 
  {
    table = (double**) malloc(sizeof(double*)*nlsize);
    for(int i=0; i<nlsize; i++)
    {
      table[i] = (double*) malloc(sizeof(double)*nasize);
    }
  }
  if (recompute_clusters(C, N))
  {    
    weighted_B1M1_nointerp(0, 1.0/amin - 1., 1); // init static vars only
    #pragma omp parallel for collapse(2)
    for (int n=0; n<nlsize; n++)
    {
      for (int i=0; i<nasize; i++)
      {
        const double a = amin + i*da;
        table[n][i] = weighted_B1M1_nointerp(n, 1./a - 1.0, 0);
      }
    }
    update_cosmopara(&C);
    update_nuisance(&N);
  }
  const double a = 1./(z + 1.);
  return interpol(table[nl], nasize, amin, amax, da, a, 0., 0.);
}

// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// cluster number counts
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// nl = lambda_obs bin, ni = cluster redshift bin

double binned_average_number_counts(int nl, double z)
{ // def: eq 3 of https://arxiv.org/pdf/1810.09456.pdf; nl = lambda_obs bin, nz = redshift bin
  double param[2] = {(double) nl, z};
  const double ln_M_min = log_M_min/LOG10_E; 
  const double ln_M_max = log_M_max/LOG10_E; 
  return int_gsl_integrate_low_precision(dndlnM_times_binned_P_lambda_obs_given_M, 
    (void*) param, ln_M_min, ln_M_max, NULL, GSL_WORKSPACE_SIZE);
}

// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// CLUSTER POWER SPECTRUM
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------

// nl{1,2} = lambda_obs bins, n{i,j} = cluster redshift bins
double binned_p_cc(double k, double a, int nl1, int nl2, int use_linear_ps)
{ // binned in lambda_obs (nl1, nl2 = lambda_obs bin (cluster 1 and 2))
  if(!(a>0) || !(a<1)) 
  {
    log_fatal("a>0 and a<1 not true");
    exit(1);
  }
  const double z = 1.0/a - 1;
  const double cluster_B1_1 = weighted_B1(nl1, z); 
  const double PK = use_linear_ps == 1 ? p_lin(k, a) : Pdelta(k, a);

  double P_1loop = 0;
  if (has_b2_galaxies() || (Cluster.nonlinear_bias > 0) && use_linear_ps != 1)
  {
    const double g = growfac(a);
    const double g4 = g*g*g*g;
    
    const double b1g = cluster_B1_1;
    const double b2g = weighted_B2(nl1, z);
    const double bs2g = bs2_from_b1(b1g);
    
    const double b1c = weighted_B1(nl2, z);
    const double b2c = weighted_B2(nl2, z);
    const double bs2c = bs2_from_b1(b1c);
    
    P_1loop = 0.5*(b1c*b2g+b2c*b1g)*PT_d1d2(k) + 0.25*b2g*b2c*PT_d2d2(k) 
      + 0.5*(b1c*bs2g+b1g*bs2c)*PT_d1s2(k) +0.25*(b2c*bs2g+b2g*bs2c)*PT_d2s2(k) 
      + 0.25*(bs2g*bs2c)*PT_s2s2(k)+ 0.5*(b1c*b3nl_from_b1(b1g) + b1g*b3nl_from_b1(b1c))*PT_d1d3(k);
    P_1loop *= g4;
  }

  if (nl1 == nl2)
  {
    return cluster_bias1*cluster_B1_1*PK + P_1loop;
  }
  else
  {
    const double cluster_B1_2 = weighted_B1(nl2, z); 
    return cluster_bias1*cluster_B1_2*PK + P_1loop;
  }
}

double int_for_int_for_binned_p_cc_incl_halo_exclusion(double lnM, void* params)
{
  double* ar = (double*) params; //n_lambda, z
  const int nl1 = ar[0]; 
  const int nl2 = ar[1];
  const double k = ar[2]; 
  const double a = ar[3];
  const double lnM1 = ar[4];

  const double z = 1.0/a - 1.0;
  const double M1 = exp(lnM1);  
  const double M2 = exp(lnM);  

  const double B1_1 = B1xBselection(M1, a);
  const double B1_2 = B1xBselection(M2, a); 

  const double tmp = 0.75/(M_PI*Ntable_cluster.delta_exclusion*cosmology.rho_crit*cosmology.Omega_m);
  const double R1 = pow(M1*tmp, 1./3.);
  const double R2 = pow(M2*tmp, 1./3.);
  
  // Cosmolike This is to remove all the small oscilliating stuffs. It makes the integral faster
  const double R = (k*(0.5*(R1 + R2)) > 1) ? 0.0 : 0.5*(R1 + R2);

  const double VexclWR = 4.0*M_PI/3.*3*(sin(k*R) - k*R*cos(k*R))/(k*k*k);
  
  const double PK = (R == 0) ?  Pdelta(k,a)*B1_1*B1_2 :
    (pk_halo_with_exclusion_Rtab(k, R, a, 1., 1., 1) + VexclWR)*B1_1*B1_2 - VexclWR;

  double res1;
  {
    double params_in[2] = {nl1, z};
    res1 = dndlnM_times_binned_P_lambda_obs_given_M(lnM1, params_in)*nuisance.cluster_selection[0];
  }
  double res2;
  {
    double params_in[2] = {nl2, z};
    res2 = dndlnM_times_binned_P_lambda_obs_given_M(lnM, params_in)*nuisance.cluster_selection[0];
  }
  return PK*res1*res2;
}

double int_for_binned_p_cc_incl_halo_exclusion(double lnM, void* params)
{
  double* ar = (double*) params;
  const int nl1 = ar[0]; 
  const int nl2 = ar[1];
  const double k = ar[2]; 
  const double a = ar[3];
  const int init_static_vars_only = (int) ar[4];
  
  const double ln_M_min = log_M_min/LOG10_E; 
  const double ln_M_max = log_M_max/LOG10_E;

  double params_in[5] = {nl1, nl2, k, a, lnM};
  return (init_static_vars_only == 1) ? 
    int_for_int_for_binned_p_cc_incl_halo_exclusion(ln_M_min, (void*) params_in) :
    int_gsl_integrate_low_precision(int_for_int_for_binned_p_cc_incl_halo_exclusion, 
      (void*) params_in, ln_M_min, ln_M_max, NULL, GSL_WORKSPACE_SIZE); 
}

double binned_p_cc_incl_halo_exclusion_nointerp(double k, double a, int nl1, int nl2, 
int init_static_vars_only)
{
  const double z = 1./a-1;
  const double ln_M_min = log_M_min/LOG10_E; 
  const double ln_M_max = log_M_max/LOG10_E; 
  
  double norm1;
  {
    double params[2] = {nl1, z};
    norm1 = (init_static_vars_only == 1) ? 
      dndlnM_times_binned_P_lambda_obs_given_M(ln_M_min, (void*) params) :
      int_gsl_integrate_low_precision(dndlnM_times_binned_P_lambda_obs_given_M, (void*) params,
        ln_M_min, ln_M_max, NULL, GSL_WORKSPACE_SIZE); 
  }
  double norm2;
  { 
    double params[2] = {nl2, z};
    norm2 = (init_static_vars_only == 1) ? 
      dndlnM_times_binned_P_lambda_obs_given_M(ln_M_min, (void*) params) :
      int_gsl_integrate_low_precision(dndlnM_times_binned_P_lambda_obs_given_M, (void*) params,
        ln_M_min, ln_M_max, NULL, GSL_WORKSPACE_SIZE);  
  }

  double params[5] = {nl1, nl2, k, a, init_static_vars_only};
  const double res = (init_static_vars_only == 1) ? 
    int_for_binned_p_cc_incl_halo_exclusion(ln_M_min, (void*) params) :
    int_gsl_integrate_low_precision(int_for_binned_p_cc_incl_halo_exclusion, (void*) params,
      ln_M_min, ln_M_max, NULL, GSL_WORKSPACE_SIZE);  
  
  return (norm1*norm2 == 0) ? 0.0 : res/norm1/norm2;
}

double binned_p_cc_incl_halo_exclusion_with_constant_lambd_ninterp(double k, double a, int nl1, 
int nl2)
{  
  const double z = 1./a-1;
    
  const double B1_1 = weighted_B1(nl1, z); 
  
  const double B1_2 = (nl1 == nl2) ? B1_1 : B1_2;

  // RedMaPPer exclusion radius (Rykoff et al. 2014  eq4) in units coverH0 [change to comoving]
  const double R = 1.5*pow(0.25*(Cluster.N_min[N_lambda1] + Cluster.N_min[N_lambda2] + 
    Cluster.N_max[N_lambda1] + Cluster.N_max[N_lambda2])/100., 0.2)/cosmology.coverH0/a; 
  
  const double VexclWR = 4*M_PI*(sin(k*R) - k*R*cos(k*R))/pow(k, 3);
  
  const double cutoff =1.;
  
  double pk;
  if(R == 0.)
  {
    pk = Pdelta(k,a)*B1_1*B1_2;
  }
  if (k*R > cutoff)
  {    
    const double PK_cutoff = Pdelta(cutoff/R, a)*B1_1*B1_2;
    
    const double VexclWR_cutoff = 4*M_PI*(sin(cutoff) - cutoff*cos(cutoff))/pow(cutoff/R, 3);
    
    double Pexclusion_cutoff = (pk_halo_with_exclusion(cutoff/R, R, a, 1, 1, 1)+VexclWR_cutoff)*
      B1_1*B1_2 - VexclWRcutoff;
    
    const double k_cutoff = cutoff/R;
    
    pk = (Pdelta(k,a)*B1_1*B1_2 - VexclWR*(1 + (PK_cutoff - VexclWR_cutoff - Pexclusioncutoff)/VexclWR_cutoff))
    pk *= pow((k/kcutoff), -0.7);
  }
  else
  {
    //Cosmolike: Check it out!! This is my cool trick!! 
    pk = (pk_halo_with_exclusion(k, R, a, 1, 1, 1) + VexclWR)*B1_1*B1_2 - VexclWR; 
  }
  return pk;
}

double binned_p_cc_incl_halo_exclusion(double k, double a, int nl1, int nl2)
{
  static cosmopara C;
  static nuisancepara N;
  static double**** table = 0;

  const int N_k = Ntable_cluster.N_k_exclusion_pk_for_cell;
  const double ln_k_min = log(1E-2); 
  const double ln_k_max = log(3E6);
  const double dlnk = (ln_k_max - ln_k_min)/((double) N_k - 1.0);
  
  const int N_a = Ntable_cluster.N_a;
  const double amin = 1.0/(1.0 + tomo.cluster_zmax[tomo.cluster_Nbin - 1]);
  const double amax = 1.0/(1.0 + tomo.cluster_zmin[0]) + 0.01;
  const double da = (amax - amin)/((double) N_a - 1.0);

  if(table == 0)
  {
    table = (double****) malloc(sizeof(double***)*Cluster.N200_Nbin);
    for(int i=0; i<Cluster.N200_Nbin; i++)
    {
      table[i] = (double***) malloc(sizeof(double**)*Cluster.N200_Nbin);
      for(int j=0; j<Cluster.N200_Nbin; j++)
      {
        table[i][j] = (double**) malloc(sizeof(double*)*N_k);
        for(int k=0; k<N_k; k++)
        {
          table[i][j][k] = (double*) malloc(sizeof(double)*N_a);
        }
      }
    }
  }
  if (recompute_DESclusters(C, N))
  { 
    if(ASSUME_CONSTANT_LAMBD_HALO_EXCLUSION == 1)
    {
      binned_p_cc_incl_halo_exclusion_with_constant_lambd_ninterp(kin, aa, i, j); // init static
                                                                                  // vars only
    }
    else
    {
      binned_p_cc_incl_halo_exclusion_nointerp(exp(logkmin), amin, 0, 0, 1); // init static 
                                                                             // vars only
    }
    #pragma omp parallel for collapse(4)
    for(int i=0; i<Cluster.N200_Nbin; i++)
    {
      for(int j=i; j<Cluster.N200_Nbin; j++)
      {
        for (int k=0; k<N_a; k++) 
        { 
          for (int p=0; p<N_k; p++) 
          { 
            const double aa = (amin + k*da > 1.0) ? 1.0 : amin + k*da;
            const double kk = exp(logkmin + p*dlnk);
            if(ASSUME_CONSTANT_LAMBD_HALO_EXCLUSION == 1)
            {
              table_P_NL[i][j][k][p] = 
                binned_p_cc_incl_halo_exclusion_with_constant_lambd_ninterp(kk, aa, i, j);
            }
            else 
            {
              table_P_NL[i][j][k][p] = binned_p_cc_incl_halo_exclusion_nointerp(kk, aa, i, j, 0);
            }
            table_P_NL[i][j][k][p] = log(kin*kin*kin*sqrt(aa)*table_P_NL[i][j][k][p]);
            table_P_NL[j][i][k][p] = table_P_NL[i][j][k][p];
          }
        }
      }
    }
    update_cosmopara(&C);
    update_nuisance(&N);
  }

  if(!(a>0) || !(a<1)) 
  {
    log_fatal("a>0 and a<1 not true");
    exit(1);
  }
  if(nl1 < 0 || nl1 > Cluster.N200_Nbin -1 || nl2 < 0 || nl2 > Cluster.N200_Nbin -1)
  {
    log_fatal("invalid bin input (nl1, nl2) = (%d, %d)", nl1, nl2);
    exit(1);
  }
  const double lnk = log(k);
  const double val = interpol2d(table[nl1][nl2], N_a, amin, amax, da, a, N_k, ln_k_min, logkmax, 
    ln_k_max, dlnk, 0., 0.);
  return exp(val)/(k*k*k*sqrt(a));
}

// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// CLUSTER CROSS SPECTRUM WITH GALAXIES 
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------

// nl = lambda_obs bin, nj = galaxy redshift bin
double binned_p_cg(double k, double a, int nl, int nj, int use_linear_ps)
{ // binend in lambda_obs
  if(!(a>0) || !(a<1)) 
  {
    log_fatal("a>0 and a<1 not true");
    exit(1);
  }
  const double z = 1.0/a - 1.0;
  const double cluster_bias = weighted_B1(nl, z);
  const double PK = use_linear_ps == 1 ? p_lin(k, a) : Pdelta(k, a);
  
  const double PCG = cluster_bias * gbias.b1_function(z, nj) * PK;

  double PCG_1L = 0.;
  if ((has_b2_galaxies() || Cluster.nonlinear_bias == 1) && use_linear_ps != 1)
  {
    const double g = growfac(a);
    const double g4 = g*g*g*g;
    
    const double b1c = weighted_B1(nl, z);
    const double b2c = weighted_B2(nl, z);
    const double bs2c = bs2_from_b1(b1c);

    const double b1g = gbias.b1_function(z, nj);
    const double b2g = (gbias.b2[nj]  == 0) ? b2_from_b1(gbias.b[nj]) : gbias.b2[nj];
    const double bs2g = gbias.bs2[nj];
    
    PCG_1L = 0.5*(b1c*b2g+b2c*b1g)*PT_d1d2(k) + 0.25*b2g*b2c*PT_d2d2(k) 
      + 0.5*(b1c*bs2g+b1g*bs2c)*PT_d1s2(k) + 0.25*(b2c*bs2g+b2g*bs2c)*PT_d2s2(k) 
      + 0.25*(bs2g*bs2c)*PT_s2s2(k) + 0.5*(b1c*b3nl_from_b1(b1g) 
      + b1g*b3nl_from_b1(b1c))*PT_d1d3(k);
    PCG_1L *= g4;
  }

  // Cosmolike: this makes the code much much faster like 1000 times faster
  return (PCG + PCG_1L) < 0 ? 0.0 : (PCG + Pcg_1l)
}

// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// CLUSTER CROSS SPECTRUM WITH MATTER)
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------

double int_for_binned_p_cm(double lnM, void* params)
{ // binned in lambda_obs (nl = lambda_obs bin)
  double* ar = (double*) params;
  const int nl = (int) ar[0];
  const double z = ar[1];
  const double k = ar[2];
  const int include_1h_term = (int) ar[3];
  const int use_linear_ps = (int) ar[4];
  
  const double M = exp(lnM); 
  const double a = 1./(1 + z);

  const double PCM_1H = (include_1h_term == 1) ? 
    M/(cosmology.rho_crit*cosmology.Omega_m)*u_nfw_c(conc(M, a), k, M, a) : 0.0; 
   
  const double B1BS = B1_BS(M, a);
  const double PCM_2H = use_linear_ps == 1 ? B1BS*p_lin(k, a) : B1BS*Pdelta(k, a);
  
  double params_in[2] = {(double) nl, z}
  return (PCM_1H + PCM_2H)*dndlnM_times_binned_P_lambda_obs_given_M(lnM, (void*) params_in);
}

double binned_p_cm_nointerp(double k, double a, int nl, int include_1h_term, int use_linear_ps, 
int init_static_vars_only)
{
  if(!(a>0) || !(a<1)) 
  {
    log_fatal("a>0 and a<1 not true");
    exit(1);
  }
  const double z = 1.0/a - 1;
  double params[5] = {(double) nl, z, k, (double) include_1h_term, (double) use_linear_ps};
  const double ln_M_min = log_M_min/LOG10_E;
  const double ln_M_max = log_M_max/LOG10_E;  
  
  const double ncounts = binned_average_number_counts(nl, z);

  return (init_static_vars_only == 1) ? int_for_binned_p_cm(ln_M_min, (void*) params) :
    int_gsl_integrate_low_precision(int_for_binned_p_cm, (void*) params, ln_M_min, ln_M_max, NULL, 
      GSL_WORKSPACE_SIZE)/ncounts; 
}

double binned_p_cm(double k, double a, int nl, int ni, int include_1h_term, int use_linear_ps)
{
  static cosmopara C;
  static nuisancepara N;
  static double**** table = 0;
  static double* amin = 0;
  static double* amax = 0;
  static double* da = 0;

  const int N_k = Ntable.N_ell;
  const double ln_k_min = log(limits.k_min_cH0);
  const double ln_k_max = log(limits.k_max_cH0);
  const double dlnk = (ln_k_max - ln_k_min)/((double) N_k - 1.0);

  const int N_a = binned_p_cm_size_a_table;

  if(table == 0)
  {
    table = (double****) malloc(sizeof(double***)*Cluster.N200_Nbin);
    for(int i=0; i<Cluster.N200_Nbin; i++)
    {
      table[i] = (double***) malloc(sizeof(double**)*tomo.cluster_Nbin);
      for(int j=0; j<tomo.cluster_Nbin; j++)
      {
        table[i][j] = (double**) malloc(sizeof(double*)*N_k);
        for(int k=0; k<N_k; k++)
        {
          table[i][j][k] = (double*) malloc(sizeof(double)*N_a);
        }
      }
    }
    amin = (double*) malloc(sizeof(double)*tomo.cluster_Nbin);
    amax = (double*) malloc(sizeof(double)*tomo.cluster_Nbin);
    da = (double*) malloc(sizeof(double)*tomo.cluster_Nbin);
    for(int i=0; i<tomo.cluster_Nbin; i++)
    {
      const double zmin = tomo.cluster_zmin[j];
      const double zmax = tomo.cluster_zmax[j];
      amin[i] = 1.0/(1.0 + zmax); 
      amax[i] = 1.0/(1.0 + zmin);
      da[i] = (amax - amin)/((double) N_a - 1.0);
    }
  }
  if (recompute_clusters(C, N))
  {
    binned_p_cm_nointerp(exp(ln_k_min), amin[0], 0, include_1h_term, use_linear_ps, 1); //init static
                                                                                        //vars only
    #pragma omp parallel for collapse(4)
    for(int i=0; i<Cluster.N200_Nbin; i++)
    {
      for(int j=0; j<tomo.cluster_Nbin; j++)
      {
        for(int k=0; k<N_k; k++)
        {
          for(int p=0; p<N_a; p++)
          { 
            table[i][j][k][p] = binned_p_cm_nointerp(exp(ln_k_min + k*dlnk), amin[i] + p*da[i], nl, 
              include_1h_term, use_linear_ps, 0); 
          }
        }
      }
    }
    update_cosmopara(&C);
    update_nuisance(&N);   
  }
  if(!(a>0) || !(a<1)) 
  {
    log_fatal("a>0 and a<1 not true");
    exit(1);
  }
  if(ni < 0 || ni > tomo.cluster_Nbin -1 || nl < 0 || nl > Cluster.N200_Nbin -1)
  {
    log_fatal("invalid bin input (nl, ni) = (%d, %d)", nl, ni);
    exit(1);
  }
  const double lnk = log(k);
  return ((lnk < ln_k_min) & (lnk > ln_k_max)) ? 0.0 : interpol2d(table[nl][ni], N_k, ln_k_min, 
    ln_k_max, dlnk, lnk, N_a, amin[ni], amax[ni], da[ni], a, 1.0, 1.0); 
}

// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// Area
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------

// Cocoa: we try to avoid reading of files in the cosmolike_core code 
// Cocoa: (reading is done in the C++/python interface)
void setup_get_area(int* io_nz, double** io_z, double** io_A, int io) 
{
  static int nz;
  static double* z = NULL;
  static double* A = NULL;

  // IO == 1 IMPLES THAT IO_A(Z) WILL COPIED TO LOCAL A(Z)
  if (io == 1)
  {
    if (z != NULL)
    {
      free(z);
    }
    if (A != NULL)
    {
      free(A);
    }

    nz = (*io_nz);
    if (!(nz > 5))
    {
      log_fatal("array to small");
      exit(1);
    }

    z = (double*) malloc(nz*sizeof(double));
    A = (double*) malloc(nz*sizeof(double));
    if (z == NULL || A == NULL)
    {
      log_fatal("fail allocation");
      exit(1);
    }

    for (int i=0; i<nz; i++)
    {
      z[i] = (*io_z)[i];
      A[i] = (*io_A)[i];
    }
  }
  else
  {
    // IO != 1 IMPLES THAT LOCAL A(Z) WILL BE COPIED TO IO_A(Z)
    if (z == NULL || A == NULL)
    {
      log_fatal("array/pointer not allocated");
      exit(1);
    }

    if((*io_z) != NULL)
    {
      free((*io_z));
      (*io_z) = NULL;
    }
    else if((*io_A) != NULL)
    {
      free((*io_A));
      (*io_A) = NULL;
    }

    (*io_z) = (double*) malloc(nz*sizeof(double));
    (*io_A) = (double*) malloc(nz*sizeof(double));
    if((*io_z) == NULL || (*io_A) == NULL) 
    {
      log_fatal("fail allocation");
      exit(1);
    }

    for (int i=0; i<nz; i++)
    {
      (*io_z)[i] = z[i];
      (*io_A)[i] = A[i];
    }
  }
}

double get_area(double zz)
{
  if(INTERPOLATE_SURVEY_AREA == 1)
  {
    static gsl_spline* fA = NULL;

    if (fA == NULL)
    {
      double** z = (double**) malloc(1*sizeof(double*));
      double** A = (double**) malloc(1*sizeof(double*));
      int nz;   
      
      if (z == NULL || A == NULL)
      {
        log_fatal("fail allocation");
        exit(1);
      }
      else 
      {
        (*z) = NULL;
        (*A) = NULL;
      }

      setup_get_area(&nz, z, A, 0);

      const gsl_interp_type* T = gsl_interp_linear;
      fA = gsl_spline_alloc(T, nz);
      if (fA == NULL)
      {
        log_fatal("fail allocation");
        exit(1);
      }

      int status = 0;
      status = gsl_spline_init(fA, (*z), (*A), nz);
      if (status) 
      {
        log_fatal(gsl_strerror(status));
        exit(1);
      }
    }

    double res;
    int status = gsl_spline_eval_e(fA, zz, NULL, &res);
    if (status) 
    {
      log_fatal(gsl_strerror(status));
      exit(1);
    }
    return res;
  }
  else
  {
    return survey.area;
  }
}