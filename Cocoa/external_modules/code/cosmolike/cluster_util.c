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
static int GSL_WORKSPACE_SIZE = 250;
static double log_M_min = 12.0;
static double log_M_max = 15.9; 

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
  const double Mpiv = 5E14; //Msun/h
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
void setup_SDSS_P_true_lambda_given_mass(int* io_nintrinsic_sigma,  double** io_intrinsic_sigma, 
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
  // intrisic scatter of mass-richness relation
  const double intrinsic_sigma = nuisance.cluster_MOR[3];
  
  // average true satellite richness given mass
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

double SDSS_int_P_lambda_obs_given_M_nointerp(double true_lambda, void* params)
{
  double* ar = (double*) params;
  
  const double M = ar[0];
  const double z = ar[1];
  const double obs_lambda = ar[2];

  const double r1 = SDSS_P_lambda_obs_given_true_lambda(obs_lambda, true_lambda, z);
  const double r2 = SDSS_P_true_lambda_given_mass(true_lambda, M, z);

  return r1*r2;
}

double SDSS_P_lambda_obs_given_M_nointerp(double obs_lambda, void* params)
{
  double* ar = (double*) params; // mass , redshift 
  
  const double M = ar[0];
  const double z = ar[1];
  
  double params_in[3] = {M, z, obs_lambda}; 
  
  const double true_lambda_min =  3.;
  const double true_lambda_max =  160.;

  return int_gsl_integrate_medium_precision(SDSS_int_P_lambda_obs_given_M_nointerp, 
    (void*) params_in, true_lambda_min, true_lambda_max, NULL, GSL_WORKSPACE_SIZE);
}

// \int_(bin_lambda_obs_min)^(bin_lambda_obs_max) \dlambda_obs P(\lambda_obs|M)
// (see for example https://arxiv.org/pdf/1810.09456.pdf - eq3 qnd 6) 
double SDSS_binned_P_lambda_obs_given_M_nointerp(int nl, double M, double z)
{
  double params[2] = {M, z};
  
  const double bin_lambda_obs_min = Cluster.N_min[nl];
  const double bin_lambda_obs_max = Cluster.N_max[nl];
  
  return int_gsl_integrate_medium_precision(SDSS_P_lambda_obs_given_M_nointerp, (void*) params, 
    bin_lambda_obs_min, bin_lambda_obs_max, NULL, GSL_WORKSPACE_SIZE);
}

// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// INTERFACE
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------

// \int_(bin_lambda_obs_min)^(bin_lambda_obs_max) \dlambda_obs P(\lambda_obs|M)
// (see for example https://arxiv.org/pdf/1810.09456.pdf - eq3 qnd 6) 
double binned_P_lambda_obs_given_M_nointerp(int nl, double M, double z)
{
  if (strcmp(Cluster.model, "SDSS") == 0) 
  {
    return SDSS_binned_P_lambda_obs_given_M_nointerp(nl, M, z);
  }
  else if (strcmp(Cluster.model, "BUZZARD") == 0)
  {
    return buzzard_binned_P_lambda_obs_given_M_nointerp(nl, M, z);
  }
  else
  {
    log_fatal("option not implemented");
    exit(1);
  }
}

// Threading this function is not allowed (unless nl is fixed)
double binned_P_lambda_obs_given_M(int nl, double M, double z) 
{
  static int BIN_OBS_LAMBDA = -1; // lambda = bin number of observed lambda
  static cosmopara C;
  static nuisancepara N;
  static double** table = 0;

  const double logM = log10(M);
  const double logmmin = log_M_min;
  const double logmmax = log_M_max; 
  const double zmin = 0.20; 
  const double zmax = 0.80;
  const int nz = 10;
  const int nm = 50;
  const double dm = (logmmax-logmmin)/(nm-1.);
  const double dz = (zmax-zmin)/(nz-1.);
  const double mmin = pow(10.0, logmmin);
  const double mmax = pow(10.0, logmmax);

  if (table == 0)
  {
    table = (double**) malloc(sizeof(double*)*nz);
    for(int i=0; i<nz; i++)
    {
      table[i] = (double*) malloc(sizeof(double)*nm);
    }
  }
  if ((BIN_OBS_LAMBDA != nl) || recompute_clusters(C, N))
  { 
    {
      const int i = 0;
      {
        const int j = 0;
        const double M = pow(10.0, logmmin + j*dm);
        table[i][j] = binned_P_lambda_obs_given_M_nointerp(nl, M, z);
      } 
      #pragma omp parallel for     
      for (int j=1; j<nm; j++) 
      {
        const double M = pow(10.0, logmmin + j*dm);
        table[i][j] = binned_P_lambda_obs_given_M_nointerp(nl, M, z);
      }
    }
    #pragma omp parallel for
    for (int i=1; i<nz; i++) 
    {
      for (int j=0; j<nm; j++) 
      {
        const double M = pow(10.0, logmmin + j*dm);
        table[i][j] = binned_P_lambda_obs_given_M_nointerp(nl, M, z);
      }
    }
    
    update_cosmopara(&C);
    update_nuisance(&N);
    BIN_OBS_LAMBDA = nl;
  }
  if (z > zmin && z < zmax && M > mmin && M < mmax)
  {
    return interpol2d(table, nz, zmin, zmax, dz, z, nm, logmmin, logmmax, dm, logM, 1, 1);
  }
  return 0.;
}

// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// CLUSTER BIAS
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------



double B1_x_Bselection(double mass, double a)
{ // cluster bias including selection bias
   double B1_in;
   if (clusterAnalysisChoice.bias_model == 0) 
    {B1_in = B1(mass,a);
   else if (clusterAnalysisChoice.bias_model==1) B1_in = tinker_emu_B1_tab(mass,a);
   else {
       printf("bias model %i not implemented\n", clusterAnalysisChoice.bias_model);
       exit(0); 
   }
   
  double Bselection;



   if (nuisance.N_cluster_selection==1) return B1_in*nuisance.cluster_selection[0]; 
   if (nuisance.N_cluster_selection==2){
    return B1_in*nuisance.cluster_selection[0]*pow((mass/5E14), nuisance.cluster_selection[1]); 
  }
  else{
      printf("not implemented\n\n");
      exit(0);
  }
}



//cluster bias (b2) including selection bias
/*
 * Return : b2(m) corrected by seleciton bias depends on paramters defined in nuisance.
 *
 */
double B2xBselection(double mass, double a){
    double B1_in;
   if (clusterAnalysisChoice.bias_model==0) B1_in = B1(mass,a);
   else if (clusterAnalysisChoice.bias_model==1) B1_in = tinker_emu_B1_tab(mass,a);
   else {
       printf("bias model %i not implemented\n", clusterAnalysisChoice.bias_model);
       exit(0); 
   }
   if (nuisance.N_cluster_selection==1) return b2_from_b1(B1_in)*nuisance.cluster_selection[0]; 
   if (nuisance.N_cluster_selection==2){
    return b2_from_b1(B1_in)*nuisance.cluster_selection[0]*pow((mass/5E14), nuisance.cluster_selection[1]); 
  }
  else{
      printf("not implemented\n\n");
      exit(0);
  }
}


double int_dndlogM_times_binned_P_lambda_obs_given_M(double logM, void* params)
{
   double *ar = (double*) params; //n_lambda, z
   
   const int nl = (int) ar[0];
   const double z = ar[1];

   const double M = exp(logM);
   const double a = 1./(1. + z);
   
   return massfunc(M, a)*M*binned_P_lambda_obs_given_M(nl, M , z);
}


//
////cluster bias
double int_dndlogM_P_lambda_obs_z(double logM, void *params){
   double *array = (double*) params; //n_lambda, z
   double m = exp(logM);

   double massfunc_in;
   if (clusterAnalysisChoice.hmf_model==0) massfunc_in = massfunc(m,1./(1+array[1]));
   else if (clusterAnalysisChoice.hmf_model==1) massfunc_in = tinker_emu_massfunc_tab(m,1./(1+array[1]));
   else {
       printf("massfunc model %i not implemented\n", clusterAnalysisChoice.bias_model);
       exit(0); 
   }

   return massfunc_in*m*probability_observed_richness_given_mass_tab((int) array[0], m , array[1]);
}



double int_weighted_bias_nointerp(double logM, void* params)
{
  double *ar = (double*) params; //nl, z
  
  const double z = ar[1];
  
  const double M = exp(logM);  
  const double a = 1./(1. + z);
  
  return B1(M, a)*int_dndlogM_times_binned_P_lambda_obs_given_M(logM, params)*
    nuisance.cluster_selection[0]; 
}

double int_weightedbias(double logM, void *params){
   double * array = (double*) params; //n_lambda, z
   double mass = exp(logM);  
   double a = 1./(1+array[1]);
   return B1xBselection(mass,a)*int_dndlogM_P_lambda_obs_z(logM, params);
}
 
double int_weightedbias2(double logM, void *params){
    //modifiedchto
   double * array = (double*) params; //n_lambda, z
   double mass = exp(logM);  
   double a = 1./(1+array[1]);
    return B2xBselection(mass,a)*int_dndlogM_P_lambda_obs_z(logM, params);
}

double weighted_bias_nointerp(int nl, double z)
{
  double param[2] = {(double) nl, z};
  
  const double mmin = log_M_min/0.4342944819; // log(pow(10.,12.)) --- 0.4342944819 = log10(e)
  const double mmax = log_M_max/0.4342944819; // log(pow(10.,15.9)) --- 0.4342944819 = log10(e)
  
  const double r1 = int_gsl_integrate_low_precision(int_weighted_bias_nointerp, (void*) param, mmin, 
    mmax, NULL, GSL_WORKSPACE_SIZE);
  
  const double r2 = int_gsl_integrate_low_precision(int_dndlogM_times_binned_P_lambda_obs_given_M, 
    (void*) param, mmin, mmax, NULL, GSL_WORKSPACE_SIZE);

  return r1/r2; 
}

double weighted_bias(int nl, double z)
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
    {
      const int n = 0;
      {
        const int i = 0;
        const double a = amin + i*da;
        table[n][i] = weighted_bias_nointerp(n, 1./a - 1.);
      }
      #pragma omp parallel for
      for (int i=1; i<nasize; i++)
      {
        const double a = amin + i*da;
        table[n][i] = weighted_bias_nointerp(n, 1./a - 1.);
      }
    }
    #pragma omp parallel for
    for (int n=1; n<nlsize; n++)
    {
      for (int i=0; i<nasize; i++)
      {
        const double a = amin + i*da;
        table[n][i] = weighted_bias_nointerp(n, 1./a - 1.);
      }
    }
    update_cosmopara(&C);
    update_nuisance(&N);
  }

  const double a = 1./(z+1.);
  return interpol(table[nl], nasize, amin, amax, da, a, 0., 0.);
}

double int_weightedbias(double logM, void *params){
   double * array = (double*) params; //n_lambda, z
   double mass = exp(logM);  
   double a = 1./(1+array[1]);
   return B1xBselection(mass,a)*int_dndlogM_P_lambda_obs_z(logM, params);
}
 
double int_weightedbias2(double logM, void *params){
    //modifiedchto
   double * array = (double*) params; //n_lambda, z
   double mass = exp(logM);  
   double a = 1./(1+array[1]);
    return B2xBselection(mass,a)*int_dndlogM_P_lambda_obs_z(logM, params);
}

double weighted_bias_exact(int nlambda, double z){
   double param[2] = {1.0*nlambda, z};
    if (Cluster.N_min[nlambda]>1E10){
        double denominator = int_gsl_integrate_low_precision(int_dndlogM_P_lambda_obs_z, (void*) param,log(Cluster.N_min[nlambda]),log(Cluster.N_max[nlambda]),NULL,1000);
        if(denominator<=0) return 0;
        else return int_gsl_integrate_low_precision(int_weightedbias, (void*) param, log(Cluster.N_min[nlambda]),log(Cluster.N_max[nlambda]), NULL,1000)/denominator; 
    }
    

   double denominator = int_gsl_integrate_low_precision(int_dndlogM_P_lambda_obs_z, (void*) param,log(pow(10.,12.)),log(pow(10.,15.9)),NULL,1000);
   if(denominator<=0) return 0;
   else return int_gsl_integrate_low_precision(int_weightedbias, (void*) param, log(pow(10.,12.)),log(pow(10.,15.9)),NULL,1000)/denominator; 
}
double weighted_bias2_exact(int nlambda, double z){
   double param[2] = {1.0*nlambda, z};
   double denominator = int_gsl_integrate_low_precision(int_dndlogM_P_lambda_obs_z, (void*) param,log(pow(10.,12.)),log(pow(10.,15.9)),NULL,1000);
   if(denominator<=0) return 0;
   else return int_gsl_integrate_low_precision(int_weightedbias2, (void*) param, log(pow(10.,12.)),log(pow(10.,15.9)),NULL,1000)/denominator; 
}
double weighted_bias(int nlambda, double z){
   static cosmopara C;
   static nuisancepara N;
   static double **table;
   static double da, amin, amax;
   if (table==0) {
      table = create_double_matrix(0, Cluster.N200_Nbin-1, 0, Ntable.N_a);
      double zmin, zmax;
      zmin = fmax(tomo.cluster_zmin[0]-0.05,0.01); zmax = tomo.cluster_zmax[tomo.cluster_Nbin-1]+0.05;
      amin = 1./(1.+zmax); amax = 1./(1.+zmin);
      da = (amax-amin)/(Ntable.N_a-1.);
      //printf("zmin = %.2f  zmax = %.2f\n",zmin, zmax);
   }
   if (recompute_DESclusters(C, N)){
      for (int n = 0; n < Cluster.N200_Nbin; n++){
         double aa = amin;
         for (int i = 0; i < Ntable.N_a; i++, aa+=da){
            table[n][i] = weighted_bias_exact(n,1./aa-1.);
         }
      }
      update_cosmopara(&C);
      update_nuisance(&N);
   }
   return interpol(table[nlambda], Ntable.N_a, amin, amax, da, 1./(z+1.), 0., 0.);
}

double weighted_bias2(int nlambda, double z){
   static cosmopara C;
   static nuisancepara N;
   static double **table;
   static double da, amin, amax;
   if (table==0) {
      table = create_double_matrix(0, Cluster.N200_Nbin-1, 0, Ntable.N_a);
      double zmin, zmax;
      zmin = fmax(tomo.cluster_zmin[0]-0.05,0.01); zmax = tomo.cluster_zmax[tomo.cluster_Nbin-1]+0.05;
      amin = 1./(1.+zmax); amax = 1./(1.+zmin);
      da = (amax-amin)/(Ntable.N_a-1.);
      //printf("zmin = %.2f  zmax = %.2f\n",zmin, zmax);
   }
   if (recompute_DESclusters(C, N)){
      for (int n = 0; n < Cluster.N200_Nbin; n++){
         double aa = amin;
         for (int i = 0; i < Ntable.N_a; i++, aa+=da){
            table[n][i] = weighted_bias2_exact(n,1./aa-1.);
         }
      }
      update_cosmopara(&C);
      update_nuisance(&N);
   }
   return interpol(table[nlambda], Ntable.N_a, amin, amax, da, 1./(z+1.), 0., 0.);
}

double int_weightedbiasminus1(double logM, void *params){
   double * array = (double*) params; //n_lambda, z
   double mass = exp(logM);  
   double a = 1./(1+array[1]);
   return B1minus1xBselection(mass,a)*int_dndlogM_P_lambda_obs_z(logM, params);
}
 

double weighted_biasminus1_exact(int nlambda, double z){
   double param[2] = {1.0*nlambda, z};
    if (Cluster.N_min[nlambda]>1E10){
        double denominator = int_gsl_integrate_low_precision(int_dndlogM_P_lambda_obs_z, (void*) param,log(Cluster.N_min[nlambda]),log(Cluster.N_max[nlambda]),NULL,1000);
        if(denominator<=0) return 0;
        else return int_gsl_integrate_low_precision(int_weightedbiasminus1, (void*) param, log(Cluster.N_min[nlambda]),log(Cluster.N_max[nlambda]), NULL,1000)/denominator; 
    }
    

   double denominator = int_gsl_integrate_low_precision(int_dndlogM_P_lambda_obs_z, (void*) param,log(pow(10.,12.)),log(pow(10.,15.9)),NULL,1000);
   if(denominator<=0) return 0;
   else return int_gsl_integrate_low_precision(int_weightedbiasminus1, (void*) param, log(pow(10.,12.)),log(pow(10.,15.9)),NULL,1000)/denominator; 
}


double weighted_biasminus1(int nlambda, double z){
   static cosmopara C;
   static nuisancepara N;
   static double **table;
   static double da, amin, amax;
   if (table==0) {
      table = create_double_matrix(0, Cluster.N200_Nbin-1, 0, Ntable.N_a);
      double zmin, zmax;
      zmin = fmax(tomo.cluster_zmin[0]-0.05,0.01); zmax = tomo.cluster_zmax[tomo.cluster_Nbin-1]+0.05;
      amin = 1./(1.+zmax); amax = 1./(1.+zmin);
      da = (amax-amin)/(Ntable.N_a-1.);
      //printf("zmin = %.2f  zmax = %.2f\n",zmin, zmax);
   }
   if (recompute_DESclusters(C, N)){
      for (int n = 0; n < Cluster.N200_Nbin; n++){
         double aa = amin;
         for (int i = 0; i < Ntable.N_a; i++, aa+=da){
            table[n][i] = weighted_biasminus1_exact(n,1./aa-1.);
         }
      }
      update_cosmopara(&C);
      update_nuisance(&N);
   }
   return interpol(table[nlambda], Ntable.N_a, amin, amax, da, 1./(z+1.), 0., 0.);
}



// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// cluster mass
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------


double int_mass_mean(double logM, void *params){
   double *array = (double*) params; //n_lambda, z
   double m = exp(logM);

   double massfunc_in;
   if (clusterAnalysisChoice.hmf_model==0) massfunc_in = massfunc(m,1./(1+array[1]));
   else if (clusterAnalysisChoice.hmf_model==1) massfunc_in = tinker_emu_massfunc_tab(m,1./(1+array[1]));
   else {
       printf("massfunc model %i not implemented\n", clusterAnalysisChoice.hmf_model);
       exit(0); 
   }
   return m*massfunc_in*m*probability_observed_richness_given_mass_tab((int) array[0], m , array[1]);
}


double mass_mean(int nlambda, double z){
   double param[2] = {1.0*nlambda, z};
   double denominator = int_gsl_integrate_low_precision(int_dndlogM_P_lambda_obs_z, (void*) param,log(pow(10.,12.)),log(pow(10.,15.9)),NULL,1000);
   if(denominator<=0) return 0;
   else return int_gsl_integrate_low_precision(int_mass_mean, (void*) param, log(pow(10.,12.)),log(pow(10.,15.9)),NULL,1000)/denominator; 
}

double int_mass_dipole(double logM, void *params){
   double *array = (double*) params; //n_lambda, z
   double m = exp(logM);

   double massfunc_in;
   if (clusterAnalysisChoice.hmf_model==0) massfunc_in = massfunc(m,1./(1+array[1]));
   else if (clusterAnalysisChoice.hmf_model==1) massfunc_in = tinker_emu_massfunc_tab(m,1./(1+array[1]));
   else {
       printf("massfunc model %i not implemented\n", clusterAnalysisChoice.hmf_model);
       exit(0); 
   }
   return m*m*massfunc_in*m*probability_observed_richness_given_mass_tab((int) array[0], m , array[1]);
}

double mass_dipole(int nlambda, double z){
   double param[2] = {1.0*nlambda, z};
   double denominator = int_gsl_integrate_low_precision(int_dndlogM_P_lambda_obs_z, (void*) param,log(pow(10.,12.)),log(pow(10.,15.9)),NULL,1000);
   if(denominator<=0) return 0;
   else return int_gsl_integrate_low_precision(int_mass_dipole, (void*) param, log(pow(10.,12.)),log(pow(10.,15.9)),NULL,1000)/denominator; 
}

double int_mass_z(int nlambda, double z){
   double param[2] = {1.0*nlambda, z};
   return int_gsl_integrate_low_precision(int_mass_mean, (void*) param, log(pow(10.,12.)),log(pow(10.,15.9)),NULL,1000);
}
double int_total_mass_lambda_lambda_mass_z(double a, void* params){
   double *array = (double*) params;
   double z = 1./a-1 ;
   double norm = get_area(z);
   //printf("%f, %f\n", z, norm);
   return pow(f_K(chi(a)),2.0)*dchi_da(a)*int_mass_z((int)array[0], z)*norm;// number(deg^2/rad^2)
}
double mean_mass_Delta_lambda_obs(int N_lambda, int nz_cluster){
    double param[1] = {1.0*N_lambda};
    return 4.*M_PI/41253.0*int_gsl_integrate_low_precision(int_total_mass_lambda_lambda_mass_z, (void*) param, 1./(1+tomo.cluster_zmax[nz_cluster]),1./(1+tomo.cluster_zmin[nz_cluster]),NULL,1000)/ N_Delta_lambda_obs(N_lambda, nz_cluster);
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
  const double mmin = 12.0/0.4342944819; // log(pow(10.,12.)) --- 0.4342944819 = log10(e)
  const double mmax = 15.9/0.4342944819; // log(pow(10.,15.9)) --- 0.4342944819 = log10(e)
  
  return int_gsl_integrate_low_precision(int_dndlogM_times_binned_P_lambda_obs_given_M, 
    (void*) param, mmin, mmax, NULL, GSL_WORKSPACE_SIZE);
}

// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// CLUSTER POWER SPECTRUM AND CROSS SPECTRUM (WITH GALAXIES AND MATTER)
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
  const double z = 1./a - 1;
  const double cluster_bias1 = weighted_bias(nl1, z); 
  const double PK = use_linear_ps == 1 ? p_lin(k, a) : Pdelta(k, a);

  if (nl1 == nl2)
  {
    return cluster_bias1*cluster_bias1*PK;
  }
  else
  {
    const double cluster_bias2 = weighted_bias(nl2, z); 
    return cluster_bias1*cluster_bias2*PK;
  }
}

//Halo exclusion routine//
 
double int_P_cluster_x_cluster_clustering_exclusion_exact_2(double logM2, void *params){
       double pk;
       double * array = (double*) params; //n_lambda, z
       double mass_1 = exp(array[4]);  
       double mass_2 = exp(logM2);  
       double a = array[3];
       double k = array[2]; 
       double param1[2] = {array[0], 1./a-1};
       double param2[2] = {array[1], 1./a-1};

       double b1 = B1xBselection(mass_1,a);
       double b2 = B1xBselection(mass_2,a); 
       double R1 = pow((3*mass_1/(M_PI*4*(Ntable_cluster.delta_exclusion*cosmology.rho_crit*cosmology.Omega_m))), (1./3.));
       double R2 = pow((3*mass_2/(M_PI*4*(Ntable_cluster.delta_exclusion*cosmology.rho_crit*cosmology.Omega_m))), (1./3.));;
       double R=0.5*(R1+R2);
       if (k*R>1.0) R=0; // This is to remove all the small oscilliating stuffs. It makes the integral faster.
       double VexclWR = 4*M_PI/3.*3*(sin(k*R) - k*R*cos(k*R))/pow(k, 3.);
       if(R==0) VexclWR=0;
       if(R==0){
            pk = Pdelta(k,a)*b1*b2;
       }else{
            pk = (pk_halo_with_exclusion_Rtab(k, R, a, 1., 1., 1)+VexclWR)*b1*b2-VexclWR;
       }
       return  pk*int_dndlogM_P_lambda_obs_z(array[4], param1)*nuisance.cluster_selection[0]*int_dndlogM_P_lambda_obs_z(logM2, param2)*nuisance.cluster_selection[0];
}

double int_P_cluster_x_cluster_clustering_exclusion_exact_1(double logM1, void *params){
       double * array = (double*) params; //n_lambda, z
       double param[5] = {array[0], array[1], array[2], array[3], logM1};
       return int_gsl_integrate_low_precision(int_P_cluster_x_cluster_clustering_exclusion_exact_2, (void*) param, log(pow(10.,12.)),log(pow(10.,15.9)),NULL,1000);
} 

double P_cluster_x_cluster_clustering_exclusion_exact(double k, double a, int N_lambda1, int N_lambda2){
       double z = 1./a-1;
       double param[4] = {N_lambda1, N_lambda2, k, a};
       double param1[2] = {N_lambda1, z};
       double param2[2] = {N_lambda2, z};
       double noormalization1 = int_gsl_integrate_low_precision(int_dndlogM_P_lambda_obs_z, (void*) param1,log(pow(10.,12.)),log(pow(10.,15.9)),NULL,1000); 
       double noormalization2 = int_gsl_integrate_low_precision(int_dndlogM_P_lambda_obs_z, (void*) param2,log(pow(10.,12.)),log(pow(10.,15.9)),NULL,1000); 
       if(noormalization1*noormalization2==0) return 0.0;
       double P_unnormalized = int_gsl_integrate_low_precision(int_P_cluster_x_cluster_clustering_exclusion_exact_1, (void*) param, log(pow(10.,12.)),log(pow(10.,15.9)),NULL,1000);
       return P_unnormalized/noormalization1/noormalization2; 
}

double P_cluster_x_cluster_clustering_exclusion_tab(double k, double a, int N_lambda1, int N_lambda2){
          static cosmopara C;
          static nuisancepara N;
          static int N_lambda1_in=-1;
          static int N_lambda2_in=-1;
          static double logkmin = 0., logkmax = 0., dk = 0., da = 0.;
          static double **table_P_NL=0;
          const double amin = 1./(1+tomo.cluster_zmax[tomo.cluster_Nbin-1]);
          const double amax = 1./(1+tomo.cluster_zmin[0])+0.01;
          double klog,val;
          int i,j;
          double kin;
          logkmin = log(1E-2);
          logkmax = log(3E6);
          dk = (logkmax - logkmin)/(Ntable_cluster.N_k_exclusion_pk_for_cell-1.);
          da = (amax - amin)/(Ntable_cluster.N_a-1.);

          if (recompute_DESclusters(C, N)|| (N_lambda1_in != N_lambda1)|| (N_lambda2_in != N_lambda2)){
            update_cosmopara(&C);
            update_nuisance(&N);
            N_lambda1_in = N_lambda1;
            N_lambda2_in = N_lambda2;
            if (table_P_NL!=0) free_double_matrix(table_P_NL,0, Ntable_cluster.N_a-1, 0, Ntable_cluster.N_k_exclusion_pk_for_cell-1);
            table_P_NL = create_double_matrix(0, Ntable_cluster.N_a-1, 0, Ntable_cluster.N_k_exclusion_pk_for_cell-1);     
            double aa = amin;
            for (i=0; i<Ntable_cluster.N_a; i++, aa +=da) { 
                printf("%d, %d, \n", i, Ntable_cluster.N_a);
                for (j=0; j<Ntable_cluster.N_k_exclusion_pk_for_cell; ++j) { 
                    //printf("in %d, %d, \n", j, Ntable_cluster.N_k_exclusion_pk_for_cell);
                    if(aa>1.0) aa=1.0;
                    kin   = exp(logkmin+j*dk);
                    table_P_NL[i][j] = log(pow(kin,3)*pow(aa,0.5)*P_cluster_x_cluster_clustering_exclusion_exact(kin, aa, N_lambda1, N_lambda2));
                }
            }

          }
          klog = log(k);
          val = interpol2d(table_P_NL, Ntable_cluster.N_a, amin, amax, da, a, Ntable_cluster.N_k_exclusion_pk_for_cell, logkmin, logkmax, dk, klog, 0.0, 0.0);
          return exp(val)/k/k/k*pow(a,-0.5);
          //return P_cluster_x_cluster_clustering_mass_given_Dlambda_obs(k, a, N_lambda1, N_lambda2); 
        //exp(val)/k/k/k;
        //exp(val);
}



double P_cluster_x_cluster_clustering_exclusion_constant_lambd_exact(double k, double a, int N_lambda1, int N_lambda2){
       double z = 1./a-1;
       double cluster_bias1, cluster_bias2;
       double pk, R, VexclWR;
       //static int count =1;
       cluster_bias1 = weighted_bias(N_lambda1, z); 
       if (N_lambda1==N_lambda2){
            cluster_bias2 = cluster_bias1;
       }
       else{
          cluster_bias2 = weighted_bias(N_lambda2, z); 
       }
       R=1.5*pow(0.25*(Cluster.N_min[N_lambda1]+Cluster.N_min[N_lambda2]+Cluster.N_max[N_lambda1]+Cluster.N_max[N_lambda2])/100., 0.2)/cosmology.coverH0/a; // RedMaPPer exclusion radius (Rykoff et al. 2014  eq4) in units coverH0 [change to comoving]
       //R = pow((3*mass_mean(N_lambda1, z)/(M_PI*4*(200*cosmology.rho_crit*cosmology.Omega_m))), (1./3.));;
       //R = pow((3*mass_mean(N_lambda1, z)/(M_PI*4*(30*cosmology.rho_crit*cosmology.Omega_m))), (1./3.))/a;
       //printf("R_compare: R_200 %e, Rperco %e \n", R2, R);
       VexclWR = 4*M_PI*(sin(k*R) - k*R*cos(k*R))/pow(k, 3.);
       double cutoff =1.;
       //if (k*R>cutoff){ // to avoid high k oscillation. 1/x**2 damping term is motivated by the envelop of 3(sinkx-kx*coskx)/kx**3
       //     pk = pow(cutoff/R,2)/pow(k,2); 
       //     k = cutoff/R;
       //}
       if (k*R> cutoff){
            pk = Pdelta(k,a)*cluster_bias1*cluster_bias2;
            double Pdeltacutoff = Pdelta(cutoff/R,a)*cluster_bias1*cluster_bias2;
            double VexclWRcutoff = (4*M_PI*(sin(cutoff) - cutoff*cos(cutoff))/pow(cutoff/R, 3.));
            double Pexclusioncutoff = (pk_halo_with_exclusion(cutoff/R, R, a, 1, 1, 1)+VexclWRcutoff)*cluster_bias1*cluster_bias2-VexclWRcutoff;
            double kcutoff = cutoff/R;
            //pk = Pexclusioncutoff*(pk-VexclWR)/(Pdeltacutoff-VexclWRcutoff)*(VexclWR)/VexclWRcutoff;
            pk = (pk-VexclWR-VexclWR*(-1*Pexclusioncutoff+(Pdeltacutoff-VexclWRcutoff))/VexclWRcutoff)*pow((k/kcutoff),-0.7);
            return pk;
       }

       if(R==0){
            pk = Pdelta(k,a)*cluster_bias1*cluster_bias2;
       }else{
            pk = (pk_halo_with_exclusion(k, R, a, 1, 1, 1)+VexclWR)*cluster_bias1*cluster_bias2-VexclWR; // Check it out!! This is my cool trick!! 
       }
       return pk;
}

double P_cluster_x_cluster_clustering_exclusion_constant_lambd_tab(double k, double a, int N_lambda1, int N_lambda2, double linear){
          if(linear>0) return  P_cluster_x_cluster_clustering_mass_given_Dlambda_obs(k, a, N_lambda1, N_lambda2, linear);
          static cosmopara C;
          static nuisancepara N;
          static int N_lambda1_in=-1;
          static int N_lambda2_in=-1;
          static double logkmin = 0., logkmax = 0., dk = 0., da = 0.;
          static double **table_P_NL=0;
          const double amin = 1./(1+tomo.cluster_zmax[tomo.cluster_Nbin-1]);
          const double amax = 1./(1+tomo.cluster_zmin[0]-1E-6);
          double klog,val;
          int i,j;
          double kin;
          logkmin = log(1E-2);
          logkmax = log(1E8);
          dk = (logkmax - logkmin)/(Ntable_cluster.N_k_exclusion_pk_for_cell-1.);
          da = (amax - amin)/(Ntable_cluster.N_a-1.);

          if (recompute_DESclusters(C, N)|| (N_lambda1_in != N_lambda1)|| (N_lambda2_in != N_lambda2)){
            update_cosmopara(&C);
            update_nuisance(&N);
            N_lambda1_in = N_lambda1;
            N_lambda2_in = N_lambda2;
            if (table_P_NL!=0) free_double_matrix(table_P_NL,0, Ntable_cluster.N_a-1, 0, Ntable_cluster.N_k_exclusion_pk_for_cell-1);
            table_P_NL = create_double_matrix(0, Ntable_cluster.N_a-1, 0, Ntable_cluster.N_k_exclusion_pk_for_cell-1);     
            double aa = amin;
            for (i=0; i<Ntable_cluster.N_a; i++, aa +=da) { 
                for (j=0; j<Ntable_cluster.N_k_exclusion_pk_for_cell; ++j) { 
                    if(aa>1.0) aa=1.0;
                    kin   = exp(logkmin+j*dk);
                    table_P_NL[i][j] = log(pow(kin,3)*pow(aa, 0.5)*P_cluster_x_cluster_clustering_exclusion_constant_lambd_exact(kin, aa, N_lambda1, N_lambda2)+1E8);
                }
            }

          }
          klog = log(k);
          val = interpol2d(table_P_NL, Ntable_cluster.N_a, amin, amax, da, a, Ntable_cluster.N_k_exclusion_pk_for_cell, logkmin, logkmax, dk, klog, 1.0, 1.0);
          return (exp(val)-1E8)/k/k/k*pow(a, -0.5);
}


double P_cluster_x_cluster_clustering_mass_given_Dlambda_obs(double k, double a, int N_lambda1, int N_lambda2, double linear){
       double z = 1./a-1;
       //static int count =1;
       double cluster_bias1 = weighted_bias(N_lambda1, z); 
       double P_1loop=0;
       if ((gbias.b2[0] || (clusterAnalysisChoice.nonlinear_bias>0))&(linear<1)){
            double g4 = pow(growfac(a)/growfac(1.0),4.);
            double b1g = cluster_bias1;
            double b2g = weighted_bias2((int)N_lambda1, 1./a-1.);
            double bs2gminus1plus1 = weighted_biasminus1((int)N_lambda1, 1./a-1.)+1;
            double b1c =  weighted_bias(N_lambda2, z);
            double b2c = weighted_bias2((int)N_lambda2, 1./a-1.);
            double bs2cminus1plus1 = weighted_biasminus1((int)N_lambda2, 1./a-1.)+1;
            double bs2g = bs2_from_b1(bs2gminus1plus1);
            double bs2c = bs2_from_b1(bs2cminus1plus1);
            P_1loop = 0.5*(b1c*b2g+b2c*b1g)*PT_d1d2(k) + 0.25*b2g*b2c*PT_d2d2(k) + 0.5*(b1c*bs2g+b1g*bs2c)*PT_d1s2(k)
                +0.25*(b2c*bs2g+b2g*bs2c)*PT_d2s2(k) + 0.25*(bs2g*bs2c)*PT_s2s2(k)
                +0.5*(b1c*b3nl_from_b1(bs2gminus1plus1) + b1g*b3nl_from_b1(bs2cminus1plus1))*PT_d1d3(k);
            P_1loop *= g4;
        }

       if (N_lambda1==N_lambda2){
         if(linear>0.1) return pow(cluster_bias1,2.0)*p_lin(k,a);
         return pow(cluster_bias1,2.0)*Pdelta(k,a)+P_1loop;
       }
       else{
          double cluster_bias2 = weighted_bias(N_lambda2, z); 
          if(linear>0) return cluster_bias1*cluster_bias2*p_lin(k,a);
          return cluster_bias1*cluster_bias2*Pdelta(k,a)+P_1loop;
       }
}


// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// CLUSTER CROSS SPECTRUM (WITH GALAXIES AND MATTER)
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
  const double cluster_bias = weighted_bias(nl, z);
  const double PK = use_linear_ps == 1 ? p_lin(k, a) : Pdelta(k, a);
  return cluster_bias * gbias.b1_function(z, nj)* PK;
}

// nl = lambda_obs bin, ni = cluster redshift bin, nj = source redshift bin
double int_binned_p_cm(double logM, void* params)
{ // binned in lambda_obs (nl = lambda_obs bin)
  double* ar = (double*) params;
   
  const int nl = (int) ar[0];
  const double z = ar[1];
  const double k = ar[2];
  const int use_linear_ps = (int) ar[3];
  
  const double M = exp(logM); 
  const double a = 1./(1 + z);

  const double PCM_1h = M/(cosmology.rho_crit*cosmology.Omega_m)*u_nfw_c(conc(M, a), k, M, a); 
   
  const double PCM_2h = use_linear_ps == 1 ? B1(M, a)*p_lin(k, a) : B1(M, a)*Pdelta(k, a);
   
  return M*massfunc(M, a)*(PCM_1h + PCM_2h)*binned_P_lambda_obs_given_M(nl, M , z);
}

double int_P_cluster_mass_given_Dlambda_obs_tab(double logM, void *params){
   double *array = (double*) params; //n_lambda, z,k
   static double Nlambda=-1, z, k;
   static cosmopara C;
   static nuisancepara N;
   static double **table_P1= 0;
   static double dm =0.,logmmin = 12.0,logmmax = 15.9;
   static int NM = 50;
   if (table_P1 ==0){
      table_P1 = create_double_matrix(0,1, 0, NM-1);
      dm = (logmmax-logmmin)/(NM-1.);
   }
   if( ((Nlambda != array[0])|| (z != array[1]) || (k != array[2]) )|| recompute_DESclusters(C,N) ){
        update_cosmopara(&C);
        update_nuisance(&N);
         Nlambda = array[0];
         z = array[1];
         k = array[2];
         int j = 0;
         for (double lgM = logmmin; lgM < logmmax; lgM+= dm){
                double mass = pow(10.,lgM);
                table_P1[0][j] = int_P_cluster_mass_given_Dlambda_obs(log(mass), params);
                j++;
         }
    }
    double mass = exp(logM);  
    if (mass <  pow(10.,logmmax) && (mass > pow(10., logmmin))){
      return interpol(table_P1[0], NM, logmmin, logmmax, dm, log10(mass), 1.0, 1.0);
    }
    else return 0;
}

double int_P_cluster_mass_given_Dlambda_obs_tab_2haloonly(double logM, void *params){
   double *array = (double*) params; //n_lambda, z,k
   static double Nlambda=-1, z;
   static cosmopara C;
   static nuisancepara N;
   static double **table_P1= 0;
   static double dm =0.,logmmin = 12.0,logmmax = 15.9;
   static int NM = 50;
   if (table_P1 ==0){
      table_P1 = create_double_matrix(0,1, 0, NM-1);
      dm = (logmmax-logmmin)/(NM-1.);
   }
   if( ((Nlambda != array[0])|| (z != array[1]))|| recompute_DESclusters(C,N) ){
        update_cosmopara(&C);
        update_nuisance(&N);
         Nlambda = array[0];
         z = array[1];
         int j = 0;
         for (double lgM = logmmin; lgM < logmmax; lgM+= dm){
                double mass = pow(10.,lgM);
                table_P1[0][j] = int_P_cluster_mass_given_Dlambda_obs_2halo_nopdeltak(log(mass), params);
                j++;
         }
    }
    double a = 1./(1+array[1]);
    double mass = exp(logM);  
    if (mass <  pow(10.,logmmax) && (mass > pow(10., logmmin))){
      return interpol(table_P1[0], NM, logmmin, logmmax, dm, log10(mass), 1.0, 1.0);
    }
    else return 0;
}



double binned_p_cm(double k, double a, int nl, int use_linear_ps)
{
  if(!(a>0) || !(a<1)) 
  {
    log_fatal("a>0 and a<1 not true");
    exit(1);
  }
  const double z = 1.0/a - 1;
  double params[4] = {(double) nl, z, k, (double) use_linear_ps};
  const double mmin = log_M_min/0.4342944819; // log(pow(10.,12.)) --- 0.4342944819 = log10(e)
  const double mmax = log_M_max/0.4342944819; // log(pow(10.,15.9)) --- 0.4342944819 = log10(e)  
  const double ncounts = binned_average_number_counts(nl, z);

  return int_gsl_integrate_low_precision(int_binned_p_cm, params, mmin, mmax, NULL, 
    GSL_WORKSPACE_SIZE)/ncounts; 
}

double P_cluster_mass_given_Dlambda_obs(double k, double a, int N_lambda, int nz_cluster, int nz_galaxy){
   double z = 1./a-1;
   double params[3] = {1.0*N_lambda, z, k};
   double result; 
   if (clusterAnalysisChoice.Ytransform>0) result = int_gsl_integrate_low_precision(int_P_cluster_mass_given_Dlambda_obs_tab_2haloonly, params,log(pow(10.,12.)),log(pow(10.,15.9)),NULL,1000)/n_lambda_obs_z_tab(N_lambda, z)*Pdelta(k,a); 
   else result = int_gsl_integrate_low_precision(int_P_cluster_mass_given_Dlambda_obs_tab,params,log(pow(10.,12.)),log(pow(10.,15.9)),NULL,1000)/n_lambda_obs_z_tab(N_lambda, z); 
   //double result = int_gsl_integrate_low_precision(int_P_cluster_mass_given_Dlambda_obs_1halo,params,log(pow(10.,12.)),log(pow(10.,15.9)),NULL,1000)/n_lambda_obs_z_tab(N_lambda, z)+
   //                int_gsl_integrate_low_precision(int_P_cluster_mass_given_Dlambda_obs_2halo,params,log(pow(10.,12.)),log(pow(10.,15.9)),NULL,1000)/n_lambda_obs_z_tab(N_lambda, z)     ; 
   if(isnan(result)) return 0;
   else return result;
}


double P_cluster_mass_given_Dlambda_obs_tab(double k, double a, int N_lambda, int nz_cluster, int nz_galaxy){
    static cosmopara C;
    static nuisancepara N;
    static double dk = 0., da = 0.;
    static int N_lambda_in, nz_cluster_in, nz_galaxy_in;
    static double logkmin, logkmax;
    static double **table_P_k_a=0;
    double klog,val, result; 
    double zmin, zmax, aa, kin;
    static double amin, amax;
    int N_a = 30;
    int i, j;
    if (recompute_DESclusters(C, N) || (N_lambda_in != N_lambda) || (nz_cluster_in != nz_cluster) || (nz_galaxy_in != nz_galaxy)){
        update_cosmopara(&C);
        update_nuisance(&N);
        N_lambda_in = N_lambda; 
        nz_cluster_in = nz_cluster; 
        nz_galaxy_in = nz_galaxy;
        zmin = tomo.cluster_zmin[nz_cluster];
        zmax = tomo.cluster_zmax[nz_cluster];
        amin = 1./(1.+zmax); amax = 1./(1.+zmin);
        da = (amax-amin)/(N_a-1.);

        if (table_P_k_a!=0) free_double_matrix(table_P_k_a,0, Ntable.N_ell-1, 0, N_a-1);
        table_P_k_a = create_double_matrix(0, Ntable.N_ell-1, 0, N_a-1);     
        logkmin = log(limits.k_min_cH0);
        logkmax = log(limits.k_max_cH0);
        dk = (logkmax - logkmin)/(Ntable.N_ell-1);     
        aa = amin;
        for (j=0; j<N_a; ++j, aa+=da) { 
            for (i=0; i<Ntable.N_ell; i++) { 
                kin   = exp(logkmin+i*dk);
                result = P_cluster_mass_given_Dlambda_obs(kin, aa, N_lambda, nz_cluster, nz_galaxy); 
                //if (result ==0) table_P_k_a[i][j] = -1E10; 
                table_P_k_a[i][j] = result; 
            }
        }
        
    }
    klog = log(k);
    if ((klog > logkmin)&(klog < logkmax)) val = interpol2d(table_P_k_a, Ntable.N_ell, logkmin, logkmax, dk, klog, N_a, amin, amax, da, a, 1.0, 1.0);
    else val=0;
    return val;
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