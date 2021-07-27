#include <assert.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../cfftlog/cfftlog.h"

#include "bias.h"
#include "basics.h"
#include "cosmo3D.h"
#include "cluster_util.h"
#include "cosmo2D_cluster.h"
#include "radial_weights.h"
#include "recompute.h"
#include "redshift_spline.h"
#include "structs.h"

#include "log.c/src/log.h"

static int GSL_WORKSPACE_SIZE = 250;
static int use_linear_ps_limber = 0; /* 0 or 1 */
static double w_l_min = 0.0001;
static double w_l_max = 5.0e6; 
static INCLUDE_MAG_IN_C_CC_NONLIMBER = 0; /* 0 or 1 */
static INCLUDE_MAG_IN_C_CG_NONLIMBER = 0; /* 0 or 1 */
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// Correlation Functions (real Space) - Full Sky - bin average
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------

double w_gammat_cluster_tomo(int nt, int nl, int ni, int nj, int limber)
{ // nt = theta bin, nl = lambda_obs bin, ni = cluster redshift bin, nj = source redshift bin
  if (like.Ntheta == 0)
  {
    log_fatal("like.Ntheta not initialized");
    exit(1);
  }
  if(nt > like.Ntheta - 1)
  {
    log_fatal("error in selecting bin number");
    exit(1); 
  }

  static cosmopara C;
  static nuisancepara N;
  static galpara G;
  static double** Pl = 0;
  static double* w_vec = 0;
  
  const int ngammat_size = tomo.cgl_Npowerspectra;
  const int NSIZE = nlsize*ngammat_size;
  const int nell = limits.LMAX;
  const int ntheta = like.Ntheta;

  if (Pl == 0)
  {    
    Pl = (double**) malloc(sizeof(double*)*ntheta);
    for(int i=0; i<ntheta; i++) 
    {
      Pl[i] = (double*) malloc(sizeof(double)*nell);
    }
    w_vec = (double*) malloc(sizeof(double)*NSIZE*ntheta); 
    double xmin[ntheta];
    double xmax[ntheta];
    // Cocoa: dont thread (init of static variables inside set_bin_average)
    for (int i=0; i<ntheta; i++)
    {
      bin_avg r = set_bin_average(i,0);
      xmin[i] = r.xmin;
      xmax[i] = r.xmax;
    }
    #pragma omp parallel for
    for (int i=0; i<ntheta; i++)
    {
      double* Pmin = (double*) malloc(sizeof(double)*(nell + 1));
      double* Pmax = (double*) malloc(sizeof(double)*(nell + 1));
      for (int l=0; l<nell; l++)
      {
        bin_avg r = set_bin_average(i,l);
        Pmin[l] = r.Pmin;
        Pmax[l] = r.Pmax;
      }
      for (int l=1; l<nell; l++)
      {
        Pl[i][l] = (2.*l+1)/(4.*M_PI*l*(l+1)*(xmin[i]-xmax[i]))
          *((l+2./(2*l+1.))*(Pmin[l-1]-Pmax[l-1])
          +(2-l)*(xmin[i]*Pmin[l]-xmax[i]*Pmax[l])
          -2./(2*l+1.)*(Pmin[l+1]-Pmax[l+1]));
      }
      free(Pmin);
      free(Pmax);
    }
  }
  
  if (recompute_cs(C, G, N))
  {
    double** Cl = malloc(NSIZE*sizeof(double*));
    for (int i = 0; i<NSIZE; i++)
    {
      Cl[i] = calloc(nell, sizeof(double));
    }
    if(limber == 1)
    {
      for(int i=0; i<nlsize; i++)
      { // loop in lambda_obs - (cannot thread loops in lambda_obs  - see cluster_utils)
        C_cs_tomo_limber(limits.P_2_s_min + 1, i, ZCL(0), ZCS(0)); // init function
        #pragma omp parallel for collapse(2)
        for(int j=0; j<ngammat_size; j++)
        { 
          for (int l=2; l<nell; l++)
          {
            const int ZC = ZCL(j);
            const int ZSC = ZCS(j);
            const int q = i*ngammat_size + j;
            Cl[q][l] = (l > limits.P_2_s_min) ? C_cs_tomo_limber(l, i, ZC, ZSC) :
              C_cs_tomo_limber_nointerp(l, i, ZC, ZSC, use_linear_ps_limber);
          }
        }    
      } 
    }
    else
    { 
      log_fatal("NonLimber not implemented");
      exit(1);     
      /*
      const int L = 1;
      const double tolerance = 0.0075;    // required fractional accuracy in C(l)
      const double dev = 10. * tolerance; // will be diff exact vs Limber init to large
                                          // value in order to start while loop
      for(int j=0; j<ngammat_size; j++)
      { // Cocoa: no threading allowed here - (fftw allocation)
        const int ZC = ZCL(j);
        const int ZSC = ZCS(j);
        const int q = i*ngammat_size + j;
        Cl[q][l] = 0.0; // TODO miss nonlimber function...
      }
      #pragma omp parallel for collapse(2)
      for(int j=0; j<ngammat_size; j++)
      { 
        for (int l=limits.LMAX_NOLIMBER+1; l<nell; l++)
        {
          const int ZC = ZCL(j);
          const int ZSC = ZCS(j);
          const int q = i*ngammat_size + j;
          Cl[q][l] = C_cs_tomo_limber(l, i, ZC, ZSC);
        }
      }
      */ 
    }
    #pragma omp parallel for collapse(3)
    for(int i=0; i<nlsize; i++)
    {
      for(int j=0; j<ngammat_size; j++)
      { 
        for (int p=0; p<ntheta; p++)
        {
          const int nz = i*ngammat_size + j;
          const int q = nz*ntheta + p;
          w_vec[q] = 0;
          for (int l=1; l<nell; l++)
          {
            w_vec[q] += Pl[p][l]*Cl[nz][l];
          }
        }
      }
    }
    for (int i=0; i<NSIZE; i++)
    {
      free(Cl[i]);
    }
    free(Cl);

    update_cosmopara(&C);
    update_galpara(&G);
    update_nuisance(&N);
  }

  if (!test_zoverlap_c(ni, nj)) 
  {
    return 0.0;
  }

  const int q = (nl*ngammat_size + N_cgl(ni, nj))*ntheta + nt;
  if(q > NSIZE*ntheta - 1)
  {
    log_fatal("error in selecting bin number");
    exit(1);
  }

  return w_vec[q];
}


double w_cc_tomo(int nt, int nl1, int nl2, int ni, int nj, int limber)
{ // nt = theta bin , nl{1,2} = lambda_obs bins, n{i,j} = cluster redshift bins
  if (like.Ntheta == 0)
  {
    log_fatal("like.Ntheta not initialized");
    exit(1);
  }
  if(nt > like.Ntheta - 1)
  {
    log_fatal("error in selecting bin number");
    exit(1); 
  } 

  static cosmopara C;
  static nuisancepara N;
  static double** Pl = 0;
  static double* w_vec = 0;

  const int nlsize = Cluster.N200_Nbin;
  const int nccl_size = tomo.cluster_Nbin; // cross redshift bin not supported so not using
                                           // tomo.cc_clustering_Npowerspectra
  const int NSIZE = nlsize*nlsize*nccl_size;
  const int nell = limits.LMAX;
  const int ntheta = like.Ntheta;

  if (Pl == 0)
  {
    Pl = (double**) malloc(sizeof(double*)*ntheta);
    for(int i=0; i<ntheta; i++) 
    {
      Pl[i] = (double*) malloc(sizeof(double)*nell);
    }
    w_vec = (double*) malloc(sizeof(double)*NSIZE*ntheta);
    double xmin[ntheta];
    double xmax[ntheta];
    // Cocoa: dont thread (init of static variables inside set_bin_average)
    for (int i=0; i<ntheta; i ++)
    {
      bin_avg r = set_bin_average(i,0);
      xmin[i] = r.xmin;
      xmax[i] = r.xmax;
    }
    double** Pmin = (double**) malloc(sizeof(double)*ntheta);
    double** Pmax = (double**) malloc(sizeof(double)*ntheta);
    for (int i=0; i<ntheta; i ++)
    {
      Pmin[i] = (double*) malloc(sizeof(double)*(nell + 1));
      Pmax[i] = (double*) malloc(sizeof(double)*(nell + 1));
    }
    #pragma omp parallel for
    for (int i=0; i<ntheta; i ++)
    {
      for (int l=0; l<nell; l++)
      {
        bin_avg r = set_bin_average(i,l);
        Pmin[i][l] = r.Pmin;
        Pmax[i][l] = r.Pmax;
      }
      Pl[i][0] = 1.0;
      const double tmp = (1.0/(xmin[i] - xmax[i]))*(1. / (4. * M_PI));
      for (int l=1; l<nell; l++)
      {
        Pl[i][l] = tmp*(Pmin[i][l + 1] - Pmax[i][l + 1] - Pmin[i][l - 1] + Pmax[i][l - 1]);
      }
    }
    for (int i=0; i<ntheta; i ++)
    {
      free(Pmin[i]);
      free(Pmax[i]);
    }
    free(Pmin);
    free(Pmax);
  }
  if (recompute_cc(C, N))
  {
    double** Cl = malloc(NSIZE*sizeof(double*));
    for (int i=0; i<NSIZE; i++)
    {
      Cl[i] = calloc(nell, sizeof(double));
    }
    for (int i=0; i<nlsize; i++) 
    { // loop in lambda_obs - (cannot thread loops in lambda_obs - see cluster_utils)
      for (int j=i; j<nlsize; j++)  
      { // loop in lambda_obs - (cannot thread loops in lambda_obs - see cluster_utils)
        if(limber == 1)
        {
          C_cc_tomo_limber(limits.P_2_s_min+1, i, j, 0, 0); // need to init that
          #pragma omp parallel for collapse(2)
          for (int k=0; k<nccl_size; k++)  
          {
            for (int l=1; l<nell; l++)
            {
              const int q = i*nlsize*nccl_size + j*nccl_size + k;
              const int qstar = j*nlsize*nccl_size + i*nccl_size + k;
              const int ZCCL1 = k; // cross redshift bin not supported so not using ZCCL1(k)
              const int ZCCL2 = k; // cross redshift bin not supported so not using ZCCL2(k)

              Cl[q][l] = (l > limits.P_2_s_min) ? C_cc_tomo_limber(l, i, j, ZCCL1, ZCCL2) :
                C_cc_tomo_limber_nointerp(l, i, j, ZCCL1, ZCCL2, use_linear_ps_limber);
              Cl[qstar][l] = Cl[q][l];
            }
          } 
        }
        else
        {
          const int L = 1;
          const double tol = 0.01; // required fractional accuracy in C(l)
          const double dev = 10. * tolerance; // will be diff  exact vs Limber init to
                                              // large value in order to start while loop
          for (int i=0; i<nlsize; i++) 
          { 
            for (int j=i; j<nlsize; j++)  
            {
              for (int k=0; k<nccl_size; k++)  
              {
                const int ZCCL1 = k; // cross redshift bin not supported so not using ZCCL1(k)
                const int ZCCL2 = k; // cross redshift bin not supported so not using ZCCL2(k)
                const int q = i*nlsize*nccl_size + j*nccl_size + k;
                C_cc_tomo(L, i, j, ZCCL1, ZCCL2, Cl[q], dev, tol);
              } 
              #pragma omp parallel for collapse(2)
              for (int k=0; k<nccl_size; k++)  
              {
                for (int l=limits.LMAX_NOLIMBER+1; l<nell; l++)
                {
                  const int ZCCL1 = k; // cross redshift bin not supported so not using ZCCL1(k)
                  const int ZCCL2 = k; // cross redshift bin not supported so not using ZCCL2(k)
                  const int q = i*nlsize*nccl_size + j*nccl_size + k;
                  Cl[q][l] = C_cc_tomo_limber(l, i, j, ZCCL1, ZCCL2)
                }
              }
              #pragma omp parallel for collapse(2)
              for (int k=0; k<nccl_size; k++)  
              {
                for (int l=1; l<nell; l++)
                {
                  const int q = i*nlsize*nccl_size + j*nccl_size + k;
                  const int qstar = j*nlsize*nccl_size + i*nccl_size + k;
                  Cl[qstar][l] = Cl[q][l];
                }
              }
            }
          } 
        }
      }
    }
    #pragma omp parallel for collapse(4)
    for (int i=0; i<nlsize; i++)
    { 
      for (int j=0; j<nlsize; j++)  
      {
        for (int k=0; k<nccl_size; k++)  
        {
          for (int p=0; p<ntheta; p++)
          {
            const int nz = i*nlsize*nccl_size + j*nccl_size + k;
            const int q = nz*ntheta + p;
            w_vec[q] = 0;
            for (int l=1; l<nell; l++)
            {
              w_vec[q] += Pl[p][l]*Cl[nz][l];
            }
          }
        }
      }
    }
    for (int i=0; i<NSIZE; i++)
    {
      free(Cl[i]);
    }
    free(Cl);

    update_cosmopara(&C);
    update_nuisance(&N);
  }
                                                        
  const int q = (nl1*nlsize*nccl_size + nl2*nccl_size + ni)*ntheta + nt; // cross redshift bin not 
                                                                         // supported so not using 
                                                                         // N_CCL(ni, nj) instead of ni
  if(q > NSIZE*ntheta - 1)
  {
    log_fatal("error in selecting bin number");
    exit(1);
  }
  return w_vec[q];
}

/*
double w_cg_tomo(int nt, int nl, int ni, int nj, int limber)
{
  if (like.Ntheta == 0)
  {
    log_fatal("like.Ntheta not initialized");
    exit(1);
  }
  if(nt > like.Ntheta - 1)
  {
    log_fatal("error in selecting bin number");
    exit(1); 
  } 

  static cosmopara C;
  static galpara G;
  static nuisancepara N;
  static double** Pl = 0;
  static double* w_vec = 0;

  const int nell = limits.LMAX;
  const int ntheta = like.Ntheta;

  if (Pl == 0)
  {
    Pl = (double**) malloc(sizeof(double*)*ntheta);
    for(int i=0; i<ntheta; i++) 
    {
      Pl[i] = (double*) malloc(sizeof(double)*nell);
    }
    w_vec = (double*) malloc(sizeof(double)*NSIZE*ntheta);
    double xmin[ntheta];
    double xmax[ntheta];
    // Cocoa: dont thread (init of static variables inside set_bin_average)
    for (int i=0; i<ntheta; i ++)
    {
      bin_avg r = set_bin_average(i,0);
      xmin[i] = r.xmin;
      xmax[i] = r.xmax;
    }
    double** Pmin = (double**) malloc(sizeof(double)*ntheta);
    double** Pmax = (double**) malloc(sizeof(double)*ntheta);
    for (int i=0; i<ntheta; i ++)
    {
      Pmin[i] = (double*) malloc(sizeof(double)*(nell + 1));
      Pmax[i] = (double*) malloc(sizeof(double)*(nell + 1));
    }
    #pragma omp parallel for
    for (int i=0; i<ntheta; i ++)
    {
      for (int l=0; l<nell; l++)
      {
        bin_avg r = set_bin_average(i,l);
        Pmin[i][l] = r.Pmin;
        Pmax[i][l] = r.Pmax;
      }
      Pl[i][0] = 1.0;
      const double tmp = (1.0/(xmin[i] - xmax[i]))*(1. / (4. * M_PI));
      for (int l=1; l<nell; l++)
      {
        Pl[i][l] = tmp*(Pmin[i][l + 1] - Pmax[i][l + 1] - Pmin[i][l - 1] + Pmax[i][l - 1]);
      }
    }
    for (int i=0; i<ntheta; i ++)
    {
      free(Pmin[i]);
      free(Pmax[i]);
    }
    free(Pmin);
    free(Pmax);
  }

  if (recompute_cg(C, G, N))
  {
    for (int i=0; i<NSIZE; i++)
    {
      free(Cl[i]);
    }
    free(Cl);

    update_cosmopara(&C);
    update_galpara(&G);
    update_nuisance(&N);
  } 

  if(q > NSIZE*ntheta - 1)
  {
    log_fatal("error in selecting bin number");
    exit(1);
  }
  return w_vec[q];
}
*/

// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// Correlation Functions (real space) - flat sky
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------

double w_gammat_cluster_tomo_flatsky(double theta, int nl, int ni, int nj, int limber) 
{ // nl = lambda_obs bin, ni = cluster redshift bin, nj = source redshift bin
  static cosmopara C;
  static nuisancepara N;
  static galpara G;
  static double** table;
 
  const int ntheta = Ntable.N_thetaH;
  const int nlsize = Cluster.N200_Nbin;
  const int ngammat_size = tomo.cgl_Npowerspectra;
  const int NSIZE = nlsize*ngammat_size;

  const double l_min = w_l_min;
  const double l_max = w_l_max;
  const double lnlmax = log(l_max);
  const double lnlmin = log(l_min);
  const double dlnl = (lnlmax-lnlmin)/(1.0*ntheta - 1.);
  const double lnrc = 0.5*(lnlmax + lnlmin);
  const double nc = ntheta/2.0 + 1;

  const double lnthetamin = (nc-ntheta+1)*dlnl-lnrc;
  const double lnthetamax = nc*dlnl-lnrc;
  const double dlntheta = (lnthetamax - lnthetamin)/((double) ntheta);
  const double lntheta = log(theta);

  if (table == 0)
  {
    table = (double**) malloc(sizeof(double*)*NSIZE);
    for(int i=0; i<NSIZE; i++)
    {
      table[i] = (double*) malloc(sizeof(double)*ntheta);
    }
  }

  if (recompute_cs(C, G, N))
  {
    typedef fftw_complex fftwZ;
    // go to log-Fourier-space
    fftwZ** flP = (fftwZ**) malloc(sizeof(fftwZ*)*NSIZE);
    for (int j=0; j<NSIZE; j++)
    {
      flP[j] = fftw_malloc((ntheta/2 + 1)*sizeof(fftwZ));
    }
    { 
      double** lP = (double**) malloc(sizeof(double*)*NSIZE);
      fftw_plan* plan = (fftw_plan*) malloc(sizeof(fftw_plan)*NSIZE);
      for (int j=0; j<NSIZE; j++)
      {
        lP[j] = (double*) malloc(ntheta*sizeof(double));
        plan[j] = fftw_plan_dft_r2c_1d(ntheta, lP[j], flP[j], FFTW_ESTIMATE);
      }

      for(int i=0; i<nlsize; i++) // Power spectrum on logarithmic bins
      { // loop in lambda_obs - (cannot thread loops in lambda_obs  - see cluster_utils)
        if (limber != 1)
        { 
          log_fatal("NonLimber not implemented");
          exit(1);
        }
        /*
        gsl_spline** fCL_NL = NULL; 
        if (limber != 1)
        { 
          fCL_NL = (gsl_spline**) malloc(sizeof(gsl_spline*)*ngammat_size);
          double** Cl_NL = (double**) malloc(sizeof(double*)*ngammat_size);
          for (int i=0; j<ngammat_size; j++)
          {
            Cl_NL[j] = calloc(limits.LMAX_NOLIMBER, sizeof(double));
          }
          double* ll = calloc(limits.LMAX_NOLIMBER, sizeof(double));
          for(int i=0; i<limits.LMAX_NOLIMBER; i++)
          {
            ll[i] = i;
          }

          const int L = 1;
          const double tolerance = 0.0075;    // required fractional accuracy in C(l)
          const double dev = 10. * tolerance; // will be diff exact vs Limber init to large
                                              // value in order to start while loop
          for (int k=0; k<ngammat_size; k++)
          { // Cocoa: no threading allowed here - (fftw allocation @C_gl_tomo)
            const int Z1 = ZL(k);
            const int Z2 = ZS(k);
            ????(L, Z1, Z2, Cl_NL[k], dev, tolerance); // no nonlimber

            const gsl_interp_type* T = gsl_interp_linear;
            fCL_NL[k] = gsl_spline_alloc(T, limits.LMAX_NOLIMBER);
            if (fCL_NL[k] == NULL)
            {
              log_fatal("fail allocation");
              exit(1);
            }
          }
          #pragma omp parallel for
          for (int k=0; k<ngammat_size; k++)
          {
            int status = gsl_spline_init(fCL_NL[k], ll, Cl_NL[k], limits.LMAX_NOLIMBER);
            if (status) 
            {
              log_fatal(gsl_strerror(status));
              exit(1);
            }
          }
          for (int i=0; j<ngammat_size; j++)
          {
            free(Cl_NL[j]);
          }
          free(Cl_NL);
          free(ll);
        }
        */
        
        C_cs_tomo_limber(limits.P_2_s_min + 1, i, ZCL(0), ZCS(0)); // need to init
        #pragma omp parallel for collapse (2);
        for(int j=0; j<ngammat_size; j++)
        { 
          for(int p=0; p<ntheta; p++)
          {
            const int ZC = ZCL(j);
            const int ZSC = ZCS(j);
            const int q = i*ngammat_size + j;
            const double l = exp(lnrc + (p - nc)*dlnl);
            if (limber == 1 || (limber != 1 && l > limits.LMAX_NOLIMBER - 1))
            {
              lP[q][p] = (l > limits.P_2_s_min) ? l*C_cs_tomo_limber(l, i, ZC, ZSC) :
                l*C_cs_tomo_limber_nointerp(l, i, ZC, ZSC, use_linear_ps_limber);
            }
            else
            {
              /*
              double CLNL;
              int status = gsl_spline_eval_e(fCL_NL[j], l, NULL, &CLNL);
              if (status) 
              {
                log_fatal(gsl_strerror(status));
                exit(1);
              }
              lP[q][p] = l*CLNL;
              */
            }
          }
        }
        /*
        if (limber != 1)
        {
          for (int i=0; i<ngammat_size; i++)
          {
             gsl_spline_free(fCL_NL[i]);
          }
          free(fCL_NL);
        }
        */   
      }

      #pragma omp parallel for
      for (int j=0; j<NSIZE; j++)
      { // Execute FFTW in parallel (thread-safe)
        fftw_execute(plan[j]); 
      }
      
      for (int j=0; j<NSIZE; j++)
      {
        fftw_free(lP[j]);
        fftw_destroy_plan(plan[j]);
      }
      free(lP);
      free(plan);
    }

    double** lP = (double**) malloc(sizeof(double*)*NSIZE);
    fftwZ** kernel = (fftwZ**) fftw_malloc(sizeof(fftwZ*)*NSIZE);
    fftwZ** conv = (fftwZ**) fftw_malloc(sizeof(fftwZ*)*NSIZE);
    fftw_plan* plan = (fftw_plan*) malloc(sizeof(fftw_plan)*NSIZE);
    double*** tab = (double***) malloc(sizeof(double**)*NSIZE);
    for (int j=0; j<NSIZE; j++)
    {
      lP[j] = (double*) malloc(ntheta*sizeof(double));
      kernel[j] = (fftwZ*) fftw_malloc(sizeof(fftwZ));
      conv[j] = (fftwZ*) fftw_malloc((ntheta/2+1)*sizeof(fftwZ));
      plan[j] = fftw_plan_dft_c2r_1d(ntheta, conv[j], lP[j], FFTW_ESTIMATE);
      tab[j] = (double**) malloc(sizeof(double*)*1);
      for(int i=0; i<1; i++) 
      {
        tab[j][i] = (double*) malloc(sizeof(double)*ntheta);
      }
    }

    #pragma omp parallel for
    for (int j=0; j<NSIZE; j++)
    { 
      double arg[2];
      arg[0] = 0; // bias
      arg[1] = 2; // order of Bessel function

      for(int i=0; i<(ntheta/2+1); i++)
      {
        const double kk = 2*M_PI*i/(dlnl*ntheta);
        hankel_kernel_FT(kk, kernel[j], arg, 2);
        conv[j][i][0] = flP[j][i][0]*kernel[j][0][0] - flP[j][i][1]*kernel[j][0][1];
        conv[j][i][1] = flP[j][i][1]*kernel[j][0][0] + flP[j][i][0]*kernel[j][0][1];
      }

      // force Nyquist- and 0-frequency-components to be double
      conv[j][0][1] = 0;
      conv[j][ntheta/2][1] = 0;

      fftw_execute(plan[j]); // Execute FFTW in parallel (thread-safe)

      for(int k=0; k<ntheta; k++)
      {
        const double t = exp((nc-k)*dlnl-lnrc); 
        tab[j][0][ntheta-k-1] = lP[j][k]/(t*2*M_PI*ntheta);
      }
      for (int k=0; k<ntheta; k++)
      {
        table[j][k] = tab[j][0][k];
      }
    }

    for (int j=0; j<NSIZE; j++)
    {
      fftw_free(flP[j]);
      fftw_free(lP[j]);
      fftw_free(conv[j]);
      fftw_free(kernel[j]);
      fftw_destroy_plan(plan[j]);
      free(tab[j]);
    }
    free(flP);
    free(lP);
    free(conv);
    free(kernel);
    free(plan);
    free(tab);

    update_galpara(&G);
    update_cosmopara(&C);
    update_nuisance(&N);
  } 
  
  if (!test_zoverlap_c(ni, nj)) 
  {
    return 0.0;
  }

  if (lntheta < lnthetamin || lntheta > lnthetamax)
  {
    const double theta = exp(lntheta);
    const double theta_min = exp(lnthetamin);
    const double theta_max = exp(lnthetamax);
    log_fatal("theta = %e outside look-up table range [%e, %e]", theta, theta_min, theta_max);
    exit(1);
  }

  const int q = nl*ngammat_size + N_cgl(ni, nj);
  if(q > NSIZE - 1)
  {
    log_fatal("error in selecting bin number");
    exit(1);
  }
  return interpol(table[q], ntheta, lnthetamin, lnthetamax, dlntheta, lntheta, 0, 0);
}

double w_cc_tomo_flatsky(double theta, int nl1, int nl2, int ni, int nj, int limber)
{ // nl{1,2} = lambda_obs bins, n{i,j} = cluster redshift bins
  static cosmopara C;
  static nuisancepara N;
  static double** table;

  const int nlsize = Cluster.N200_Nbin;
  const int nccl_size = tomo.cluster_Nbin; // cross redshift bin not supported so not using
                                           // tomo.cc_clustering_Npowerspectra
  const int NSIZE = nlsize*nlsize*nccl_size;
  const int ntheta = Ntable.N_thetaH;

  const double l_min = w_l_min;
  const double l_max = w_l_max;
  const double lnlmax = log(l_max);
  const double lnlmin = log(l_min);
  const double dlnl = (lnlmax-lnlmin)/(1.0*ntheta-1.);
  const double lnrc = 0.5*(lnlmax + lnlmin);
  const double nc = ntheta/2+1;

  const double lnthetamin = (nc-ntheta+1)*dlnl-lnrc;
  const double lnthetamax = nc*dlnl-lnrc;
  const double dlntheta = (lnthetamax - lnthetamin)/((double) ntheta);
  const double lntheta = log(theta);

  if (table == 0)
  {
    table = (double**) malloc(sizeof(double*)*NSIZE);
    for(int i=0; i<NSIZE; i++)
    {
      table[i] = (double*) malloc(sizeof(double)*ntheta);
    }
  } 
  if (recompute_cc(C, N))
  {
    typedef fftw_complex fftwZ;

    // go to log-Fourier-space
    fftwZ** flP = (fftwZ**) malloc(sizeof(fftwZ*)*NSIZE);
    for (int j=0; j<NSIZE; j++)
    {
      flP[j] = fftw_malloc((ntheta/2+1)*sizeof(fftwZ));
    }
    { 
      double** lP = (double**) malloc(sizeof(double*)*NSIZE);
      fftw_plan* plan = (fftw_plan*) malloc(sizeof(fftw_plan)*NSIZE);
      for (int j=0; j<NSIZE; j++)
      {
        lP[j] = (double*) malloc(ntheta*sizeof(double));
        plan[j] = fftw_plan_dft_r2c_1d(ntheta, lP[j], flP[j], FFTW_ESTIMATE);
      }
      for (int i=0; i<nlsize; i++) // Power spectrum on logarithmic bins begins
      { // loop in lambda_obs - (cannot thread loops in lambda_obs - see cluster_utils)
        for (int j=i; j<nlsize; j++)  
        { // loop in lambda_obs - (cannot thread loops in lambda_obs - see cluster_utils)
          gsl_spline** fCL_NL = NULL; 
          if (limber != 1)
          { 
            fCL_NL = (gsl_spline**) malloc(sizeof(gsl_spline*)*nccl_size);
            double** Cl_NL = (double**) malloc(sizeof(double*)*nccl_size);
            for (int i=0; j<nccl_size; j++)
            {
              Cl_NL[j] = calloc(limits.LMAX_NOLIMBER, sizeof(double));
            }
            double* ll = calloc(limits.LMAX_NOLIMBER, sizeof(double));
            for(int i=0; i<limits.LMAX_NOLIMBER; i++)
            {
              ll[i] = i;
            }

            const int L = 1;
            const double tol = 0.0075;    // required fractional accuracy in C(l)
            const double dev = 10. * tolerance; // will be diff exact vs Limber init to large
                                                // value in order to start while loop
            for (int k=0; k<nccl_size; k++)
            { // Cocoa: no threading allowed here - (fftw allocation @C_gl_tomo)
              const int ZCCL1 = k; // cross redshift bin not supported so not using ZCCL1(k)
              const int ZCCL2 = k; // cross redshift bin not supported so not using ZCCL1(k)
              C_cc_tomo(L, i, j, ZCCL1, ZCCL2, Cl[q], dev, tol);
              
              const gsl_interp_type* T = gsl_interp_linear;
              fCL_NL[k] = gsl_spline_alloc(T, limits.LMAX_NOLIMBER);
              if (fCL_NL[k] == NULL)
              {
                log_fatal("fail allocation");
                exit(1);
              }
            }
            #pragma omp parallel for
            for (int k=0; k<nccl_size; k++)
            {
              int status = gsl_spline_init(fCL_NL[k], ll, Cl_NL[k], limits.LMAX_NOLIMBER);
              if (status) 
              {
                log_fatal(gsl_strerror(status));
                exit(1);
              }
            }
            for (int i=0; j<nccl_size; j++)
            {
              free(Cl_NL[j]);
            }
            free(Cl_NL);
            free(ll);
          }

          C_cc_tomo_limber(limits.P_2_s_min+1, i, j, 0, 0); // need to init that (if limber)
          #pragma omp parallel for collapse(2)
          for (int k=0; k<nccl_size; k++)  
          {
            for(int p=0; p<ntheta; p++)
            {
              const int q = i*nlsize*nccl_size + j*nccl_size + k;
              const int qstar = j*nlsize*nccl_size + i*nccl_size + k;
              const int ZCCL1 = k; // cross redshift bin not supported so not using ZCCL1(k)
              const int ZCCL2 = k; // cross redshift bin not supported so not using ZCCL2(k)
              const double l = exp(lnrc + (p - nc)*dlnl);
              if (limber == 1 || (limber != 1 && l > limits.LMAX_NOLIMBER - 1))
              {
                lP[q][p] = (l > limits.P_2_s_min) ? l*C_cc_tomo_limber(l, i, j, ZCCL1, ZCCL2) :
                  l*C_cc_tomo_limber_nointerp(l, i, j, ZCCL1, ZCCL2, use_linear_ps_limber); 
              }
              else
              {
                double CLNL;
                int status = gsl_spline_eval_e(fCL_NL[k], l, NULL, &CLNL);
                if (status) 
                {
                  log_fatal(gsl_strerror(status));
                  exit(1);
                }
                lP[q][p] = l*CLNL;
              }
              lP[qstar][p] = lP[q][p];
            }
          }
        }
        if (limber != 1)
        {
          for (int i=0; i<nccl_size; i++)
          {
             gsl_spline_free(fCL_NL[i]);
          }
          free(fCL_NL);
        }
      } // Power spectrum on logarithmic bins ends

      #pragma omp parallel for
      for (int j=0; j<NSIZE; j++)
      {
        fftw_execute(plan[j]); // Execute FFTW in parallel (thread-safe)
      }

      for (int j=0; j<NSIZE; j++)
      {
        fftw_free(lP[j]);
        fftw_destroy_plan(plan[j]);
      }
      free(lP);
      free(plan);
    }
    
    double** lP = (double**) malloc(sizeof(double*)*NSIZE);
    fftwZ** kernel = (fftwZ**) fftw_malloc(sizeof(fftwZ*)*NSIZE);
    fftwZ** conv = (fftwZ**) fftw_malloc(sizeof(fftwZ*)*NSIZE);
    fftw_plan* plan = (fftw_plan*) malloc(sizeof(fftw_plan)*NSIZE);
    double*** tab = (double***) malloc(sizeof(double**)*NSIZE);
    for (int j=0; j<NSIZE; j++)
    {
      lP[j] = (double*) malloc(ntheta*sizeof(double));
      kernel[j] = (fftwZ*) fftw_malloc(sizeof(fftwZ));
      conv[j] = (fftwZ*) fftw_malloc((ntheta/2+1)*sizeof(fftwZ));
      plan[j] = fftw_plan_dft_c2r_1d(ntheta, conv[j], lP[j], FFTW_ESTIMATE);
      tab[j] = (double**) malloc(sizeof(double*)*1);
      for(int i=0; i<1; i++) 
      {
        tab[j][i] = (double*) malloc(sizeof(double)*ntheta);
      }
    }
    #pragma omp parallel for
    for (int j=0; j<NSIZE; j++)
    {
      double arg[2];
      arg[0] = 0; // bias
      arg[1] = 0; // order of Bessel function

      for(int i=0; i<(ntheta/2+1); i++)
      {
        const double kk = 2*M_PI*i/(dlnl*ntheta);
        hankel_kernel_FT(kk, kernel[j], arg, 2);
        conv[j][i][0] = flP[j][i][0]*kernel[j][0][0] - flP[j][i][1]*kernel[j][0][1];
        conv[j][i][1] = flP[j][i][1]*kernel[j][0][0] + flP[j][i][0]*kernel[j][0][1];
      }

      // force Nyquist- and 0-frequency-components to be double
      conv[j][0][1] = 0;
      conv[j][ntheta/2][1] = 0;

      fftw_execute(plan[j]); // Execute FFTW in parallel (thread-safe)

      for(int k=0; k<ntheta; k++)
      {
        const double t = exp((nc-k)*dlnl-lnrc); 
        tab[j][0][ntheta-k-1] = lP[j][k]/(t*2*M_PI*ntheta);
      }
      for (int k=0; k<ntheta; k++)
      {
        table[j][k] = tab[j][0][k];
      }
    }
    for (int j=0; j<NSIZE; j++)
    {
      fftw_free(flP[j]);
      fftw_free(lP[j]);
      fftw_free(conv[j]);
      fftw_free(kernel[j]);
      fftw_destroy_plan(plan[j]);
      free(tab[j]);
    }
    free(flP);
    free(lP);
    free(conv);
    free(kernel);
    free(plan);
    free(tab);

    update_cosmopara(&C);
    update_nuisance(&N);
  }

  if (ni != nj) 
  {
    log_fatal("cross-tomography not supported");
    exit(1);
  }

  if (lntheta < lnthetamin || lntheta > lnthetamax)
  {
    const double theta = exp(lntheta);
    const double theta_min = exp(lnthetamin);
    const double theta_max = exp(lnthetamax);
    log_fatal("theta = %e outside look-up table range [%e, %e]", theta, theta_min, theta_max);
    exit(1);
  }

  const int q = i*nlsize*nccl_size + j*nccl_size + k; // cross redshift bin not supported so not
                                                      // using N_CCL(ni, nj) instead of ni
  if(q > NSIZE - 1)
  {
    log_fatal("error in selecting bin number");
    exit(1);
  }
  return interpol(table[q], ntheta, lnthetamin, lnthetamax, dlntheta, lntheta, 0, 0);
}

/*
double w_cg_tomo_flatsky(double theta, int nl, int ni, int nj, int limber)
{ // nl = lambda_obs bin, ni = cluster bin, nj = galaxy bin
  static cosmopara C;
  static galpara G;
  static nuisancepara N;
  static double** table;

  const int nlsize = Cluster.N200_Nbin;
  const int nisize = tomo.cluster_Nbin;
  const int njsize = tomo.clustering_Nbin;
  const int NSIZE = nlsize*nisize*njsize;
  const int ntheta = Ntable.N_thetaH;

  const double l_min = w_l_min;
  const double l_max = w_l_max;
  const double lnlmax = log(l_max);
  const double lnlmin = log(l_min);
  const double dlnl = (lnlmax-lnlmin)/(1.0*ntheta-1.);
  const double lnrc = 0.5*(lnlmax + lnlmin);
  const double nc = ntheta/2+1;

  const double lnthetamin = (nc-ntheta+1)*dlnl-lnrc;
  const double lnthetamax = nc*dlnl-lnrc;
  const double dlntheta = (lnthetamax - lnthetamin)/((double) ntheta);
  const double lntheta = log(theta);

  if (table == 0)
  {
    table = (double**) malloc(sizeof(double*)*NSIZE);
    for(int i=0; i<NSIZE; i++)
    {
      table[i] = (double*) malloc(sizeof(double)*ntheta);
    }
  } 
  if (recompute_cg(C, G, N))
  {
    if (limber == 1)
    {
      typedef fftw_complex fftwZ;

      // go to log-Fourier-space
      fftwZ** flP = (fftwZ**) malloc(sizeof(fftwZ*)*NSIZE);
      for (int j=0; j<NSIZE; j++)
      {
        flP[j] = fftw_malloc((ntheta/2+1)*sizeof(fftwZ));
      }
      { 
        double** lP = (double**) malloc(sizeof(double*)*NSIZE);
        fftw_plan* plan = (fftw_plan*) malloc(sizeof(fftw_plan)*NSIZE);
        for (int j=0; j<NSIZE; j++)
        {
          lP[j] = (double*) malloc(ntheta*sizeof(double));
          plan[j] = fftw_plan_dft_r2c_1d(ntheta, lP[j], flP[j], FFTW_ESTIMATE);
        }

        // Power spectrum on logarithmic bins (begins)
        for (int i=0; i<nlsize; i++)  
        { // cannot thread loop on lambda_obs (see cluster_util)
          {
            const int j = 0;
            {
              const int k = 0;
              const int q = i*nlsize*nisize + j*nisize + k;
              {
                const int p = 0;
                const double l = exp(lnrc + (p - nc)*dlnl);
                lP[q][p] = l*C_cg_tomo_limber_nointerp(l, i, j, k, use_linear_ps_limber);
                C_cg_tomo_limber(l, i, j, k); // need to init that
              }
              #pragma omp parallel for
              for(int p=1; p<ntheta; p++)
              {
                const double l = exp(lnrc + (p - nc)*dlnl);
                lP[q][p] = (l > limits.P_2_s_min) ? l*C_cg_tomo_limber(l, i, j, k) :
                  l*C_cg_tomo_limber_nointerp(l, i, j, k, use_linear_ps_limber);
              }
            }
            #pragma omp parallel for
            for (int k=1; k<njsize; k++)  
            {
              for(int p=0; p<ntheta; p++)
              {
                const int q = i*nlsize + j*nlsize*nisize + k;
                const double l = exp(lnrc + (p - nc)*dlnl);
                lP[q][p] = (l > limits.P_2_s_min) ? l*C_cg_tomo_limber(l, i, j, k) :
                  l*C_cg_tomo_limber_nointerp(l, i, j, k, use_linear_ps_limber);
              }
            }
          }
          #pragma omp parallel for
          for (int j=1; j<nisize; j++)  
          {
            for (int k=0; k<njsize; k++)  
            {
              const int q = i*nlsize + nlsize*nisize*j + k;
              for(int p=0; p<ntheta; p++)
              {
                const double l = exp(lnrc + (p - nc)*dlnl);
                lP[q][p] = (l > limits.P_2_s_min) ? l*C_cg_tomo_limber(l, i, j, k) :
                  l*C_cg_tomo_limber_nointerp(l, i, j, k, use_linear_ps_limber);
              }
            }
          }
        }

        #pragma omp parallel for
        for (int j=0; j<NSIZE; j++)
        {
          fftw_execute(plan[j]); // Execute FFTW in parallel (thread-safe)
        }
        
        for (int j=0; j<NSIZE; j++)
        {
          fftw_free(lP[j]);
          fftw_destroy_plan(plan[j]);
        }
        free(lP);
        free(plan);
      }
      
      double** lP = (double**) malloc(sizeof(double*)*NSIZE);
      fftwZ** kernel = (fftwZ**) fftw_malloc(sizeof(fftwZ*)*NSIZE);
      fftwZ** conv = (fftwZ**) fftw_malloc(sizeof(fftwZ*)*NSIZE);
      fftw_plan* plan = (fftw_plan*) malloc(sizeof(fftw_plan)*NSIZE);
      double*** tab = (double***) malloc(sizeof(double**)*NSIZE);
      for (int j=0; j<NSIZE; j++)
      {
        lP[j] = (double*) malloc(ntheta*sizeof(double));
        kernel[j] = (fftwZ*) fftw_malloc(sizeof(fftwZ));
        conv[j] = (fftwZ*) fftw_malloc((ntheta/2+1)*sizeof(fftwZ));
        plan[j] = fftw_plan_dft_c2r_1d(ntheta, conv[j], lP[j], FFTW_ESTIMATE);
        tab[j] = (double**) malloc(sizeof(double*)*1);
        for(int i=0; i<1; i++) 
        {
          tab[j][i] = (double*) malloc(sizeof(double)*ntheta);
        }
      }

      #pragma omp parallel for
      for (int j=0; j<NSIZE; j++)
      {
        double arg[2];
        arg[0] = 0; // bias
        arg[1] = 0; // order of Bessel function

        for(int i=0; i<(ntheta/2+1); i++)
        {
          const double kk = 2*M_PI*i/(dlnl*ntheta);
          hankel_kernel_FT(kk, kernel[j], arg, 2);
          conv[j][i][0] = flP[j][i][0]*kernel[j][0][0] - flP[j][i][1]*kernel[j][0][1];
          conv[j][i][1] = flP[j][i][1]*kernel[j][0][0] + flP[j][i][0]*kernel[j][0][1];
        }

        // force Nyquist- and 0-frequency-components to be double
        conv[j][0][1] = 0;
        conv[j][ntheta/2][1] = 0;

        fftw_execute(plan[j]); // Execute FFTW in parallel (thread-safe)

        for(int k=0; k<ntheta; k++)
        {
          const double t = exp((nc-k)*dlnl-lnrc); 
          tab[j][0][ntheta-k-1] = lP[j][k]/(t*2*M_PI*ntheta);
        }
        for (int k=0; k<ntheta; k++)
        {
          table[j][k] = tab[j][0][k];
        }
      }
      
      for (int j=0; j<NSIZE; j++)
      {
        fftw_free(flP[j]);
        fftw_free(lP[j]);
        fftw_free(conv[j]);
        fftw_free(kernel[j]);
        fftw_destroy_plan(plan[j]);
        free(tab[j]);
      }
      free(flP);
      free(lP);
      free(conv);
      free(kernel);
      free(plan);
      free(tab);
    }
    else
    {
      log_fatal("NonLimber not implemented");
      exit(1);
    }

    update_galpara(&G);
    update_cosmopara(&C);
    update_nuisance(&N);
  }

  if (lntheta < lnthetamin || lntheta > lnthetamax)
  {
    const double theta = exp(lntheta);
    const double theta_min = exp(lnthetamin);
    const double theta_max = exp(lnthetamax);
    log_fatal("theta = %e outside look-up table range [%e, %e]", theta, theta_min, theta_max);
    exit(1);
  }

  const int q = nlsize*nl + nlsize*nisize*ni + nj;
  if(q > NSIZE-1)
  {
    log_fatal("error in selecting bin number");
    exit(1);
  }
  return interpol(table[q], ntheta, lnthetamin, lnthetamax, dlntheta, lntheta, 0, 0);
}
*/

// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// Limber Approximation (Angular Power Spectrum)
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// cluster lensing 
// ----------------------------------------------------------------------------
// nl = lambda_obs bin, ni = cluster redshift bin, nj = source redshift bin

double int_for_C_cs_tomo_limber_nointerp(double a, void* params)
{ 
  if(!(a>0) || !(a<1)) 
  {
    log_fatal("a>0 and a<1 not true");
    exit(1);
  }
  double* ar = (double *) params;

  const int nl = (int) ar[0];
  const int ni = (int) ar[1];
  const int nj = (int) ar[2];
  const double ell = ar[3] + 0.5;
  const int use_linear_ps = (int) ar[4];

  struct chis chidchi = chi_all(a);
  const double hoverh0 = hoverh0v2(a, chidchi.dchida);
  
  const double fK = f_K(chidchi.chi);
  const double k = ell/fK;
   
  const double Pcm = binned_p_cm(k, a, nl, use_linear_ps);
  const double tmp = W_cluster(ni, a, chidchi.chi, hoverh0)*W_kappa(a, fK, nj);
  return tmp*Pcm*chidchi.dchida/(fK*fK);
}


double C_cs_tomo_limber_nointerp(double l, int nl, int ni, int nj, int use_linear_ps)
{
  double ar[5] = {(double) nl, (double) ni, (double) nj, l, (double) use_linear_ps};
  
  const double zmin = tomo.cluster_zmin[ni];
  const double zmax = tomo.cluster_zmax[ni];
  const double amin = 1./(1. + zmax);
  const double amax = 1./(1. + zmin);
  
  if(zmin > zmax) 
  {
    return 0.;
  }
  else 
  {
    return int_gsl_integrate_low_precision(int_for_C_cs_tomo_limber_nointerp, ar, amin, amax, 
      NULL, GSL_WORKSPACE_SIZE);
  }
}

double C_cs_tomo_limber(double l, int nl, int ni, int nj)
{
  static cosmopara C;
  static galpara G;
  static nuisancepara N;
  static double** table = 0;

  const int nell = Ntable.N_ell;
  const int nlsize = Cluster.N200_Nbin;
  const int ngammat_size = tomo.cgl_Npowerspectra;
  const int NSIZE = nlsize*ngammat_size;

  const double lnlmin = log(limits.P_2_s_min);
  const double lnlmax = log(limits.P_2_s_max);
  const double dlnl = (lnlmax - lnlmin)/(nell-1);

  if (table == 0)
  {
    table = (double**) malloc(sizeof(double*)*NSIZE);
    for(int i=0; i<NSIZE; i++)
    {
      table[i] = (double*) malloc(sizeof(double)*nell);
    }
  }

  if (recompute_cs(C, G, N))
  {
    for(int i=0; i<nlsize; i++)
    { // loop in lambda_obs - (cannot thread loops in lambda_obs  - see cluster_utils)
      {
        const int j = 0;
        const int ZC = ZCL(j);
        const int ZS = ZCS(j);
        const int q = nlsize*i + j;
        {
          const int p=0;
          const double lnl = lnlmin + p*dlnl;
          const double l = exp(lnl);
          table[q][p] = log(C_cs_tomo_limber_nointerp(l, i, ZC, ZS, use_linear_ps_limber));
        }
        #pragma omp parallel for
        for (int p=1; p<nell; p++)
        {
          const double lnl = lnlmin + p*dlnl;
          const double l = exp(lnl);
          table[q][p] = log(C_cs_tomo_limber_nointerp(l, i, ZC, ZS, use_linear_ps_limber));
        }
      }
      #pragma omp parallel for
      for(int j=1; j<ngammat_size; j++)
      {
        const int ZC = ZCL(j);
        const int ZS = ZCS(j);
        const int q = nlsize*i + j;
        for (int p=0; p<nell; p++)
        {
          const double lnl = lnlmin + p*dlnl;
          const double l = exp(lnl);
          table[q][p] = log(C_cs_tomo_limber_nointerp(l, i, ZC, ZS, use_linear_ps_limber));
        }
      }
    }
    update_cosmopara(&C);
    update_nuisance(&N);
  }
  
  if (!test_zoverlap_c(ni, nj)) 
  {
    return 0.0;
  }

  const double lnl = log(l);
  if (lnl < lnlmin || lnl > lnlmax)
  {
    log_fatal("l = %e outside look-up table range [%e, %e]", l, exp(lnlmin), exp(lnlmax));
    exit(1);
  }

  const int q = nlsize*nl + N_cgl(ni, nj);
  if(q > NSIZE - 1)
  {
    log_fatal("error in selecting bin number");
    exit(1);
  }

  const double f1 = exp(interpol(table[q], nell, lnlmin, lnlmax, dlnl, lnl, 1, 1));
  if (isnan(f1)) 
  {
    return 0.0;
  }
  else
  {
    return f1;
  }
}

// ---------------------------------------------------------------------------------------------
// Cluster clustering 
// ---------------------------------------------------------------------------------------------
// nl{1,2} = lambda_obs bins, n{i,j} = cluster redshift bins

double int_for_C_cc_tomo_limber_nointerp(double a, void* params)
{
  if(!(a>0) || !(a<1)) 
  {
    log_fatal("a>0 and a<1 not true");
    exit(1);
  }
  double* ar = (double*) params;
   
  const int nl1 = (int) ar[0];
  const int nl2 = (int) ar[1];
  const int ni = (int) ar[2];
  const int nj = (int) ar[3];
  const double ell = ar[4] + 0.5;
  const int use_linear_ps = (int) ar[5];
  
  struct chis chidchi = chi_all(a);
  const double hoverh0 = hoverh0v2(a, chidchi.dchida);
  
  const double fK = f_K(chidchi.chi);
  const double k  = ell/fK;
  
  double res = W_cluster(ni, a, chidchi.chi, hoverh0)*W_cluster(nj, a, chidchi.chi, hoverh0);
  if (res != 0)
  {
    res = res*binned_p_cc(k, a, nl1, nl2, use_linear_ps);
  }

  return res*chidchi.dchida/(fK*fK);
}

double C_cc_tomo_limber_nointerp(double l, int nl1, int nl2, int ni, int nj, int use_linear_ps)
{ 
  double ar[6] = {(double) nl1, (double) nl2, (double) ni, (double) nj, l, (double) use_linear_ps};  
  const double zmin = fmax(tomo.cluster_zmin[ni], tomo.cluster_zmin[nj]);
  const double zmax = fmin(tomo.cluster_zmax[ni], tomo.cluster_zmax[nj]);
  const double amin = 1./(1. + zmax);
  const double amax = 1./(1. + zmin);

  if(zmin > zmax) 
  {
    return 0.;
  }
  else
  {
    return int_gsl_integrate_low_precision(int_for_C_cc_tomo_limber_nointerp, (void*) ar, 
      amin, amax, NULL, GSL_WORKSPACE_SIZE);
  }
}

double C_cc_tomo_limber(double l, int nl1, int nl2, int ni, int nj) 
{
  static cosmopara C;
  static nuisancepara N;
  static double** table = 0;

  if (ni != nj)
  {
    log_fatal("ni != nj tomography not supported");
    exit(1);
  }  

  const int nell = Ntable.N_ell;
  const int nlsize = Cluster.N200_Nbin;
  const int nccl_size = tomo.cc_clustering_Npowerspectra; 
  const int NSIZE = nccl_size*nlsize*nlsize;
  
  const double lnlmin = log(limits.P_2_s_min);
  const double lnlmax = log(limits.P_2_s_max);
  const double dl = (lnlmax - lnlmin)/(nell - 1); 

  if (table == 0)
  { 
    table = (double**) malloc(sizeof(double*)*NSIZE);
    for(int i=0; i<NSIZE; i++)
    {
      table[i] = (double*) malloc(sizeof(double)*nell);
    }   
  }

  if (recompute_cc(C, N))
  {
    for (int i=0; i<nlsize; i++) 
    { // loop in lambda_obs - (cannot thread loops in lambda_obs  - see cluster_utils)
      for (int j=0; j<i+1; j++) 
      { // loop in lambda_obs - (cannot thread loops in lambda_obs  - see cluster_utils)
        {
          const int k = 0;     
          const int ZCCL1 = k; // cross redshift bin not supported so not using ZCCL1(k)
          const int ZCCL2 = k; // cross redshift bin not supported so not using ZCCL2(k)
          const int q = nlsize*i + nlsize*nlsize*j + k;
          const int qstar = nlsize*j + nlsize*nlsize*i + k;
          {
            const int p = 0;
            const double lnl = lnlmin + p*dl;
            const double l = exp(lnl); 
            table[q][p] = log(C_cc_tomo_limber_nointerp(l, i, j, ZCCL1, ZCCL2, use_linear_ps_limber));
            table[qstar][p] = table[qstar][p];
          }
          #pragma omp parallel for
          for (int p=1; p<nell; ++p)
          {
            const double lnl = lnlmin + p*dl;
            const double l = exp(lnl); 
            table[q][p] = log(C_cc_tomo_limber_nointerp(l, i, j, ZCCL1, ZCCL2, use_linear_ps_limber));
            table[qstar][p] = table[q][p];
          }
        }
        #pragma omp parallel for
        for (int k=1; k<nccl_size; k++) // loop in ni (= nj)
        {
          const int ZCCL1 = k; // cross redshift bin not supported so not using ZCCL1(k)
          const int ZCCL2 = k; // cross redshift bin not supported so not using ZCCL2(k)
          const int q = nlsize*i + nlsize*nlsize*j + k;
          const int qstar = nlsize*j + nlsize*nlsize*i + k;
          for (int p=0; p<nell; ++p)
          {
            const double lnl = lnlmin + p*dl;
            const double l = exp(lnl); 
            table[q][p] = log(C_cc_tomo_limber_nointerp(l, i, j, ZCCL1, ZCCL2, use_linear_ps_limber));
            table[qstar][p] = table[q][p];
          }
        }
      }
    }

    update_cosmopara(&C);
    update_nuisance(&N);
  }

  const double lnl = log(l);
  if (lnl < lnlmin || lnl > lnlmax)
  {
    log_fatal("l = %e outside look-up table range [%e,%e]", l, exp(lnlmin), exp(lnlmax));
    exit(1);
  } 

  const int q = nlsize*nl1 + nlsize*nlsize*nl2 + ni; // cross redshift bin not supported so not
                                                     // using N_CCL(ni, nj) instead of ni
  if(q > NSIZE-1)
  {
    log_fatal("error in selecting bin number");
    exit(1);
  }

  const double f1 = exp(interpol(table[q], nell, lnlmin, lnlmax, dl, lnl, 1, 1));
  if (isnan(f1)) 
  {
    return 0.0;
  }
  return f1;
}


// ---------------------------------------------------------------------------------------------
// cluster x galaxy clustering
// ---------------------------------------------------------------------------------------------
// nl = lambda_obs bin, ni = cluster redshift bin, nj = galaxy redshift bin

double int_for_C_cg_tomo_limber_nointerp(double a, void* params)
{
  if(!(a>0) || !(a<1)) 
  {
    log_fatal("a>0 and a<1 not true");
    exit(1);
  }
  double* ar = (double*) params;

  const int nl = (int) ar[0];
  const int ni = (int) ar[1];
  const int nj = (int) ar[2];
  const double ell = ar[3] + 0.5;
  const int use_linear_ps = (int) ar[4];

  struct chis chidchi = chi_all(a);
  const double hoverh0 = hoverh0v2(a, chidchi.dchida);
  
  const  double fK  = f_K(chidchi.chi);
  const double k = ell/fK;

  const double Pcg = binned_p_cg(k, a, nl, nj, use_linear_ps);

  double Pcg_1l = 0.;
  if (gbias.b2[(int)nz2] || (clusterAnalysisChoice.nonlinear_bias>0))
  {
    // Cosmolike: if non-linear bias of galaxy sample is set, include 1-loop terms for P_cg
    // Cosmolike: using non-linear biases of cluster and galaxy sample
    const double g = growfac(a);
    const double g4 = g*g*g*g;
    
    const double double b1g = gbias.b1_function(z, nj);
    double b2g = gbias.b2[nj];
    if(b2g == 0) 
    {
      b2g = b2_from_b1(gbias.b[nj]);
    }
    double bs2g = gbias.bs2[nj];
    
    const double b1c = weighted_bias(nl, z);
    const double b2c = weighted_bias2(nl, z);
    const double bs2c = bs2_from_b1(b1c);
    
    Pcg_1l = 0.5*(b1c*b2g+b2c*b1g)*PT_d1d2(k) + 0.25*b2g*b2c*PT_d2d2(k) 
      + 0.5*(b1c*bs2g+b1g*bs2c)*PT_d1s2(k) +0.25*(b2c*bs2g+b2g*bs2c)*PT_d2s2(k) 
      + 0.25*(bs2g*bs2c)*PT_s2s2(k) +0.5*(b1c*b3nl_from_b1(b1g) 
      + b1g*b3nl_from_b1(b1c))*PT_d1d3(k);
    Pcg_1l *= g4;
  }

  const double tmp = W_cluster(ni, a, chidchi.chi, hoverh0)*W_HOD(a, nj, hoverh0);
  if (tmp != 0)
  {
    if (Pcg + Pcg_1l < 0) 
    { // Cosmolike: this makes the code much much faster like 1000 times faster
      return 0; 
    }
    return tmp*(Pcg + Pcg_1l)*chidchi.dchida/(fK*fK);
  }
  else
  {
    return 0.0; 
  }  
}

double C_cg_tomo_limber_nointerp(double l, int nl, int ni, int nj, int use_linear_ps) 
{
   double ar[5] = {(double) nl, (double) ni, (double) nj, l, (double) use_linear_ps};
   const double zmin = fmax(tomo.cluster_zmin[ni], tomo.clustering_zmin[nj]);
   const double zmax = fmin(tomo.cluster_zmax[ni], tomo.clustering_zmax[nj]);
   const double amin = 1./(1. + zmax);
   const double amax = 1./(1 + zmin);
   
   if (zmin > zmax) 
   {
      return 0.;
   }
   else
   {
      return int_gsl_integrate_low_precision(int_for_C_cg_tomo_limber_nointerp, (void*) ar, 
        amin, amax, NULL, GSL_WORKSPACE_SIZE);
   }
}

double C_cg_tomo_limber(double l, int nl, int ni, int nj)
{ // TODO: external ni, nj pairs
  static cosmopara C;
  static galpara G;
  static nuisancepara N;
  static double** table = 0;

  const int nell = Ntable.N_ell;
  const int nlsize = Cluster.N200_Nbin;
  const int nisize = tomo.cluster_Nbin;
  const int njsize = tomo.clustering_Nbin;
  const int NSIZE = nlsize*nisize*njsize;

  const double lnlmin = log(limits.P_2_s_min);
  const double lnlmax = log(limits.P_2_s_max);
  const double dl = (lnlmax - lnlmin)/(nell-1);

  if (table == 0)
  { 
    table = (double**) malloc(sizeof(double*)*NSIZE);
    for(int i=0; i<NSIZE; i++)
    {
      table[i] = (double*) malloc(sizeof(double)*nell);
    }      
  }
  if (recompute_cg(C, G, N))
  {
    for (int i=0; i<nlsize; i++) 
    { // loop in lambda_obs - (cannot thread loops in lambda_obs  - see cluster_utils)
      {
        const int j = 0;
        {
          const int k = 0;
          const int q = i*nlsize + nlsize*nisize*j + k;
          {
            const int p = 0;
            const double lnl = lnlmin + p*dl;
            const double l = exp(lnl);
            table[q][p] = log(C_cg_tomo_limber_nointerp(l, i, j, k, use_linear_ps_limber));
          }
          #pragma omp parallel for
          for (int p=1; p<nell; ++p)
          {
            const double lnl = lnlmin + p*dl;
            const double l = exp(lnl);
            table[q][p] = log(C_cg_tomo_limber_nointerp(l, i, j, k, use_linear_ps_limber));
          }
        }
        #pragma omp parallel for
        for (int k=1; k<njsize; k++)  
        {
          const int q = i*nlsize + nlsize*nisize*j + k;
          for (int p=0; p<nell; ++p)
          {
            const double lnl = lnlmin + p*dl;
            const double l = exp(lnl);
            table[q][p] = log(C_cg_tomo_limber_nointerp(l, i, j, k, use_linear_ps_limber));
          }
        }
      }
      #pragma omp parallel for
      for (int j=1; j<nisize; j++)  
      {
        for (int k=0; k<njsize; k++)  
        {
          const int q = i*nlsize + nlsize*nisize*j + k;
          for (int p=0; p<nell; ++p)
          {
            const double lnl = lnlmin + p*dl;
            const double l = exp(lnl);
            table[q][p] = log(C_cg_tomo_limber_nointerp(l, i, j, k, use_linear_ps_limber));
          }
        }
      }
    }

    update_galpara(&G);
    update_cosmopara(&C);
    update_nuisance(&N);
  }
   
  const double lnl = log(l);
  if (lnl < lnlmin || lnl > lnlmax)
  {
    log_fatal("l = %e outside look-up table range [%e, %e]", l, exp(lnlmin), exp(lnlmax));
    exit(1);
  } 

  const int q = nlsize*nl + nlsize*nisize*ni + nj;
  if(q > NSIZE-1)
  {
    log_fatal("error in selecting bin number");
    exit(1);
  }
  const double f1 = exp(interpol(table[q], nell, lnlmin, lnlmax, dl, lnl, 1, 1));
  if (isnan(f1)) 
  {
    return 0.0;
  }
  
  return f1;
}

// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// Non Limber (Angular Power Spectrum)
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// cluster x cluster clustering
// ----------------------------------------------------------------------------

void f_chi_for_Psi_cluster_cl(double* chi, int Nchi, double* fchi, int ni, int nl, double zmin,
double zmax)
{
  const double real_coverH0 = cosmology.coverH0 / cosmology.h0; // unit Mpc
  {
    const int i = 0;
    const double a = a_chi(chi[i]/real_coverH0 /* convert unit to c/H0 */);
    const double z = 1./a - 1.;
    const double tmp1 = zdistr_cluster(z, ni, nl);
    const double pf = (tmp1 < 0.) ? 0 : tmp1; ; // get rid of unphysical negatives
    fchi[i] = chi[i]*pf*growfac(a)*weighted_bias(nl, z)*hoverh0(a)/real_coverH0;
    if ((z < zmin) || (z > zmax))
    {
      fchi[i] = 0.;
    }
  }
  #pragma omp parallel for
  for (int i=1; i<Nchi; i++) 
  {
    const double a = a_chi(chi[i]/real_coverH0 /* convert unit to c/H0 */);
    const double z = 1./a - 1.;
    const double tmp1 = zdistr_cluster(z, ni, nl);
    const double pf = (tmp1 < 0.) ? 0 : tmp1; // get rid of unphysical negatives
    fchi[i] = chi[i]*pf*growfac(a)*weighted_bias(nl, z)*hoverh0(a)/real_coverH0;
    if ((z < zmin) || (z > zmax))
    {
      fchi[i] = 0.;
    }
  }
}

void f_chi_for_Psi_cluster_cl_RSD(double* chi, int Nchi, double* fchi, int ni, int nl, double zmin, 
double zmax)
{
  const double real_coverH0 = cosmology.coverH0 / cosmology.h0;
  {
    const int i = 0;
    const double a = a_chi(chi[i]/real_coverH0 /* convert unit to c/H0 */);
    const double z = 1./a - 1.;
    const double tmp1 = zdistr_cluster(z,ni, nlambda);
    const double pf = (tmp1 < 0.) ? 0 : tmp1; // get rid of unphysical negatives
    struct growths tmp2 = growfac_all(a);
    fchi[i] = -chi[i]*pf*tmp2.D*tmp2.f*hoverh0(a)/real_coverH0;
    if ((z < zmin) || (z > zmax))
    {
      fchi[i] = 0.;
    }
  }
  #pragma omp parallel for
  for(int i=1; i<Nchi; i++) 
  {
    const double a = a_chi(chi[i]/real_coverH0 /* convert unit to c/H0 */);
    const double z = 1./a - 1.;
    const double tmp1 = zdistr_cluster(z,ni, nlambda);
    const double pf = (tmp1 < 0.) ? 0 : tmp1; // get rid of unphysical negatives
    struct growths tmp2 = growfac_all(a);
    fchi[i] = -chi[i]*pf*tmp2.D*tmp2.f*hoverh0(a)/real_coverH0;
    if ((z < zmin) || (z > zmax))
    {
      fchi[i] = 0.;
    }
  }
}

void f_chi_for_Psi_cluster_cl_Mag(double* chi, int Nchi, double* fchi, int ni, int nl, double zmax) 
{
  const double real_coverH0 = cosmology.coverH0 / cosmology.h0;
  {
    const int i = 0;
    const double a = a_chi(chi[i]/real_coverH0 /* convert unit to c/H0 */);
    const double z = 1.0/a - 1.0;
    const double fK = f_K(chi[i]/real_coverH0);
    const double wmag = W_mag_cluster(a, fK, ni, nl);
    const double window_M = wmag/fK/(real_coverH0*real_coverH0);
    
    fchi[i] = window_M*growfac(a); // unit [Mpc^-2]    
    if(z > zmax)
    {
      fchi[i] = 0.;
    }
  }
  #pragma omp parallel for
  for(int i=1; i<Nchi; i++) 
  {
    const double a = a_chi(chi[i]/real_coverH0 /* convert unit to c/H0 */);
    const double z = 1.0/a - 1.0;
    const double fK = f_K(chi[i]/real_coverH0);
    const double wmag = W_mag_cluster(a, fK, ni, nl);
    const double window_M = wmag/fK/(real_coverH0*real_coverH0);
    
    fchi[i] = window_M*growfac(a); // unit [Mpc^-2]    
    if(z > zmax)
    {
      fchi[i] = 0.;
    }
  }
}

void C_cc_tomo(int L, int nl1, int nl2, int ni, int nj, double* Cl, double dev, double tol)
{ // nl{1,2} = lambda_obs bins, n{i,j} = cluster redshift bins
  if (ni != nj)
  {
    log_fatal("cross-spectrum not supported");
    exit(1);
  }

  static double** k1;
  static double** k2;
  static double** Fk1;
  static double** Fk2;
  static double** Fk1_Mag;
  static double** Fk2_Mag;
  static double* chi_ar;

  const int Nell_block = Ntable.NL_Nell_block;
  const int Nchi = Ntable.NL_Nchi;  
  int ell_ar[Nell_block];
  double f1_chi[Nchi];
  double f2_chi[Nchi];
  double f1_chi_RSD[Nchi];
  double f2_chi_RSD[Nchi];
  double f1_chi_Mag[Nchi];
  double f2_chi_Mag[Nchi];

  const double real_coverH0 = cosmology.coverH0 / cosmology.h0;
  const double chi_min = chi(1./(1.+0.002))*real_coverH0; // DIMENSIONELESS
  const double chi_max = chi(1./(1.+4.))*real_coverH0; // DIMENSIONELESS
  const double dlnchi = log(chi_max/chi_min) / (Nchi - 1.);
  const double dlnk = dlnchi;

  if(k1 == 0)
  { // COCOA: no need to create/destroy arrays with same size at every call
    k1 = (double**) malloc(Nell_block * sizeof(double*));
    k2 = (double**) malloc(Nell_block * sizeof(double*));
    Fk1 = (double**) malloc(Nell_block * sizeof(double*));
    Fk2 = (double**) malloc(Nell_block * sizeof(double*));
    Fk1_Mag = (double**) malloc(Nell_block * sizeof(double*));
    Fk2_Mag = (double**) malloc(Nell_block * sizeof(double*));
    for (int i=0; i<Nell_block; i++)
    {
      k1[i] = (double*) malloc(Nchi * sizeof(double));
      k2[i] = (double*) malloc(Nchi * sizeof(double));
      Fk1[i] = (double*) malloc(Nchi * sizeof(double));
      Fk2[i] = (double*) malloc(Nchi * sizeof(double));
      Fk1_Mag[i] = (double*) malloc(Nchi * sizeof(double));
      Fk2_Mag[i] = (double*) malloc(Nchi * sizeof(double));
    }
    chi_ar = (double*) malloc(Nchi * sizeof(double));
  }
  for (int i=0; i<Nchi; i++)
  {
    // chi_min and chi_max may be cosmology dependent
    chi_ar[i] = chi_min * exp(dlnchi * i);
  }
  #pragma omp parallel for
  for (int i=0; i<Nell_block; i++)
  {
    for (int j=0; j<Nchi; j++)
    {
      k1[i][j] = 0.0;
      k2[i][j] = 0.0;
      Fk1[i][j] = 0.0;
      Fk2[i][j] = 0.0;
      Fk1_Mag[i][j] = 0.0;
      Fk2_Mag[i][j] = 0.0;
    }
  }

  const double zmin = tomo.cluster_zmin[ni];
  const double zmax = tomo.cluster_zmax[ni];
  f_chi_for_Psi_cluster_cl(chi_ar, Nchi, f1_chi, ni, nl1, zmin, zmax);
  f_chi_for_Psi_cluster_cl_RSD(chi_ar, Nchi, f1_chi_RSD, ni, nl1, zmin, zmax);
  if (INCLUDE_MAG_IN_C_CC_NONLIMBER == 1)
  {
    f_chi_for_Psi_cluster_cl_Mag(chi_ar, Nchi, f1_chi_Mag, ni, nl1, zmax);
  }
  if (nl1 != nl2) 
  {
    f_chi_for_Psi_cluster_cl(chi_ar, Nchi, f2_chi, ni, nl2, zmin, zmax);
    f_chi_for_Psi_cluster_cl_RSD(chi_ar, Nchi, f1_chi_RSD, ni, nl2, zmin, zmax);
    if (INCLUDE_MAG_IN_C_CC_NONLIMBER == 1)
    {
      f_chi_for_Psi_cluster_cl_Mag(chi_ar, Nchi, f1_chi_Mag, ni, nl2, zmax);
    }
  }
  
  config cfg;
  my_config.nu = 1.;
  my_config.c_window_width = 0.25;
  my_config.derivative = 0;
  my_config.N_pad = 200;
  my_config.N_extrap_low = 0;
  my_config.N_extrap_high = 0;

  config cfg_RSD
  my_config_RSD.nu = 1.01;
  my_config_RSD.c_window_width = 0.25;
  my_config_RSD.derivative = 2;
  my_config_RSD.N_pad = 500;
  my_config_RSD.N_extrap_low = 0;
  my_config_RSD.N_extrap_high = 0;

  config cfg_Mag;
  my_config_Mag.nu = 1.;
  my_config_Mag.c_window_width = 0.25;
  my_config_Mag.derivative = 0;
  my_config_Mag.N_pad = 500;
  my_config_Mag.N_extrap_low = 0;
  my_config_Mag.N_extrap_high = 0;

  int i_block = 0;
    
  while ((fabs(dev) > tolerance) & (L < limits.LMAX_NOLIMBER))
  { 
    for(int i=0; i<Nell_block; i++) 
    {
      ell_ar[i]=i+i_block*Nell_block;
    }

    cfftlog_ells(chi_ar, f1_chi, Nchi, &cfg, ell_ar, Nell_block, k1, Fk1);
    cfftlog_ells_increment(chi_ar, f1_chi_RSD, Nchi, &cfg_RSD, ell_ar, Nell_block, k1, Fk1);
    if(nl1 != nl2) 
    {
      cfftlog_ells(chi_ar, f2_chi, Nchi, &cfg, ell_ar, Nell_block, k2, Fk2);
      cfftlog_ells_increment(chi_ar, f2_chi_RSD, Nchi, &cfg_RSD, ell_ar, Nell_block, k2, Fk2);
    }
    if(INCLUDE_MAG_IN_C_CC_NONLIMBER == 1)
    {
      cfftlog_ells(chi_ar, f1_chi_Mag, Nchi, &cfg_Mag, ell_ar, Nell_block, k1, Fk1_Mag);
      if (nl1 != nl2) 
      {
        cfftlog_ells(chi_ar, f2_chi_Mag, Nchi, &cfg_Mag, ell_ar, Nell_block, k2, Fk2_Mag);
      }
      #pragma omp parallel for collapse(2)
      for(int i=0; i<Nell_block; i++) 
      {
        for(int j=0; j<Nchi; j++) 
        { // magnification bias set to 1 for now.
          const double ell_prefactor = ell_ar[i]*(ell_ar[i]+1.);
          Fk1[i][j] += ell_prefactor*Fk1_Mag[i][j]/(k1[i][j]*k1[i][j]); 
          if (nl1 != nl2) 
          {
            Fk2[i][j] += ell_prefactor*Fk2_Mag[i][j]/(k2[i][j]*k2[i][j]) ;
          }
        }
      }
    }

    {
      const int i = 0;
      cl_temp = 0.;
      for(int j=0; j<Nchi; j++) 
      {
        const double  k1_cH0 = k1[i][j] * real_coverH0;
        const double PK = p_lin(k1_cH0, 1.0);
        cl_temp += (nl1 == nl2) ? Fk1[i][j]*Fk1[i][j]*k1_cH0*k1_cH0*k1_cH0*PK :
          Fk1[i][j]*Fk2[i][j]*k1_cH0*k1_cH0*k1_cH0*PK;
      }
      Cl[ell_ar[i]] = cl_temp * dlnk * 2./M_PI + 
        C_cc_tomo_limber_nointerp((double) ell_ar[i], nl1, nl2, ni, ni, use_linear_ps_limber) 
       -C_cc_tomo_limber_nointerp((double) ell_ar[i], nl1, nl2, ni, ni, 1);
    }
    #pragma omp parallel for
    for(int i=1; i<Nell_block; i++) 
    {
      cl_temp = 0.;
      for(int j=0; j<Nchi; j++) 
      {
        const double  k1_cH0 = k1[i][j] * real_coverH0;
        const double PK = p_lin(k1_cH0, 1.0);
        cl_temp += (nl1 == nl2) ? Fk1[i][j]*Fk1[i][j]*k1_cH0*k1_cH0*k1_cH0 *PK :
          Fk1[i][j]*Fk2[i][j]*k1_cH0*k1_cH0*k1_cH0 *PK;
      }
      Cl[ell_ar[i]] = cl_temp * dlnk * 2./M_PI + 
        C_cc_tomo_limber_nointerp((double) ell_ar[i], nl1, nl2, ni, ni, use_linear_ps_limber) 
       -C_cc_tomo_limber_nointerp((double) ell_ar[i], nl1, nl2, ni, ni, 1);
    }

    i_block++;
    
    if(L >= limits.LMAX_NOLIMBER - Nell_block)
    { //Xiao: break before memory leak in next iteration
      break;
    }
    L = i_block*Nell_block - 1 ;
    dev = Cl[L]/C_cc_tomo_limber_nointerp((double) L, nl1, nl2, ni, ni, use_linear_ps_limber) - 1;
  }   
  L++;
  
  Cl[limits.LMAX_NOLIMBER+1] = C_cc_tomo_limber((double) limits.LMAX_NOLIMBER+1, nl1, nl2, ni, ni);
  #pragma omp parallel for
  for (int l=L; l<limits.LMAX_NOLIMBER; l++)
  {
    Cl[l] = (l > limits.P_2_s_min) ? C_cc_tomo_limber((double) l, nl1, nl2, ni, ni) :
      C_cc_tomo_limber_nointerp((double) l, nl1, nl2, ni, ni, use_linear_ps_limber);
  }
}

// ----------------------------------------------------------------------------
// cluster x galaxy clustering
// ----------------------------------------------------------------------------

void C_cg_tomo(int L, int nl, int ni, int nj, double* Cl, double dev, double tol)
{ // nl = lambda_obs bin, ni = cluster redshift bin, nj = galaxy redshift bin
  if (ni != nj)
  {
    log_fatal("cross-spectrum not supported");
    exit(1);
  }

  static double** k1;
  static double** k2;
  static double** Fk1;
  static double** Fk2;
  static double** Fk1_Mag;
  static double** Fk2_Mag;
  static double* chi_ar;

  const int Nell_block = Ntable.NL_Nell_block;
  const int Nchi = Ntable.NL_Nchi;  
  int ell_ar[Nell_block];
  double f1_chi[Nchi];
  double f2_chi[Nchi];
  double f1_chi_RSD[Nchi];
  double f2_chi_RSD[Nchi];
  double f1_chi_Mag[Nchi];
  double f2_chi_Mag[Nchi];

  const double real_coverH0 = cosmology.coverH0 / cosmology.h0;
  const double chi_min = chi(1./(1.+0.002))*real_coverH0; // DIMENSIONELESS
  const double chi_max = chi(1./(1.+4.))*real_coverH0; // DIMENSIONELESS
  const double dlnchi = log(chi_max/chi_min) / (Nchi - 1.);
  const double dlnk = dlnchi;

  if(k1 == 0)
  { // COCOA: no need to create/destroy arrays with same size at every call
    k1 = (double**) malloc(Nell_block * sizeof(double*));
    k2 = (double**) malloc(Nell_block * sizeof(double*));
    Fk1 = (double**) malloc(Nell_block * sizeof(double*));
    Fk2 = (double**) malloc(Nell_block * sizeof(double*));
    Fk1_Mag = (double**) malloc(Nell_block * sizeof(double*));
    Fk2_Mag = (double**) malloc(Nell_block * sizeof(double*));
    for (int i=0; i<Nell_block; i++)
    {
      k1[i] = (double*) malloc(Nchi * sizeof(double));
      k2[i] = (double*) malloc(Nchi * sizeof(double));
      Fk1[i] = (double*) malloc(Nchi * sizeof(double));
      Fk2[i] = (double*) malloc(Nchi * sizeof(double));
      Fk1_Mag[i] = (double*) malloc(Nchi * sizeof(double));
      Fk2_Mag[i] = (double*) malloc(Nchi * sizeof(double));
    }
    chi_ar = (double*) malloc(Nchi * sizeof(double));
  }
  for (int i=0; i<Nchi; i++)
  {
    // chi_min and chi_max may be cosmology dependent
    chi_ar[i] = chi_min * exp(dlnchi * i);
  }
  #pragma omp parallel for
  for (int i=0; i<Nell_block; i++)
  {
    for (int j=0; j<Nchi; j++)
    {
      k1[i][j] = 0.0;
      k2[i][j] = 0.0;
      Fk1[i][j] = 0.0;
      Fk2[i][j] = 0.0;
      Fk1_Mag[i][j] = 0.0;
      Fk2_Mag[i][j] = 0.0;
    }
  }

  const double zmin = fmax(tomo.cluster_zmin[ni], tomo.clustering_zmin[nj]);
  const double zmax = fmin(tomo.cluster_zmax[ni], tomo.clustering_zmax[nj]);

  f_chi_for_Psi_cluster_cl(chi_ar, Nchi, f1_chi, ni, Nlambda1, zmin, zmax);
  f_chi_for_Psi_cluster_cl_RSD(chi_ar, Nchi, f1_chi_RSD_ar, ni, Nlambda1, zmin, zmax);
  if(INCLUDE_MAG_IN_C_CC_NONLIMBER)
  {
    f_chi_for_Psi_cluster_cl_Mag(chi_ar, Nchi, f1_chi_Mag, ni, Nlambda1, zmin, zmax);
  }
  
  f_chi_for_Psi_cl(chi_ar, Nchi, f2_chi, nj, zmin, zmax);
  f_chi_for_Psi_cl_RSD(chi_ar, Nchi, f2_chi_RSD_ar, nj, zmin, zmax);
  if(INCLUDE_MAG_IN_C_CC_NONLIMBER)
  {
    f_chi_for_Psi_cl_Mag(chi_ar, Nchi, f2_chi_Mag, nj, zmin, zmax);
  }

  config cfg;
  my_config.nu = 1.;
  my_config.c_window_width = 0.25;
  my_config.derivative = 0;
  my_config.N_pad = 200;
  my_config.N_extrap_low = 0;
  my_config.N_extrap_high = 0;

  config cfg_RSD
  my_config_RSD.nu = 1.01;
  my_config_RSD.c_window_width = 0.25;
  my_config_RSD.derivative = 2;
  my_config_RSD.N_pad = 500;
  my_config_RSD.N_extrap_low = 0;
  my_config_RSD.N_extrap_high = 0;

  config cfg_Mag;
  my_config_Mag.nu = 1.;
  my_config_Mag.c_window_width = 0.25;
  my_config_Mag.derivative = 0;
  my_config_Mag.N_pad = 500;
  my_config_Mag.N_extrap_low = 0;
  my_config_Mag.N_extrap_high = 0;

  int i_block = 0;
    
  while ((fabs(dev) > tolerance) & (L < limits.LMAX_NOLIMBER))
  { 
    for(int i=0; i<Nell_block; i++) 
    {
      ell_ar[i]=i+i_block*Nell_block;
    } 

    cfftlog_ells(chi_ar, f1_chi, Nchi, &cfg, ell_ar, Nell_block, k1, Fk1);
    cfftlog_ells_increment(chi_ar, f1_chi_RSD_ar, Nchi, &cfg_RSD, ell_ar, Nell_block, k1, Fk1);
    cfftlog_ells(chi_ar, f2_chi, Nchi, &cfg, ell_ar, Nell_block, k2, Fk2);
    cfftlog_ells_increment(chi_ar, f2_chi_RSD_ar, Nchi, &cfg_RSD, ell_ar, Nell_block, k2, Fk2);
    if(INCLUDE_MAG_IN_C_CC_NONLIMBER)
    {
      cfftlog_ells(chi_ar, f1_chi_Mag, Nchi, &cfg_Mag, ell_ar, Nell_block, k1, Fk1_Mag);
      cfftlog_ells(chi_ar, f2_chi_Mag, Nchi, &cfg_Mag, ell_ar, Nell_block, k2, Fk2_Mag);
      #pragma omp parallel for collapse(2)
      for(int i=0;i<Nell_block;i++) 
      {
        for(int j=0; j<Nchi; j++) 
        {  // magnification bias set to 1 for now (cluster).
          const double ell_prefactor = ell_ar[i]*(ell_ar[i]+1.);
          Fk1[i][j] += (ell_prefactor/(k1[i][j]*k1[i][j])*Fk1_Mag[i][j]);
          Fk2[i][j] += gbias.b_mag[nj]*(ell_prefactor/(k2[i][j]*k2[i][j])*Fk2_Mag[i][j]);
        }
      }
    }

    {
      const int i = 0;
      double cl_temp = 0.;
      for(int j=0; j<Nchi; j++) 
      {
        const double k1_cH0 = k1[i][j] * real_coverH0;
        cl_temp += Fk1[i][j]*Fk2[i][j]*k1_cH0*k1_cH0*k1_cH0 *p_lin(k1_cH0, 1.0);
      }
      Cl[ell_ar[i]] = cl_temp * dlnk * 2./M_PI + 
        C_cg_tomo_limber_nointerp((double) ell_ar[i], nl, ni, nj, use_linear_ps_limber) 
       -C_cg_tomo_limber_nointerp((double) ell_ar[i], nl, ni, nj, 1);
    }
    #pragma omp parallel for
    for(int i=1; i<Nell_block; i++) 
    {
      double cl_temp = 0.;
      for(int j=0; j<Nchi; j++) 
      {
        const double k1_cH0 = k1[i][j] * real_coverH0;
        cl_temp += Fk1[i][j]*Fk2[i][j]*k1_cH0*k1_cH0*k1_cH0 *p_lin(k1_cH0, 1.0);
      }
      Cl[ell_ar[i]] = cl_temp * dlnk * 2./M_PI + 
        C_cg_tomo_limber_nointerp((double) ell_ar[i], nl, ni, nj, use_linear_ps_limber) 
       -C_cg_tomo_limber_nointerp((double) ell_ar[i], nl, ni, nj, 1);
    }

    i_block++;  
    if(L >= limits.LMAX_NOLIMBER - Nell_block)
    { //Xiao: break before memory leak in next iteration
      break;
    }
    L = i_block*Nell_block - 1 ;
    dev = Cl[L]/C_cg_tomo_limber_nointerp((double) L, nl, ni, nj, use_linear_ps_limber) - 1;
  }   
  L++;

  Cl[limits.LMAX_NOLIMBER + 1] = C_cg_tomo_limber((double) limits.LMAX_NOLIMBER + 1, nl, ni, nj);
  #pragma omp parallel for
  for (int l=L; l<limits.LMAX_NOLIMBER; l++)
  {
    Cl[l] = (l > limits.P_2_s_min) ? C_cg_tomo_limber((double) l, nl, ni, nj) :
      C_cg_tomo_limber_nointerp((double) l, nl, ni, nj, use_linear_ps_limber);
  }
}

// ----------------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------------
// cluster number counts
// ----------------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------------
// nl = lambda_obs bin, ni = cluster redshift bin

double int_projected_average_number_counts(double a, void* params)
{
  double* ar = (double*) params;   

  if(!(a>0)) 
  {
    log_fatal("a > 0 not true");
    exit(1);
  }
  const double z = 1./a-1 ;
  const int nl = (int) ar[0];
  const int interpolate_survey_area = (int) ar[1];
  const double norm = interpolate_survey_area > 0 ? get_area(z) : survey.area;  
  struct chis chidchi = chi_all(a);  
  const double fK = f_K(chidchi.chi);
  
  return fK*fK*chidchi.dchida*binned_average_number_counts(nl, z)*norm;
}

double projected_average_number_counts(int nl, int ni)
{ // nl = lambda_obs bin, ni = cluster redshift bin
  const int interpolate_survey_area = survey.area < 0 ? 1 : 0;
  double params[2] = {(double) nl, interpolate_survey_area};

  const double amin = 1./(1 + tomo.cluster_zmax[ni]);
  const double amax = 1./(1 + tomo.cluster_zmin[ni]);
  const double tmp = (4.0*M_PI/41253.0);

  return tmp*int_gsl_integrate_low_precision(int_projected_average_number_counts, 
    (void*) params, amin, amax, NULL, GSL_WORKSPACE_SIZE);
}

