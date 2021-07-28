#include <assert.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../cfftlog/cfftlog.h"
#include <fftw3.h>

#include "bias.h"
#include "basics.h"
#include "cosmo3D.h"
#include "cosmo2D.h"
#include "IA.h"
#include "pt_cfastpt.h"
#include "radial_weights.h"
#include "recompute.h"
#include "redshift_spline.h"
#include "structs.h"

#include "log.c/src/log.h"

static int GSL_WORKSPACE_SIZE = 250;
static int use_linear_ps_limber = 0; /* 0 or 1 */
static int include_RSD_GS = 0; /* 0 or 1 */
static int include_RSD_GG = 1; /* 0 or 1 */
static int include_RSD_GK = 0; /* 0 or 1 */
static double w_l_min = 0.0001;
static double w_l_max = 5.0e6;

double beam_cmb(const double l)
{
  const double sigma = cmb.fwhm/sqrt(8.0*M_LN2);
  return exp(-0.5*l*l*sigma*sigma);
}

double C_gk_tomo_limber_nointerp_wrapper(double l, int ni, int use_linear_ps)
{
  return C_gk_tomo_limber_nointerp(l, ni, use_linear_ps)*beam_cmb(l);
}

double C_gk_tomo_limber_wrapper(double l, int ni)
{
  return C_gk_tomo_limber(l, ni)*beam_cmb(l);
}

double C_ks_tomo_limber_nointerp_wrapper(double l, int ni, int use_linear_ps)
{
  return C_ks_tomo_limber_nointerp(l, ni, use_linear_ps)*beam_cmb(l);
}

double C_ks_tomo_limber_wrapper(double l, int ni)
{
  return C_ks_tomo_limber(l, ni)*beam_cmb(l);
}

// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// Correlation Functions (real Space) - Full Sky - bin average
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------

double xi_pm_tomo(int pm, int nt, int ni, int nj, int limber)
{
  if(ni < -1 || ni > tomo.shear_Nbin -1 || nj < -1 || nj > tomo.shear_Nbin -1)
  {
    log_fatal("invalid bin input (ni, nj) = (%d, %d)", ni, nj);
    exit(1);
  }
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

  static double** Glplus = 0;
  static double** Glminus = 0;
  static double* xi_vec_plus = 0;
  static double* xi_vec_minus = 0;
  static cosmopara C;
  static nuisancepara N;
  
  const int NSIZE = tomo.shear_Npowerspectra;
  const int nell = limits.LMAX;
  const int ntheta = like.Ntheta;

  if (Glplus == 0)
  {
    Glplus = (double**) malloc(sizeof(double*)*ntheta);
    for(int i=0; i<ntheta; i++) 
    {
      Glplus[i] = (double*) malloc(sizeof(double)*nell);
    }
    Glminus = (double**) malloc(sizeof(double*)*ntheta);
    for(int i=0; i<ntheta; i++) 
    {
      Glminus[i] = (double*) malloc(sizeof(double)*nell);
    }
    xi_vec_plus = (double*) malloc(sizeof(double)*NSIZE*ntheta);
    xi_vec_minus = (double*) malloc(sizeof(double)*NSIZE*ntheta);

    double xmin[ntheta];
    double xmax[ntheta];

    for (int i=0; i<ntheta; i++)
    { // Cocoa: dont thread (init of static variables inside set_bin_average)
      bin_avg r = set_bin_average(i,0);
      xmin[i] = r.xmin;
      xmax[i] = r.xmax;
    }

    double** Pmin = (double**) malloc(sizeof(double)*ntheta);
    double** Pmax = (double**) malloc(sizeof(double)*ntheta);
    double** dPmin = (double**) malloc(sizeof(double)*ntheta);
    double** dPmax = (double**) malloc(sizeof(double)*ntheta);
    for (int i=0; i<ntheta; i ++)
    {
      Pmin[i] = (double*) malloc(sizeof(double)*(nell + 1));
      Pmax[i] = (double*) malloc(sizeof(double)*(nell + 1));
      dPmin[i] = (double*) malloc(sizeof(double)*(nell + 1));
      dPmax[i] = (double*) malloc(sizeof(double)*(nell + 1));
    }
    #pragma omp parallel for collapse(2)
    for (int i=0; i<ntheta; i ++)
    {
      for (int l=0; l<nell; l++)
      {
        bin_avg r = set_bin_average(i, l);
        Pmin[i][l] = r.Pmin;
        Pmax[i][l] = r.Pmax;
        dPmin[i][l] = r.dPmin;
        dPmax[i][l] = r.dPmax;
      }
    }
    #pragma omp parallel for collapse(2)
    for (int i=0; i<ntheta; i ++)
    {
      for (int l=3; l<nell; l++)
      {
        Glplus[i][l] = (2.*l+1)/(2.*M_PI*l*l*(l+1)*(l+1))*(
          -l*(l-1.)/2*(l+2./(2*l+1)) * (Pmin[i][l-1]-Pmax[i][l-1])
          -l*(l-1.)*(2.-l)/2 * (xmin[i]*Pmin[i][l]-xmax[i]*Pmax[i][l])
          +l*(l-1.)/(2.*l+1) * (Pmin[i][l+1]-Pmax[i][l+1])
          +(4-l)*(dPmin[i][l]-dPmax[i][l])
          +(l+2)*(xmin[i]*dPmin[i][l-1] - xmax[i]*dPmax[i][l-1] - Pmin[i][l-1] + Pmax[i][l-1])
          +2*(l-1)*(xmin[i]*dPmin[i][l] - xmax[i]*dPmax[i][l] - Pmin[i][l] + Pmax[i][l])
          -2*(l+2)*(dPmin[i][l-1]-dPmax[i][l-1])
        )/(xmin[i]-xmax[i]);

        Glminus[i][l] = (2.*l+1)/(2.*M_PI*l*l*(l+1)*(l+1))*(
          -l*(l-1.)/2*(l+2./(2*l+1)) * (Pmin[i][l-1]-Pmax[i][l-1])
          -l*(l-1.)*(2.-l)/2 * (xmin[i]*Pmin[i][l]-xmax[i]*Pmax[i][l])
          +l*(l-1.)/(2.*l+1)* (Pmin[i][l+1]-Pmax[i][l+1])
          +(4-l)*(dPmin[i][l]-dPmax[i][l])
          +(l+2)*(xmin[i]*dPmin[i][l-1] - xmax[i]*dPmax[i][l-1] - Pmin[i][l-1] + Pmax[i][l-1])
          -2*(l-1)*(xmin[i]*dPmin[i][l] - xmax[i]*dPmax[i][l] - Pmin[i][l] + Pmax[i][l])
          +2*(l+2)*(dPmin[i][l-1]-dPmax[i][l-1])
          )/(xmin[i]-xmax[i]);
      }
    }
    for (int i=0; i<ntheta; i++)
    {
      Pmin[i];
      Pmax[i];
      dPmin[i];
      dPmax[i];
    }
    free(Pmin);
    free(Pmax);
    free(dPmin);
    free(dPmax);
  }

  if (recompute_shear(C, N))
  {
    if(limber == 1)
    {
      if (like.IA == 5 || like.IA == 6)
      { // NEW TATT MODELING
        double** Cl_EE = (double**) malloc(NSIZE*sizeof(double*));
        double** Cl_BB = (double**) malloc(NSIZE*sizeof(double*));
        for (int nz = 0; nz<NSIZE; nz++)
        {
          Cl_EE[nz] = calloc(nell, sizeof(double));
          Cl_BB[nz] = calloc(nell, sizeof(double));
        }
        
        { // init the functions C_ss_tomo_TATT_EE/BB_limber
          // only compute BB if the TATT parameters allow for B-mode terms
          const int BM = (nuisance.b_ta_z[0] || nuisance.b_ta_z[Z1NZ] ||
                          nuisance.b_ta_z[Z2NZ] || nuisance.A2_ia ||
                          nuisance.A2_z[Z1NZ] || nuisance.A2_z[Z2NZ]) ? 1 : 0;
          double x = C_ss_tomo_TATT_EE_limber(limits.LMIN_tab + 1, Z1(0), Z2(0)); 
          x = (BM == 1) ? C_ss_tomo_TATT_BB_limber(limits.LMIN_tab + 1, Z1(0), Z2(0)) : 0.0;
        }
        #pragma omp parallel for collapse(2)
        for (int nz=0; nz<NSIZE; nz++) 
        {
          for (int l=2; l<nell; l++)
          {
            const int Z1NZ = Z1(nz);
            const int Z2NZ = Z2(nz);
          
            // only compute BB if the TATT parameters allow for B-mode terms
            const int BM = (nuisance.b_ta_z[0] || nuisance.b_ta_z[Z1NZ] ||
                          nuisance.b_ta_z[Z2NZ] || nuisance.A2_ia ||
                          nuisance.A2_z[Z1NZ] || nuisance.A2_z[Z2NZ]) ? 1 : 0;

            Cl_EE[nz][l] = (l < limits.LMIN_tab + 1) ? 
              C_ss_tomo_TATT_EE_limber_nointerp((double) l, Z1NZ, Z2NZ) :
              C_ss_tomo_TATT_EE_limber((double) l,Z1NZ, Z2NZ);

            Cl_BB[nz][l] = (BM == 1) ? (l < limits.LMIN_tab + 1) ? 
              C_ss_tomo_TATT_BB_limber_nointerp((double) l, Z1NZ, Z2NZ) : 
              C_ss_tomo_TATT_BB_limber((double) l, Z1NZ, Z2NZ) : 0.0;
          }
        }
        #pragma omp parallel for collapse(2)
        for (int nz=0; nz<NSIZE; nz++)
        {
          for (int i=0; i<ntheta; i++)
          {
            const int q = nz*ntheta + i;
            xi_vec_plus[q] = 0;
            xi_vec_minus[q] = 0;
            for (int l=2; l<nell; l++)
            {
              xi_vec_plus[q] += Glplus[i][l] * (Cl_EE[nz][l] + Cl_BB[nz][l]);
              xi_vec_minus[q] += Glminus[i][l] * (Cl_EE[nz][l] - Cl_BB[nz][l]);
            }
          }
        }
        for (int nz=0; nz<NSIZE; nz++)
        {
          free(Cl_EE[nz]);
          free(Cl_BB[nz]);
        }
        free(Cl_EE);
        free(Cl_BB);
      }
      else
      {
        double** Cl = malloc(NSIZE*sizeof(double*));
        for (int nz=0; nz<NSIZE; nz++)
        {
          Cl[nz] = calloc(nell, sizeof(double));
        }

        C_ss_tomo_limber(limits.LMIN_tab, Z1(0), Z1(0)); // init the function
        #pragma omp parallel for collapse(2)
        for (int nz=0; nz<NSIZE; nz++)
        {
          for (int l=2; l<nell; l++)
          {
            const int Z1NZ = Z1(nz);
            const int Z2NZ = Z2(nz);
            Cl[nz][l] = (l < limits.LMIN_tab + 1) ? 
              C_ss_tomo_limber_nointerp((double) l, Z1NZ, Z2NZ, use_linear_ps_limber) :
              C_ss_tomo_limber((double) l, Z1NZ, Z2NZ);
          }
        }
        #pragma omp parallel for collapse(2)
        for (int nz=0; nz<NSIZE; nz++)
        {
          for (int i=0; i<ntheta; i++)
          {
            const int q = nz*ntheta + i;
            xi_vec_plus[q] = 0;
            xi_vec_minus[q] = 0;
            for (int l=2; l<nell; l++)
            {
              xi_vec_plus[q]  += Glplus[i][l]*Cl[nz][l];
              xi_vec_minus[q] += Glminus[i][l]*Cl[nz][l];
            }
          }
        }
        for (int nz=0; nz<NSIZE; nz++)
        {
          free(Cl[nz]);
        }
        free(Cl);
      }
    }
    else
    {
      log_fatal("NonLimber not implemented");
      exit(1);
    }
    update_cosmopara(&C);
    update_nuisance(&N);
  }

  const int q = N_shear(ni, nj)*ntheta + nt;
  if(q > NSIZE*ntheta - 1)
  {
    log_fatal("error in selecting bin number");
    exit(1);
  }
  if (pm > 0)
  {
    return xi_vec_plus[q];
  }
  else
  {
    return xi_vec_minus[q];
  }
}

// ---------------------------------------------------------------------------

double w_gammat_tomo(int nt, int ni, int nj, int limber)
{
  if(ni < -1 || ni > tomo.clustering_Nbin -1 || nj < -1 || nj > tomo.shear_Nbin -1)
  {
    log_fatal("invalid bin input (ni, nj) = (%d, %d)", ni, nj);
    exit(1);
  }
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

  static double** Pl = 0;
  static double* w_vec = 0;
  static cosmopara C;
  static nuisancepara N;
  static galpara G;
    
  const int nell = limits.LMAX;
  const int ntheta = like.Ntheta;
  const int NSIZE = tomo.ggl_Npowerspectra;

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

  if (recompute_gs(C, G, N))
  {
    double** Cl = malloc(NSIZE*sizeof(double*));
    for (int nz = 0; nz<NSIZE; nz++)
    {
      Cl[nz] = calloc(nell, sizeof(double));
    }
    if(limber == 1)
    {
      C_gs_tomo_limber(limits.LMIN_tab, ZL(0), ZS(0)); // init static vars
      #pragma omp parallel for collapse(2)
      for (int nz=0; nz<NSIZE; nz++)
      {
        for (int l=2; l<nell; l++)
        {
          const int ZLNZ = ZL(nz);
          const int ZSNZ = ZS(nz);
          Cl[nz][l] = (l < limits.LMIN_tab + 1) ?
           C_gs_tomo_limber_nointerp((double) l, ZLNZ, ZSNZ, use_linear_ps_limber, 0) :
           Cl[nz][l] = C_gs_tomo_limber((double) l, ZLNZ, ZSNZ);
        }
      }
    }
    else
    {
      const int L = 1;
      const double tolerance = 0.0075;    // required fractional accuracy in C(l)
      const double dev = 10. * tolerance; // will be diff exact vs Limber init to large
                                          // value in order to start while loop
      for (int nz=0; nz<NSIZE; nz++)
      { // Cocoa: no threading allowed here - (fftw allocation)
        const int Z1 = ZL(nz);
        const int Z2 = ZS(nz);
        C_gl_tomo(L, Z1, Z2, Cl[nz], dev, tolerance);
      }
      #pragma omp parallel for collapse(2)
      for (int nz=0; nz<NSIZE; nz++)
      {
        for (int l=limits.LMAX_NOLIMBER+1; l<nell; l++)
        {
          const int Z1 = ZL(nz);
          const int Z2 = ZS(nz);
          Cl[nz][l] = C_gg_tomo_limber((double) l, Z1, Z2);
        }
      }
    }
    #pragma omp parallel for collapse(2)
    for (int nz=0; nz<NSIZE; nz++)
    {
      for (int i=0; i<ntheta; i++)
      {
        const int q = nz*ntheta+i;
        w_vec[q] = 0;
        for (int l=1; l<nell; l++)
        {
          w_vec[q] += Pl[i][l]*Cl[nz][l];
        }
      }
    }
    for (int nz=0; nz<NSIZE; nz++)
    {
      free(Cl[nz]);
    }
    free(Cl);
    update_cosmopara(&C);
    update_galpara(&G);
    update_nuisance(&N);
  }

  if (!test_zoverlap(ni, nj)) 
  {
    return 0.0;
  }

  const int q = N_ggl(ni, nj)*ntheta + nt;
  if(q > NSIZE*ntheta - 1)
  {
    log_fatal("error in selecting bin number");
    exit(1);
  }

  return w_vec[q];
}

// ---------------------------------------------------------------------------

double w_gg_tomo(int nt, int ni, int nj, int limber)
{
  if(ni < -1 || ni > tomo.clustering_Nbin -1 || nj < -1 || nj > tomo.clustering_Nbin -1)
  {
    log_fatal("invalid bin input (ni, nj) = (%d, %d)", ni, nj);
    exit(1);
  }
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

  static double** Pl = 0;
  static double* w_vec = 0;
  static cosmopara C;
  static nuisancepara N;
  static galpara G;

  const int nell = limits.LMAX;
  const int ntheta = like.Ntheta;
  const int NSIZE = tomo.clustering_Npowerspectra;

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

  if (recompute_gg(C, G, N))
  {    
    double** Cl = malloc(NSIZE*sizeof(double*));
    for (int nz=0; nz<NSIZE; nz++)
    {
      Cl[nz] = calloc(nell, sizeof(double));
    }
    if(limber == 1)
    {
      C_gg_tomo_limber(limits.LMIN_tab, 0, 0); // init static vars
      #pragma omp parallel for collapse(2)
      for (int nz=0; nz<NSIZE; nz++)
      {
        for (int l=1; l<nell; l++)
        {
          const int q = nz;
          const int Z1 = nz; // cross redshift bin not supported so not using ZCL1(k)
          const int Z2 = nz; // cross redshift bin not supported so not using ZCL2(k)
          Cl[q][l] = (l < limits.LMIN_tab + 1) ?
            C_gg_tomo_limber_nointerp((double) l, Z1, Z2, use_linear_ps_limber, 0) :
            C_gg_tomo_limber((double) l, Z1, Z2);
        }
      }
    }
    else
    {
      const int L = 1;
      const double tolerance = 0.0075;    // required fractional accuracy in C(l)
      const double dev = 10. * tolerance; // will be diff  exact vs Limber init to
                                          // large value in order to start while loop
      for (int nz=0; nz<NSIZE; nz++)
      { // Cocoa: no threading allowed here - (fftw allocation)
        const int Z1 = nz; // cross redshift bin not supported so not using ZCL1(k)
        const int Z2 = nz; // cross redshift bin not supported so not using ZCL2(k)
        C_cl_tomo(L, Z1, Z2, Cl[nz], dev, tolerance);
      }
      #pragma omp parallel for collapse(2)
      for (int nz=0; nz<NSIZE; nz++)
      {
        for (int l=limits.LMAX_NOLIMBER+1; l<nell; l++)
        {
          const int q = nz;
          const int Z1 = nz; // cross redshift bin not supported so not using ZCL1(k)
          const int Z2 = nz; // cross redshift bin not supported so not using ZCL2(k)
          Cl[q][l] = C_gg_tomo_limber((double) l, Z1, Z2);
        }
      }
    }
    #pragma omp parallel for collapse(2)
    for (int nz=0; nz<NSIZE; nz++)
    {
      for (int i=0; i<ntheta; i++)
      {
        const int q = nz*ntheta + i;
        w_vec[q] = 0;
        for (int l=1; l<nell; l++)
        {
          w_vec[q] += Pl[i][l]*Cl[nz][l];
        }
      }
    }
    for (int nz=0; nz<NSIZE; nz++)
    {
      free(Cl[nz]);
    }
    free(Cl);
    update_cosmopara(&C);
    update_galpara(&G);
    update_nuisance(&N);
  }
  
  if (ni != nj)
  {
    log_fatal("ni != nj tomography not supported");
    exit(1);
  }

  const int q = ni * ntheta + nt;
  if(q > NSIZE*ntheta - 1)
  {
    log_fatal("error in selecting bin number");
    exit(1);
  }  
  return w_vec[q];
}

// ---------------------------------------------------------------------------

double w_gk_tomo(int nt, int ni, int limber)
{
  if(ni < -1 || ni > tomo.clustering_Nbin -1)
  {
    log_fatal("invalid bin input ni = %d", ni);
    exit(1);
  } 
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

  static double** Pl =0;
  static double* w_vec = 0;
  static cosmopara C;
  static nuisancepara N;
  static galpara G;
  
  const int nell = limits.LMAX;
  const int ntheta = like.Ntheta;
  const int NSIZE = tomo.clustering_Nbin;
  
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
        Pl[i][l] = 1./(4.*M_PI)*(Pmin[l+1]-Pmax[l+1]-Pmin[l-1]+Pmax[l-1])/
          (xmin[i]-xmax[i]);
      }
      free(Pmin);
      free(Pmax);
    }
  }
  
  if (recompute_gk(C, G, N))
  {
    if(limber == 1)
    {
      double** Cl = malloc(NSIZE*sizeof(double*));
      for (int nz = 0; nz<NSIZE; nz++)
      {
        Cl[nz] = calloc(nell, sizeof(double));
      }

      C_gk_tomo_limber_wrapper(limits.LMIN_tab + 1, 0); // init the function
      #pragma omp parallel for collapse(2)
      for (int nz=0; nz<NSIZE; nz++)
      {
        for (int l=1; l<nell; l++)
        {
          Cl[nz][l] = (l < limits.LMIN_tab + 1) ? 
            C_gk_tomo_limber_nointerp_wrapper((double) l, nz, use_linear_ps_limber) :
             C_gk_tomo_limber_wrapper((double) l, nz);
        }
      }
      #pragma omp parallel for collapse(2)
      for (int nz=0; nz<NSIZE; nz++)
      {
        for (int i=0; i<ntheta; i++)
        {
          w_vec[nz*ntheta+i] = 0;
          for (int l=1; l<nell; l++)
          {
            w_vec[nz*ntheta+i] += Pl[i][l]*Cl[nz][l];
          }
        }
      }
      for (int nz=0; nz<NSIZE; nz++)
      {
        free(Cl[nz]);
      }
      free(Cl);
    }
    else
    {
      log_fatal("NonLimber not implemented");
      exit(1);
    }
    update_cosmopara(&C);
    update_galpara(&G);
    update_nuisance(&N);
  }
  
  const int q = ni * ntheta + nt;
  if(q > NSIZE*ntheta - 1)
  {
    log_fatal("error in selecting bin number");
    exit(1);
  }  
  return w_vec[q];
}

// ---------------------------------------------------------------------------

double w_ks_tomo(int nt, int ni, int limber)
{
  if(ni < -1 || ni > tomo.shear_Nbin -1)
  {
    log_fatal("invalid bin input ni = %d", ni);
    exit(1);
  }
  if (like.Ntheta ==0)
  {
    log_fatal("like.Ntheta not initialized");
    exit(1);
  }
  if(nt > like.Ntheta - 1)
  {
    log_fatal("error in selecting bin number");
    exit(1); 
  }

  static double** Pl =0;
  static double* w_vec =0;
  static cosmopara C;
  static nuisancepara N;
  static galpara G;
  
  const int nell = limits.LMAX;
  const int ntheta = like.Ntheta;
  const int NSIZE = tomo.shear_Nbin;

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
  
  if (recompute_ks(C, G, N))
  {
    if(limber == 1)
    {
      double** Cl = malloc(tomo.shear_Nbin*sizeof(double*));
      for (int nz=0; nz<tomo.shear_Nbin; nz++)
      {
        Cl[nz] = calloc(nell, sizeof(double));
      }

      C_ks_tomo_limber_wrapper(limits.LMIN_tab + 1, 0.0); // init the function
      #pragma omp parallel for collapse(2)
      for (int nz=0; nz<tomo.shear_Nbin; nz++)
      {
        for (int l=1; l<nell; l++)
        {
          Cl[nz][l] = (l < limits.LMIN_tab + 1) ? 
            C_ks_tomo_limber_nointerp_wrapper(l, nz, use_linear_ps_limber) :
            Cl[nz][l] = C_ks_tomo_limber_wrapper(l, nz);
        }
      }
      #pragma omp parallel for collapse(2)
      for (int nz=0; nz<tomo.shear_Nbin; nz++)
      {
        for (int i=0; i<ntheta; i++)
        {
          w_vec[nz*ntheta+i] = 0;
          for (int l=2; l<nell; l++)
          {
            w_vec[nz*ntheta+i] += Pl[i][l]*Cl[nz][l];
          }
        }
      }
      for (int nz=0; nz<tomo.shear_Nbin; nz++)
      {
        free(Cl[nz]);
      }
      free(Cl);
    } 
    else
    {
      log_fatal("NonLimber not implemented");
      exit(1);
    }
    update_cosmopara(&C);
    update_galpara(&G);
    update_nuisance(&N);
  }

  const int q = ni * ntheta + nt;
  if(q > NSIZE*ntheta - 1)
  {
    log_fatal("error in selecting bin number");
    exit(1);
  }  
  return w_vec[q];
}

// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// Correlation Functions (real space) - flat sky
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------

double xi_pm_tomo_flatsky(int pm, double theta, int ni, int nj, int limber)
{
  if(ni < -1 || ni > tomo.shear_Nbin -1 || nj < -1 || nj > tomo.shear_Nbin -1)
  {
    log_fatal("invalid bin input (ni, nj) = (%d, %d)", ni, nj);
    exit(1);
  }

  static cosmopara C;
  static nuisancepara N;
  static double** table;

  const int NSIZE = tomo.shear_Npowerspectra;
  const int ntheta = Ntable.N_thetaH;
  const double l_min = w_l_min;
  const double l_max = w_l_max;
  const double lnlmax = log(l_max);
  const double lnlmin = log(l_min);
  const double dlnl = (lnlmax - lnlmin)/(1.0*ntheta);
  const double lnrc = 0.5*(lnlmax + lnlmin);
  const double nc = ntheta/2 + 1;
  
  const double lnthetamin = (nc-ntheta+1)*dlnl-lnrc;
  const double lnthetamax = nc*dlnl-lnrc;
  const double dlntheta = (lnthetamax - lnthetamin)/((double) ntheta);
  const double lntheta = log(theta);

  if (table == 0)
  {
    table = (double**) malloc(sizeof(double*)*2*NSIZE);  // 2 NSIZE = {xi+, xi-}
    for(int i=0; i<NSIZE; i++) 
    {
      table[i] = (double*) malloc(sizeof(double)*ntheta);
    }
  }
  if (recompute_shear(C, N))
  {
    typedef fftw_complex fftwZ;

    fftwZ** flP = (fftwZ**) malloc(sizeof(fftwZ*)*NSIZE);
    for (int j=0; j<NSIZE; j++)
    {
      flP[j] = (fftwZ*) fftw_malloc((Ntable.N_thetaH/2+1)*sizeof(fftwZ));
    }
    {
      double** lP = (double**) malloc(sizeof(double*)*NSIZE);
      fftw_plan* plan = (fftw_plan*) malloc(sizeof(fftw_plan)*NSIZE);
      for (int j=0; j<NSIZE; j++)
      {
        lP[j] = (double*) malloc(ntheta*sizeof(double));
        plan[j] = fftw_plan_dft_r2c_1d(ntheta,lP[j],flP[j],FFTW_ESTIMATE);
      }
      
      // Power spectrum on logarithmic bins (begins)        
      {
        if (like.IA == 5 || like.IA == 6)
        { // NEW TATT MODELING
          log_fatal("Limber && TATT not implemented");
          exit(1);
        }
        else
        {
          C_ss_tomo_limber(limits.LMIN_tab+1, Z1(0), Z2(0)); // init the function
          #pragma omp parallel for collapse(2)
          for (int k=0; k<NSIZE; k++)
          {
            for(int p=0; p<ntheta; p++)
            {
              const int Z1NZ = Z1(k);
              const int Z2NZ = Z2(k);
              const double l = exp(lnrc+(p - nc)*dlnl);
              lP[k][p] = (l > limits.LMIN_tab) ? l*C_ss_tomo_limber(l, Z1NZ, Z2NZ) :
                l*C_ss_tomo_limber_nointerp(l, Z1NZ, Z2NZ, use_linear_ps_limber);
            }
          }
        }
      } // Power spectrum on logarithmic bins (ends) 
      
      #pragma omp parallel for
      for (int j=0; j<NSIZE; j++)
      {
        fftw_execute(plan[j]); // Execute FFTW in parallel (thread-safe)
      }
      
      for (int j = 0; j < NSIZE; j++)
      {
        fftw_free(lP[j]);
        fftw_destroy_plan(plan[j]);
      }
      free(lP);
      free(plan);
    }
    
    double*** lP = (double***) malloc(sizeof(double**)*NSIZE);
    fftwZ*** kernel = (fftwZ***) malloc(sizeof(fftwZ**)*NSIZE);
    fftwZ*** conv = (fftwZ***) malloc(sizeof(fftwZ**)*NSIZE);
    fftw_plan** plan = (fftw_plan**) malloc(sizeof(fftw_plan*)*NSIZE);
    double*** tab = (double***) malloc(sizeof(double**)*NSIZE);
    for (int j=0; j<NSIZE; j++)
    {
      lP[j] = (double**) malloc(sizeof(double*)*2);
      kernel[j] = (fftwZ**) malloc(sizeof(fftwZ*)*2);
      conv[j] = (fftwZ**) malloc(sizeof(fftwZ*)*2);
      plan[j] = (fftw_plan*) malloc(sizeof(fftw_plan)*2);
      tab[j] = (double**) malloc(sizeof(double*)*1);
      for(int i=0; i<1; i++) 
      {
        tab[j][i] = (double*) malloc(sizeof(double)*ntheta);
      }
      for (int m=0; m<2; m++)
      {
        const int COVSZ = (ntheta/2+1);
        lP[j][m] = (double*) malloc(ntheta*sizeof(double));
        kernel[j][m] = (fftwZ*) fftw_malloc(sizeof(fftwZ));
        conv[j][m] = (fftwZ*) fftw_malloc(COVSZ*sizeof(fftwZ));
        plan[j][m] = fftw_plan_dft_c2r_1d(ntheta,conv[j][m],lP[j][m],FFTW_ESTIMATE);
      }
    }
    
    #pragma omp parallel for
    for (int j=0; j<NSIZE; j++)
    {
      for (int m=0; m<2; m++)
      {
        double arg[2];
        arg[0] = 0; // bias
        arg[1] = (m == 0 ? 0 : 4); // order of Bessel function

        // perform the convolution, negative sign for kernel (complex conj.)
        for(int i=0; i<(ntheta/2+1); i++)
        {
          const double k = 2*M_PI*i/(dlnl*ntheta);
          hankel_kernel_FT(k, kernel[j][m], arg, 2);
          conv[j][m][i][0] = flP[j][i][0]*(kernel[j][m][0][0])-flP[j][i][1]*(kernel[j][m][0][1]);
          conv[j][m][i][1] = flP[j][i][1]*(kernel[j][m][0][0])+flP[j][i][0]*(kernel[j][m][0][1]);
        }

        // force Nyquist- and 0-frequency-components to be double
        conv[j][m][0][1] = 0;
        conv[j][m][ntheta/2][1] = 0;

        fftw_execute(plan[j][m]);

        for(int k=0; k<ntheta; k++)
        {
          const double t = exp((nc-k)*dlnl-lnrc); // theta=1/l
          tab[j][m][ntheta-k-1] = lP[j][m][k]/(t*2*M_PI*ntheta);
        }
      }
      for (int k=0; k<ntheta; k++)
      {
        table[2*j][k] = tab[j][0][k];
        table[2*j+1][k] = tab[j][1][k];
      }
    }

    for (int j=0; j<NSIZE; j++)
    {
      for (int m=0; m<2; m++)
      {
        fftw_free(lP[j][m]);
        fftw_free(kernel[j][m]);
        fftw_free(conv[j][m]);
        fftw_destroy_plan(plan[j][m]);
      }
      fftw_free(flP[j]);
      free(lP[j]);
      free(kernel[j]);
      free(conv[j]);
      free(plan[j]);
      free(tab[j]);
    }
    free(flP);
    free(lP);
    free(kernel);
    free(conv);
    free(plan);
    free(tab);
    // --------------------------------------------------------------------
    // Cocoa: code extracted (& adapted) from xipm_via_hankel (ends)
    // --------------------------------------------------------------------
  }
  update_cosmopara(&C);
  update_nuisance(&N);

  const int q = 2*N_shear(ni, nj) + (1 - pm)/2;
  if(q > NSIZE - 1)
  {
    log_fatal("error in selecting bin number");
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

  return interpol(table[q], ntheta, lnthetamin, lnthetamax, dlntheta, lntheta, 0, 0);
}

// ---------------------------------------------------------------------------

double w_gammat_tomo_flatsky(double theta, int ni, int nj, int limber)
{
  static cosmopara C;
  static nuisancepara N;
  static galpara G;
  static double** table;

  const int NSIZE = tomo.ggl_Npowerspectra;
  const int ntheta = Ntable.N_thetaH;

  const double l_min = w_l_min;
  const double l_max = w_l_max;
  const double lnlmax = log(l_max);
  const double lnlmin = log(l_min);
  const double dlnl = (lnlmax-lnlmin)/(ntheta - 1.);
  const double lnrc = 0.5*(lnlmax+lnlmin);
  const double nc = ntheta/2 + 1;
  
  const double lnthetamin = (nc-ntheta + 1)*dlnl-lnrc;
  const double lnthetamax = nc*dlnl-lnrc;
  const double dlntheta = (lnthetamax - lnthetamin)/((double)ntheta);
  const double lntheta = log(theta);

  if (table == 0)
  {
    table = (double**) malloc(sizeof(double*)*NSIZE);
    for(int i=0; i<NSIZE; i++) 
    {
      table[i] = (double*) malloc(sizeof(double)*ntheta);
    }
  }
  if (recompute_gs(C, G, N)) 
  {
    typedef fftw_complex fftwZ;

    fftwZ** flP = (fftwZ**) malloc(sizeof(fftwZ*)*NSIZE); // go to log-Fourier-space
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

      { // Power spectrum on logarithmic bins (begins) 
        gsl_spline** fCL_NL = NULL; 
        if (limber != 1)
        { 
          fCL_NL = (gsl_spline**) malloc(sizeof(gsl_spline*)*NSIZE);
          double** Cl_NL = (double**) malloc(sizeof(double*)*NSIZE);
          for (int i=0; j<NSIZE; j++)
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
          for (int k=0; k<NSIZE; k++)
          { // Cocoa: no threading allowed here - (fftw allocation @C_gl_tomo)
            const int Z1 = ZL(k);
            const int Z2 = ZS(k);
            C_gl_tomo(L, Z1, Z2, Cl_NL[k], dev, tolerance);

            const gsl_interp_type* T = gsl_interp_linear;
            fCL_NL[k] = gsl_spline_alloc(T, limits.LMAX_NOLIMBER);
            if (fCL_NL[k] == NULL)
            {
              log_fatal("fail allocation");
              exit(1);
            }
          }
          #pragma omp parallel for
          for (int k=0; k<NSIZE; k++)
          {
            int status = gsl_spline_init(fCL_NL[k], ll, Cl_NL[k], limits.LMAX_NOLIMBER);
            if (status) 
            {
              log_fatal(gsl_strerror(status));
              exit(1);
            }
          }
          for (int i=0; j<NSIZE; j++)
          {
            free(Cl_NL[j]);
          }
          free(Cl_NL);
          free(ll);
        }
        C_gs_tomo_limber(limits.LMIN_tab, ZL(0), ZS(0)); // init static vars
        #pragma omp parallel for collapse(2)
        for (int k=0; k<NSIZE; k++)
        {
          for(int p=0; p<ntheta; p++)
          {
            const int ZLNZ = ZL(k);
            const int ZSNZ = ZS(k);
            const double l = exp(lnrc + (p - nc)*dlnl);
            if (limber == 1 || (limber != 1 && l > limits.LMAX_NOLIMBER - 1))
            {
              lP[q][p] = (l > limits.LMIN_tab) ? l*C_gs_tomo_limber(l, ZLNZ, ZSNZ) :
                l*C_gs_tomo_limber_nointerp(l, ZLNZ, ZSNZ, use_linear_ps_limber, 0);
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
              lP[k][p] = l*CLNL;
            }
          }
        }
        if (limber != 1)
        {
          for (int i=0; i<NSIZE; i++)
          {
             gsl_spline_free(fCL_NL[i]);
          }
          free(fCL_NL);
        }
      } // Power spectrum on logarithmic bins (ends)
      
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
      const int COVSZ = (ntheta/2+1);
      lP[j] = (double*) malloc(ntheta*sizeof(double));
      kernel[j] = (fftwZ*) fftw_malloc(sizeof(fftwZ));
      conv[j] = (fftwZ*) fftw_malloc(COVSZ*sizeof(fftwZ));
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
        const double t = exp((nc-k)*dlnl-lnrc); // theta=1/l
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
    update_galpara(&G);
    update_nuisance(&N);
  }

  if(ni < -1 || ni > tomo.clustering_Nbin -1 || nj < -1 || nj > tomo.shear_Nbin -1)
  {
    log_fatal("invalid bin input (ni, nj) = (%d, %d)", ni, nj);
    exit(1);
  }

  if(test_zoverlap(ni, nj))
  {
    const int q = N_ggl(ni, nj);
    if(q > NSIZE - 1)
    {
      log_fatal("error in selecting bin number");
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

    return interpol(table[q], ntheta, lnthetamin, lnthetamax, dlntheta, lntheta, 0, 0);
  } 
  else 
  {
    return 0.0;
  }
}

// ---------------------------------------------------------------------------

double w_gg_tomo_flatsky(double theta, int ni, int nj, int limber)
{
  static cosmopara C;
  static nuisancepara N;
  static galpara G;
  static double** table;

  const int NSIZE = tomo.clustering_Npowerspectra;
  const int ntheta = Ntable.N_thetaH;

  const double lmin = w_l_min;
  const double lmax = w_l_max;
  const double lnlmax = log(lmax);
  const double lnlmin = log(lmin);
  const double dlnl = (lnlmax - lnlmin)/(1.0*ntheta - 1.);
  const double lnrc = 0.5*(lnlmax + lnlmin);
  const double nc = ntheta/2 + 1;

  const double lntheta = log(theta);
  const double lnthetamin = (nc-ntheta+1)*dlnl-lnrc;
  const double lnthetamax = nc*dlnl-lnrc;
  const double dlntheta = (lnthetamax-lnthetamin)/((double) ntheta);

  if (table == 0)
  {
    table = (double**) malloc(sizeof(double*)*NSIZE);
    for(int i=0; i<NSIZE; i++) 
    {
      table[i] = (double*) malloc(sizeof(double)*ntheta);
    }
  }  
  if (recompute_gg(C, G, N)) 
  {
    typedef fftw_complex fftwZ;

    fftwZ** flP = (fftwZ**) malloc(sizeof(fftwZ*)*NSIZE); // go to log-Fourier-space
    for (int j=0; j<NSIZE; j++)
    {
      flP[j] = fftw_malloc((ntheta/2+1)*sizeof(fftwZ));
    }
    {
      double** lP = (double**) malloc(sizeof(double*)*NSIZE);
      fftw_plan* plan = (fftw_plan*) malloc(sizeof(fftw_plan)*NSIZE);
      for (int j=0; j<NSIZE; j++)
      {
        int ARRAYSZ = ntheta;
        lP[j] = (double*) malloc(ARRAYSZ*sizeof(double));
        plan[j] = fftw_plan_dft_r2c_1d(ARRAYSZ, lP[j], flP[j], FFTW_ESTIMATE);
      }

      { // Power spectrum on logarithmic bins (begins)
        gsl_spline** fCL_NL = NULL; 
        if (limber != 1)
        { 
          fCL_NL = (gsl_spline**) malloc(sizeof(gsl_spline*)*NSIZE);
          double** Cl_NL = (double**) malloc(sizeof(double*)*NSIZE);
          for (int i=0; j<NSIZE; j++)
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
          for (int k=0; k<NSIZE; k++)
          { // Cocoa: no threading allowed here - (fftw allocation @C_gl_tomo)
            const int ZCL1 = k; // cross redshift bin not supported so not using ZCL1(k)
            const int ZCL2 = k; // cross redshift bin not supported so not using ZCL2(k) 
            C_cc_tomo(L, ZCL1, ZCL1, Cl_NL[k], dev, tolerance);

            const gsl_interp_type* T = gsl_interp_linear;
            fCL_NL[k] = gsl_spline_alloc(T, limits.LMAX_NOLIMBER);
            if (fCL_NL[k] == NULL)
            {
              log_fatal("fail allocation");
              exit(1);
            }
          }
          #pragma omp parallel for
          for (int k=0; k<NSIZE; k++)
          {
            int status = gsl_spline_init(fCL_NL[k], ll, Cl_NL[k], limits.LMAX_NOLIMBER);
            if (status) 
            {
              log_fatal(gsl_strerror(status));
              exit(1);
            }
          } 
          for (int i=0; j<NSIZE; j++)
          {
            free(Cl_NL[j]);
          }
          free(Cl_NL);
          free(ll);
        } 

        C_gg_tomo_limber(limits.LMIN_tab, 0, 0); // init static vars
        #pragma omp parallel for collapse(2)
        for (int k=0; k<NSIZE; k++) 
        { 
          for(int p=0; p<ntheta; p++)
          {
            const int ZCL1 = k; // cross redshift bin not supported so not using ZCL1(k)
            const int ZCL2 = k; // cross redshift bin not supported so not using ZCL2(k) 
            const double l = exp(lnrc + (p - nc)*dlnl);
            if (limber == 1 || (limber != 1 && l > limits.LMAX_NOLIMBER - 1))
            {
              lP[k][p] = (l > limits.LMIN_tab) ? l*C_gg_tomo_limber(l, ZCL1, ZCL2) :
                l*C_gg_tomo_limber_nointerp(l, ZCL1, ZCL2, use_linear_ps_limber, 0);
            }
            else
            {
              double CLNL;
              int status = gsl_spline_eval_e(fCL_NL[k], l, NULL, NULL, &CLNL);
              if (status) 
              {
                log_fatal(gsl_strerror(status));
                exit(1);
              }
              lP[k][p] = l*CLNL;
            }
          }
        }
        if (limber != 1)
        {
          for (int i=0; i<NSIZE; i++)
          {
             gsl_spline_free(fCL_NL[i]);
          }
          free(fCL_NL);
        }
      } // Power spectrum on logarithmic bins (ends)

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
      const int ARRAYSZ = ntheta;
      const int COVSZ = (ntheta/2+1);
      lP[j] = (double*) malloc(ARRAYSZ*sizeof(double));
      kernel[j] = (fftwZ*) fftw_malloc(sizeof(fftwZ));
      conv[j] = (fftwZ*) fftw_malloc(COVSZ*sizeof(fftwZ));
      plan[j] = fftw_plan_dft_c2r_1d(ARRAYSZ, conv[j], lP[j], FFTW_ESTIMATE);
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
        const double t = exp((nc-k)*dlnl-lnrc); // theta=1/l
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
    update_galpara(&G);
    update_nuisance(&N);
  }
 
  if(ni < -1 || ni > tomo.clustering_Nbin -1 || nj < -1 || nj > tomo.clustering_Nbin -1)
  {
    log_fatal("invalid bin input (ni, nj) = (%d, %d)", ni, nj);
    exit(1);
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

  const int q = ni; // cross redshift bin not supported so not using N_CL(ni, nj) instead of ni
  if(q > NSIZE - 1)
  {
    log_fatal("error in selecting bin number");
    exit(1);
  }

  return interpol(table[q], ntheta, lnthetamin, lnthetamax, dlntheta, lntheta, 1, 1);
}

// ---------------------------------------------------------------------------

double w_gk_tomo_flatsky(double theta, int ni, int limber)
{
  static cosmopara C;
  static nuisancepara N;
  static galpara G;
  static double** table;

  const int nisize = tomo.clustering_Nbin;
  const int ntheta = Ntable.N_thetaH;
  const int NSIZE = nisize;

  const double lmin = w_l_min;
  const double lmax = w_l_max;
  const double lnlmax = log(lmax);
  const double lnlmin = log(lmin);
  const double dlnl = (lnlmax - lnlmin)/(1.0*ntheta - 1.);
  const double lnrc = 0.5*(lnlmax + lnlmin);
  const double nc = ntheta/2 + 1;
  
  const double lntheta = log(theta);
  const double lnthetamin = (nc-ntheta+1)*dlnl-lnrc;
  const double lnthetamax = nc*dlnl-lnrc;
  const double dlntheta = (lnthetamax-lnthetamin)/((double)ntheta);

  if (table == 0)
  {
    table = (double**) malloc(sizeof(double*)*NSIZE);
    for(int i=0; i<NSIZE; i++) 
    {
      table[i] = (double*) malloc(sizeof(double)*ntheta);
    }
  }
  if (recompute_gk(C, G, N))
  {
    if (limber == 1)
    {
      typedef fftw_complex fftwZ;

      // go to log-Fourier-space
      fftwZ** flP = (fftwZ**) malloc(sizeof(fftwZ*)*NSIZE);
      for (int i=0; i<NSIZE; i++)
      {
        flP[i] = fftw_malloc((ntheta/2+1)*sizeof(fftwZ));
      }
      {
        double** lP = (double**) malloc(sizeof(double*)*NSIZE);
        fftw_plan* plan = (fftw_plan*) malloc(sizeof(fftw_plan)*NSIZE);
        for (int i=0; i<NSIZE; i++)
        {
          lP[i] = (double*) malloc(ntheta*sizeof(double));
          plan[i] = fftw_plan_dft_r2c_1d(ntheta, lP[i], flP[i], FFTW_ESTIMATE);
        }

        { // Power spectrum on logarithmic bins (begins)
          C_gk_tomo_limber_wrapper(limits.LMIN_tab + 1, 0); // init the function
          #pragma omp parallel for collapse(2)
          for (int i=0; i<NSIZE; i++)
          {
            for(int p=0; p<ntheta; p++)
            {
              const double l = exp(lnrc + (p - nc)*dlnl);
              P[i][p] = (l > limits.LMIN_tab) ? l*C_gk_tomo_limber_wrapper(l, i) :
                l*C_gk_tomo_limber_nointerp_wrapper(l, i, use_linear_ps_limber);
            }
          }
        } // Power spectrum on logarithmic bins (ends)
        
        #pragma omp parallel for
        for (int i=0; i<NSIZE; i++)
        {
          fftw_execute(plan[i]); // Execute FFTW in parallel (thread-safe)
        }
        
        for (int i=0; i<NSIZE; i++)
        {
          fftw_free(lP[i]);
          fftw_destroy_plan(plan[i]);
        } 
        free(lP);
        free(plan);
      }

      double** lP = (double**) fftw_malloc(sizeof(double**)*NSIZE);
      fftwZ** kernel = (fftwZ**) malloc(sizeof(fftwZ*)*NSIZE);
      fftwZ** conv = (fftwZ**) malloc(sizeof(fftwZ*)*NSIZE);
      fftw_plan* plan = (fftw_plan*) malloc(sizeof(fftw_plan)*NSIZE);
      double*** tab = (double***) malloc(sizeof(double**)*NSIZE);
      for (int i=0; i<NSIZE; i++)
      {
        lP[i] = (double*) fftw_malloc(ntheta*sizeof(double));
        kernel[i] = (fftwZ*) fftw_malloc(sizeof(fftwZ));
        conv[i] = (fftwZ*) fftw_malloc((ntheta/2+1)*sizeof(fftwZ));
        plan[i] = fftw_plan_dft_c2r_1d(ntheta, conv[i], lP[i], FFTW_ESTIMATE);
        tab[i] = (double**) malloc(sizeof(double*)*1);
        for(int p=0; p<1; p++) 
        {
          tab[i][p] = (double*) malloc(sizeof(double)*ntheta);
        }
      }
      
      #pragma omp parallel for
      for (int i=0; i<NSIZE; i++)
      {
        double arg[2];
        arg[0] = 0; // bias
        arg[1] = 0; // order of Bessel function

        for(int k=0; k<(ntheta/2+1); k++)
        {
          const double kk = 2*M_PI*i/(dlnl*ntheta);
          hankel_kernel_FT(kk, kernel[i], arg, 2);
          
          conv[i][k][0] = flP[i][k][0]*kernel[i][0][0] - flP[i][k][1]*kernel[i][0][1];
          conv[i][k][1] = flP[i][k][1]*kernel[i][0][0] + flP[i][k][0]*kernel[i][0][1];
        }

        // force Nyquist- and 0-frequency-components to be double
        conv[i][0][1] = 0;
        conv[i][ntheta/2][1] = 0;

        fftw_execute(plan[i]); // Execute FFTW in parallel (thread-safe)

        for(int k=0; k<ntheta; k++)
        {
          const double t = exp((nc-k)*dlnl-lnrc); // theta=1/l
          tab[i][0][ntheta-k-1] = lP[i][k]/(t*2*M_PI*ntheta);
        }
        for (int k=0; k<ntheta; k++)
        {
          table[i][k] = tab[i][0][k];
        }
      }
      
      for (int i=0; i<NSIZE; i++)
      {
        fftw_free(flP[i]);
        fftw_free(lP[i]);
        fftw_free(kernel[i]);
        fftw_free(conv[i]);
        fftw_destroy_plan(plan[i]);
        free(tab[i]);
      }
      free(flP);
      free(lP);
      free(kernel);
      free(conv);
      free(plan);
      free(tab);
    }
    else
    {
      log_fatal("NonLimber not implemented");
      exit(1);
    }
    update_cosmopara(&C);
    update_galpara(&G);
    update_nuisance(&N);
  }

  if(ni < -1 || ni > tomo.clustering_Nbin -1)
  {
    log_fatal("invalid bin input ni = %d", ni);
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

  return interpol(table[ni], ntheta, lnthetamin, lnthetamax, dlntheta, lntheta, 1, 1);
}

// ---------------------------------------------------------------------------

double w_ks_tomo_flatsky(double theta, int ni, int limber)
{
  static cosmopara C;
  static nuisancepara N;
  static double** table;

  const int nisize = tomo.shear_Nbin;
  const int ntheta = Ntable.N_thetaH;
  const int NSIZE = nisize;

  const double l_min = w_l_min;
  const double l_max = w_l_max;
  const double lnlmax = log(l_max);
  const double lnlmin = log(l_min);
  const double dlnl = (lnlmax-lnlmin)/(ntheta - 1.);
  const double lnrc = 0.5*(lnlmax+lnlmin);
  const double nc = ntheta/2 + 1;
  
  const double lnthetamin = (nc-ntheta+1)*dlnl-lnrc;
  const double lnthetamax = nc*dlnl-lnrc;
  const double dlntheta = (lnthetamax-lnthetamin)/((double) ntheta);
  const double lntheta = log(theta);

  if (table == 0) 
  {
    table = (double**) malloc(sizeof(double*)*NSIZE);
    for(int i=0; i<NSIZE; i++) 
    {
      table[i] = (double*) malloc(sizeof(double)*ntheta);
    }
  }
  if (recompute_shear(C, N))
  {
    if (limber == 1)
    {
      typedef fftw_complex fftwZ;

      // go to log-Fourier-space
      fftwZ** flP = (fftwZ**) malloc(sizeof(fftwZ*)*NSIZE);
      for (int i=0; i<NSIZE; i++)
      {
        flP[i] = fftw_malloc((ntheta/2+1)*sizeof(fftwZ));
      }
      {
        double** lP = (double**) malloc(sizeof(double*)*NSIZE);
        fftw_plan* plan = (fftw_plan*) malloc(sizeof(fftw_plan)*NSIZE);
        for (int i=0; i<NSIZE; i++)
        {
          lP[i] = (double*) malloc(ntheta*sizeof(double));
          plan[i] = fftw_plan_dft_r2c_1d(ntheta, lP[i], flP[i], FFTW_ESTIMATE);
        }
        
        { // Power spectrum on logarithmic bins (begins)
          C_ks_tomo_limber_wrapper(limits.LMIN_tab + 1, 0) // init the function;
          #pragma omp parallel for collapse(2)
          for (int i=0; i<NSIZE; i++)
          {
            for(int p=0; p<ntheta; p++)
            {
              const double l = exp(lnrc + (p - nc)*dlnl);
              lP[i][p] = (l > limits.LMIN_tab) ? l*C_ks_tomo_limber_wrapper(l, i) :
                l*C_ks_tomo_limber_nointerp_wrapper(l, i, use_linear_ps_limber);
            }
          }
        } // Power spectrum on logarithmic bins (ends)
        
        #pragma omp parallel for
        for (int i=0; i<NSIZE; i++)
        {
          fftw_execute(plan[i]); // Execute FFTW in parallel (thread-safe)
        }
        
        for (int i=0; i<NSIZE; i++)
        {
          fftw_free(lP[i]);
          fftw_destroy_plan(plan[i]);
        }
        free(lP);
        free(plan);
      }

      double** lP = (double**) fftw_malloc(sizeof(double**)*NSIZE);
      fftwZ** kernel = (fftwZ**) malloc(sizeof(fftwZ*)*NSIZE);
      fftwZ** conv = (fftwZ**) malloc(sizeof(fftwZ*)*NSIZE);
      fftw_plan* plan = (fftw_plan*) malloc(sizeof(fftw_plan)*NSIZE);
      double*** tab = (double***) malloc(sizeof(double**)*NSIZE);
      for (int i=0; i<NSIZE; i++)
      {
        lP[i] = (double*) fftw_malloc(ntheta*sizeof(double));
        kernel[i] = (fftwZ*) fftw_malloc(sizeof(fftwZ));
        conv[i] = (fftwZ*) fftw_malloc((ntheta/2+1)*sizeof(fftwZ));
        plan[i] = fftw_plan_dft_c2r_1d(ntheta, conv[i], lP[i], FFTW_ESTIMATE);
        tab[i] = (double**) malloc(sizeof(double*)*1);
        for(int p=0; p<1; p++) 
        {
          tab[i][p] = (double*) malloc(sizeof(double)*ntheta);
        }
      }
      
      #pragma omp parallel for
      for (int i=0; i<NSIZE; i++)
      {
        double arg[2];
        arg[0] = 0; // bias
        arg[1] = 0; // order of Bessel function

        for(int k=0; k<(ntheta/2+1); k++)
        {
          const double kk = 2*M_PI*i/(dlnl*ntheta);
          hankel_kernel_FT(kk, kernel[i], arg, 2);
          conv[i][k][0] = flP[i][k][0]*kernel[i][0][0] - flP[i][k][1]*kernel[i][0][1];
          conv[i][k][1] = flP[i][k][1]*kernel[i][0][0] + flP[i][k][0]*kernel[i][0][1];
        }

        // force Nyquist- and 0-frequency-components to be double
        conv[i][0][1] = 0;
        conv[i][ntheta/2][1] = 0;

        fftw_execute(plan[i]); // Execute FFTW in parallel (thread-safe)

        for(int k=0; k<ntheta; k++)
        {
          const double t = exp((nc-k)*dlnl-lnrc); // theta=1/l
          tab[i][0][ntheta-k-1] = lP[i][k]/(t*2*M_PI*ntheta);
        }
        for (int k=0; k<ntheta; k++)
        {
          table[i][k] = tab[i][0][k];
        }
      }
      
      for (int i=0; i<NSIZE; i++)
      {
        fftw_free(flP[i]);
        fftw_free(lP[i]);
        fftw_free(kernel[i]);
        fftw_free(conv[i]);
        fftw_destroy_plan(plan[i]);
        free(tab[i]);
      }  
      free(flP);
      free(lP);
      free(kernel);
      free(conv);
      free(plan);
      free(tab);
      }
    else
    {
      log_fatal("NonLimber not implemented");
      exit(1);
    }
    update_cosmopara(&C);
    update_nuisance(&N);
  }

  if(ni < -1 || ni > tomo.shear_Nbin -1)
  {
    log_fatal("invalid bin input ni = %d", ni);
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

  return interpol(table[ni], ntheta, lnthetamin, lnthetamax, dlntheta, lntheta, 1, 1);
}


// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// Limber Approximation (Angular Power Spectrum)
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
// SS ANGULAR CORRELATION FUNCTION - TATT
// -----------------------------------------------------------------------------

// NLA/TA amplitude C1, nz argument only need if per-bin amplitude
double C1_TA(double a, int nz, double growfac_a)
{
  // per-bin IA parameters
  if (like.IA == 3 || like.IA == 5)
  {
    return -nuisance.A_z[nz]*cosmology.Omega_m*nuisance.c1rhocrit_ia/ growfac_a;
  }
  // power law evolution
  return -cosmology.Omega_m * nuisance.c1rhocrit_ia /
         growfac_a * nuisance.A_ia *
         pow(1. / (a * nuisance.oneplusz0_ia), nuisance.eta_ia);
}

// TA source bias parameter, nz argument only need if per-bin amplitude
double b_TA(double a __attribute__((unused)), int nz)
{
  // per-bin IA parameters
  if (like.IA == 5) 
  {
    return nuisance.b_ta_z[nz];
  }
  // power law evolution
  return nuisance.b_ta_z[0];
}

// TT amplitude C2, nz argument only need if per-bin amplitude
double C2_TT(double a, int nz, double growfac_a)
{
  // per-bin IA parameters
  if (like.IA == 5)
  {
    return 5. * nuisance.A2_z[nz] * cosmology.Omega_m *
           nuisance.c1rhocrit_ia * pow(1.0/growfac_a, 2.0);
  }
  // power law evolution
  return 5. * nuisance.A2_ia * cosmology.Omega_m * nuisance.c1rhocrit_ia *
         (1.0 /(growfac_a*growfac_a)) *
         pow(1. / (a * nuisance.oneplusz0_ia), nuisance.eta_ia_tt);
}

double int_for_C_ss_tomo_TATT_EE_limber(double a, void* params)
{
  if(!(a>0) || !(a<1)) 
  {
    log_fatal("a>0 and a<1 not true");
    exit(1);
  }
  double* ar = (double* ) params;

  const double growfac_a = growfac(a);
  struct chis chidchi = chi_all(a);
  double hoverh0 = hoverh0v2(a, chidchi.dchida);

  const double ell = ar[2] + 0.5;
  const double fK = f_K(chidchi.chi);
  const double k = ell / fK;

  // radial n_z weight for first source bin (for use with IA term)
  const double ws1 = W_source(a, ar[0], hoverh0);

  // radial n_z weight for second source bin (for use with IA term)
  const double ws2 = W_source(a, ar[1], hoverh0);
  // radial lens efficiency for first source bin
  const double wk1 = W_kappa(a, fK, ar[0]);
  // radial lens efficiency for second source bin
  const double wk2 = W_kappa(a, fK, ar[1]);

  // IA parameters for first source bin
  const double C1 = C1_TA(a, (int) ar[0], growfac_a);
  const double b_ta = b_TA(a, (int) ar[0]);
  const double C2 = C2_TT(a, (int) ar[0], growfac_a);

  // IA parameters for second source bin
  const double C1_2 = C1_TA(a, (int) ar[1], growfac_a);
  const double b_ta_2 = b_TA(a,(int) ar[1]);
  const double C2_2 = C2_TT(a, (int) ar[1], growfac_a);

  // GG cosmic shear
  const double pdelta_ak = Pdelta(k, a);
  double res = wk1 * wk2 * pdelta_ak;
  
  //COCOA: Took these evaluations of the parenthesis - to force them to update
  //COCOA: the static variables in the first call that is done outside OpenMP loop
  const double tmp1 = TATT_II_EE(k, a, C1, C2, b_ta, C1_2, C2_2, b_ta_2, growfac_a, pdelta_ak);
  const double tmp2 = TATT_GI_E(k, a, C1, C2, b_ta, growfac_a, pdelta_ak);
  const double tmp3 = TATT_GI_E(k, a, C1_2, C2_2, b_ta_2, growfac_a, pdelta_ak);
  if (C1 || C1_2 || C2 || C2_2) 
  {
    // II contribution
    res += ws1 * ws2 * tmp1;
    // GI contribution
    res += ws1 * wk2 * tmp2 + ws2 * wk1 * tmp3;
  }
  return res * chidchi.dchida / (fK * fK);
}

double int_for_C_ss_tomo_TATT_BB_limber(double a, void* params)
{
  if(!(a>0) || !(a<1)) 
  {
    log_fatal("a>0 and a<1 not true");
    exit(1);
  }
  double* ar = (double* ) params;

  const double growfac_a = growfac(a);
  struct chis chidchi = chi_all(a);
  const double hoverh0 = hoverh0v2(a, chidchi.dchida);

  const double ell = ar[2] + 0.5;
  const double fK = f_K(chidchi.chi);
  const double k = ell / fK;

  // radial n_z weight for first source bin (for use with IA term)
  const double ws1 = W_source(a, ar[0], hoverh0);

  // radial n_z weight for second source bin (for use with IA term)
  const double ws2 = W_source(a, ar[1], hoverh0);

  // IA parameters for first source bin
  const double C1 = C1_TA(a, (int) ar[0], growfac_a);
  const double b_ta = b_TA(a, (int) ar[0]);
  const double C2 = C2_TT(a, (int) ar[0], growfac_a);

  // IA parameters for second source bin
  const double C1_2 = C1_TA(a, (int) ar[1], growfac_a);
  const double b_ta_2 = b_TA(a, (int) ar[1]);
  const double C2_2 = C2_TT(a, (int) ar[1], growfac_a);

  //COCOA: Took these evaluations of the parenthesis - to force them to update
  //COCOA: the static variables in the first call that is done outside OpenMP loop
  double res = 0.;
  const double tmp1 = TATT_II_BB(k, a, C1, C2, b_ta, C1_2, C2_2, b_ta_2, growfac_a);
  if ((b_ta || C2) && (b_ta_2 || C2_2))
  {
    res = ws1 * ws2 * tmp1;
  }
  return res * chidchi.dchida / (fK * fK);
}


double C_ss_tomo_TATT_EE_limber_nointerp(double l, int ni, int nj)
{
  if(ni < -1 || ni > tomo.shear_Nbin -1 || nj < -1 || nj > tomo.shear_Nbin -1)
  {
    log_fatal("invalid bin input (ni, nj) = (%d, %d)", ni, nj);
    exit(1);
  }

  double array[3] = {(double) ni, (double) nj, l};

  return int_gsl_integrate_low_precision(int_for_C_ss_tomo_TATT_EE_limber, (void*) array, 
    fmax(amin_source(ni), amin_source(nj)), amax_source(ni), NULL, GSL_WORKSPACE_SIZE);
}

double C_ss_tomo_TATT_BB_limber_nointerp(double l, int ni, int nj)
{
  if(ni < -1 || ni > tomo.shear_Nbin -1 || nj < -1 || nj > tomo.shear_Nbin -1)
  {
    log_fatal("invalid bin input (ni, nj) = (%d, %d)", ni, nj);
    exit(1);
  }

  double array[3] = {(double) ni, (double) nj, l};

  return int_gsl_integrate_low_precision(int_for_C_ss_tomo_TATT_BB_limber, (void*) array,
    //fmax(amin_source(ni), amin_source(nj)), fmin(amax_source_IA(ni), amax_source_IA(nj)),
    fmax(amin_source(ni), amin_source(nj)), fmin(amax_source(ni), amax_source(nj)),
    NULL, GSL_WORKSPACE_SIZE);
}

double C_ss_tomo_TATT_EE_limber(double l, int ni, int nj)
{
  if(ni < -1 || ni > tomo.shear_Nbin -1 || nj < -1 || nj > tomo.shear_Nbin -1)
  {
    log_fatal("invalid bin input (ni, nj) = (%d, %d)", ni, nj);
    exit(1);
  }

  static cosmopara C;
  static nuisancepara N;
  static double** table;
  static double* sig;
  static int osc[100];

  const int NSIZE = tomo.shear_Npowerspectra;
  const int nell = Ntable.N_ell_TATT;
  const double lnlmin = log(fmax(limits.LMIN_tab - 1., 1.0));
  const double lnlmax = log(limits.LMAX + 1);
  const double dlnl = (lnlmax - lnlmin) / (nell - 1.);

  if (table == 0)
  {
    table = (double**) malloc(sizeof(double*)*NSIZE);
    for(int i=0; i<NSIZE; i++) 
    {
      table[i] = (double*) malloc(sizeof(double)*nell);
    }
    sig = (double*) malloc(sizeof(double)*NSIZE);
  }

  if (recompute_shear(C, N))
  {
    {
      const int k = 0;
      const int Z1NZ = Z1(k);
      const int Z2NZ = Z2(k);
      sig[k] = 1.;
      osc[k] = 0;
      if (C_ss_tomo_TATT_EE_limber_nointerp(500., Z1NZ, Z2NZ) < 0)
      {
        sig[k] = -1.;
      }
      #pragma omp parallel for
      for (int i=0; i<nell; i++)
      {
        const double llog = lnlmin + i*dlnl;
        table[k][i] = C_ss_tomo_TATT_EE_limber_nointerp(exp(llog), Z1NZ, Z2NZ);
      }
      for (int i=0; i<nell; i++)
      {
        if (table[k][i] * sig[k] < 0.)
        {
          osc[k] = 1;
        }
      }
      if (osc[k] == 0)
      {
        #pragma omp parallel for
        for (int i=0; i<nell; i++)
        {
          table[k][i] = log(sig[k] * table[k][i]);
        }
      }
    }
    #pragma omp parallel for collapse(2)
    for (int k=1; k<NSIZE; k++)
    {
      for (int i=0; i<nell; i++)
      {
        const int Z1NZ = Z1(k);
        const int Z2NZ = Z2(k);
        const double llog = lnlmin + i*dlnl;
        table[k][i] = C_ss_tomo_TATT_EE_limber_nointerp(exp(llog), Z1NZ, Z2NZ);
      }
    }
    #pragma omp parallel for
    for (int k=1; k<NSIZE; k++)
    {
      const int Z1NZ = Z1(k);
      const int Z2NZ = Z2(k);
      sig[k] = 1.;
      osc[k] = 0;
      if (C_ss_tomo_TATT_EE_limber_nointerp(500., Z1NZ, Z2NZ) < 0)
      {
        sig[k] = -1.;
      }
      for (int i=0; i<nell; i++)
      {
        if (table[k][i] * sig[k] < 0.)
        {
          osc[k] = 1;
        }
      }
      if (osc[k] == 0)
      {
        #pragma omp parallel for
        for (int i = 0; i<nell; i++)
        {
          table[k][i] = log(sig[k] * table[k][i]);
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

  int k = N_shear(ni, nj);
  if(k > NSIZE - 1)
  {
    log_fatal("error in selecting bin number");
    exit(1);
  }

  double f1;
  if (osc[k] == 0)
  {
    f1 = sig[k] * exp(interpol(table[k], nell, lnlmin, lnlmax, dlnl, lnl, 0., 0.));
  }
  else
  {
    f1 = interpol(table[k], nell, lnlmin, lnlmax, dlnl, lnl, 0., 0.);
  }
  if (isnan(f1))
  {
    f1 = 0.;
  }
  return f1;
}

double C_ss_tomo_TATT_BB_limber(double l, int ni, int nj)
{
  if(ni < -1 || ni > tomo.shear_Nbin -1 || nj < -1 || nj > tomo.shear_Nbin -1)
  {
    log_fatal("invalid bin input (ni, nj) = (%d, %d)", ni, nj);
    exit(1);
  }

  static cosmopara C;
  static nuisancepara N;
  static double** table;
  static double* sig;
  static int osc[100];

  const int nell = Ntable.N_ell_TATT;
  const int NSIZE = tomo.shear_Npowerspectra;
  const double lnlmin = log(fmax(limits.LMIN_tab - 1., 1.0));
  const double lnlmax = log(limits.LMAX + 1);
  const double dlnl = (lnlmax - lnlmin) / (nell - 1.);

  if (table == 0)
  {
    table = (double**) malloc(sizeof(double*)*NSIZE);
    for(int i=0; i<NSIZE; i++) 
    {
      table[i] = (double*) malloc(sizeof(double)*nell);
    }
    sig = (double*) malloc(sizeof(double)*NSIZE);
  }
  if (recompute_shear(C, N))
  {
    {
      const int k = 0;
      const int Z1NZ = Z1(k);
      const int Z2NZ = Z2(k);
      sig[k] = 1.;
      osc[k] = 0;
      if (C_ss_tomo_TATT_BB_limber_nointerp(500., Z1NZ, Z2NZ) < 0)
      {
        sig[k] = -1.;
      }
      #pragma omp parallel for
      for (int i=0; i<nell; i++)
      {
        const double llog = lnlmin + i*dlnl;
        table[k][i] = C_ss_tomo_TATT_BB_limber_nointerp(exp(llog), Z1NZ, Z2NZ);
      }
      for (int i=0; i<nell; i++)
      {
        if (table[k][i] * sig[k] < 0.)
        {
          osc[k] = 1;
        }
      }
      if (osc[k] == 0)
      {
        #pragma omp parallel for
        for (int i=0; i<nell; i++)
        {
          table[k][i] = log(sig[k] * table[k][i]);
        }
      }
    }
    #pragma omp parallel for collapse(2)
    for (int k=1; k<NSIZE; k++)
    {
      for (int i=0; i<nell; i++)
      {
        const int Z1NZ = Z1(k);
        const int Z2NZ = Z2(k);
        const double llog = lnlmin + i*dlnl;
        table[k][i] = C_ss_tomo_TATT_BB_limber_nointerp(exp(llog), Z1NZ, Z2NZ);
      }
    }
    #pragma omp parallel for
    for (int k=1; k<NSIZE; k++)
    {
      const int Z1NZ = Z1(k);
      const int Z2NZ = Z2(k);
      sig[k] = 1.;
      osc[k] = 0;
      if (C_ss_tomo_TATT_BB_limber_nointerp(500., Z1NZ, Z2NZ) < 0)
      {
        sig[k] = -1.;
      }
      for (int i=0; i<nell; i++)
      {
        if (table[k][i] * sig[k] < 0.)
        {
          osc[k] = 1;
        }
      }
      if (osc[k] == 0)
      {
        #pragma omp parallel for
        for (int i=0; i<nell; i++)
        {
          table[k][i] = log(sig[k] * table[k][i]);
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

  const int k = N_shear(ni, nj);
  if(k > NSIZE - 1)
  {
    log_fatal("error in selecting bin number");
    exit(1);
  }

  double f1;
  if (osc[k] == 0)
  {
    f1 = sig[k] * exp(interpol(table[k], nell, lnlmin, lnlmax, dlnl, lnl, 0, 0));
  }
  else
  {
    f1 = interpol(table[k], nell, lnlmin, lnlmax, dlnl, lnl, 0, 0);
  }
  if (isnan(f1)) 
  {
    f1 = 0.;
  }
  return f1;
}

// -----------------------------------------------------------------------------
// SS ANGULAR CORRELATION FUNCTION - NON TATT
// -----------------------------------------------------------------------------

double int_for_C_ss_tomo_limber(double a, void* params)
{
  if(!(a>0) || !(a<1)) 
  {
    log_fatal("a>0 and a<1 not true");
    exit(1);
  }
  double* ar = (double* ) params;

  // Cocoa: added extra options to reduce code duplication
  const int use_linear_ps = ar[3];

  // prefactor correction (1812.05995 eqs 74-79)
  const double ell_prefactor = (ar[2] - 1.)*(ar[2])*(ar[2] + 1.)*(ar[2] + 2.);

  const double growfac_a = growfac(a);
  struct chis chidchi = chi_all(a);
  const double hoverh0 = hoverh0v2(a, chidchi.dchida);

  const double ell = ar[2] + 0.5;
  const double fK = f_K(chidchi.chi);
  const double k = ell/fK;
  const double ell4 = ell*ell*ell*ell;

  const double ws1 = W_source(a, ar[0], hoverh0);
  const double ws2 = W_source(a, ar[1], hoverh0);
  const double wk1 = W_kappa(a, fK, ar[0]);
  const double wk2 = W_kappa(a, fK, ar[1]);

  double res = 0.0;

  switch(like.IA)
  {
    case 0:
    {
      res = wk1*wk2;

      break;
    }
    case 1:
    {
      const double norm =
        A_IA_Joachimi(a)*cosmology.Omega_m*nuisance.c1rhocrit_ia/growfac_a;

      res = ws1*ws2*norm*norm - (ws1*wk2+ws2*wk1)*norm + wk1*wk2;

      break;
    }
    case 3:
    {
      const double norm = cosmology.Omega_m*nuisance.c1rhocrit_ia/growfac_a;

      res = ws1*ws2*norm*norm - (ws1*wk2+ws2*wk1)*norm + wk1*wk2;

      break;
    }
    case 4:
    {
      const double norm = cosmology.Omega_m*nuisance.c1rhocrit_ia/growfac_a*
        nuisance.A_ia*pow(1./(a*nuisance.oneplusz0_ia), nuisance.eta_ia);

      res = ws1*ws2*norm*norm - (ws1*wk2+ws2*wk1)*norm + wk1*wk2;

      break;
    }
    default:
      log_fatal("like.IA = %d not supported", like.IA);
      exit(1);
  }

  if(use_linear_ps == 1)
  {
    res *= p_lin(k,a);
  }
  else
  {
    res *= Pdelta(k,a);
  }

  return res*(chidchi.dchida/(fK*fK))*ell_prefactor/ell4;
}

double C_ss_tomo_limber_nointerp(double l, int ni, int nj, int use_linear_ps)
{
  if(ni < -1 || ni > tomo.shear_Nbin -1 || nj < -1 || nj > tomo.shear_Nbin -1)
  {
    log_fatal("invalid bin input (ni, nj) = (%d, %d)", ni, nj);
    exit(1);
  }

  double array[4] = {(double) ni, (double) nj, l, (double) use_linear_ps};

  int j,k;
  if (ni <= nj)
  {
    j = nj;
    k = ni;
  }
  else
  {
    j = ni;
    k = nj;
  }

  switch(like.IA)
  { // different IA might require different integrator precision
    case 0:
    {
      return int_gsl_integrate_low_precision(int_for_C_ss_tomo_limber, (void*) array,
        amin_source(j), 0.99999, NULL, GSL_WORKSPACE_SIZE);
      break;
    }
    case 1:
    {
      return int_gsl_integrate_low_precision(int_for_C_ss_tomo_limber, (void*) array,
        amin_source(j), amax_source(k), NULL, GSL_WORKSPACE_SIZE);
      break;
    }
    case 3:
    {
      return int_gsl_integrate_low_precision(int_for_C_ss_tomo_limber,
        (void*) array, amin_source(j), amax_source(k), NULL, GSL_WORKSPACE_SIZE);
      break;
    }
    case 4:
    {
      return int_gsl_integrate_low_precision(int_for_C_ss_tomo_limber, (void*) array,
        amin_source(j), amax_source(k), NULL, GSL_WORKSPACE_SIZE);
      break;
    }
    default:
      log_fatal("like.IA = %d not supported", like.IA);
      exit(1);
  }
}

double C_ss_tomo_limber(double l, int ni, int nj)
{
  if(ni < -1 || ni > tomo.shear_Nbin -1 || nj < -1 || nj > tomo.shear_Nbin -1)
  {
    log_fatal("invalid bin input (ni, nj) = (%d, %d)", ni, nj);
    exit(1);
  }

  static cosmopara C;
  static nuisancepara N;
  static double **table;  

  const int nell = Ntable.N_ell;
  const int NSIZE = tomo.shear_Npowerspectra;
  const double lnlmin = log(fmax(limits.LMIN_tab - 1., 1.0));
  const double lnlmax = log(limits.LMAX + 1);
  const double dlnl = (lnlmax - lnlmin)/(nell - 1.);

  if (table == 0)
  {
    table = (double**) malloc(sizeof(double*)*NSIZE);
    for(int i=0; i<NSIZE; i++) 
    {
      table[i] = (double*) malloc(sizeof(double)*nell);
    }
  }
  if (recompute_shear(C, N))
  {
    {
      const int k = 0;
      const int Z1NZ = Z1(k);
      const int Z2NZ = Z2(k);
      {
        const int i = 0;
        const double lnl = lnlmin + i*dlnl;
        const double l = exp(lnl);
        table[k][i]= log(C_ss_tomo_limber_nointerp(l, Z1NZ, Z2NZ, use_linear_ps_limber));
      }
      #pragma omp parallel for
      for (int i=1; i<nell; i++)
      {
        const double lnl = lnlmin + i*dlnl;
        const double l = exp(lnl);
        table[k][i]= log(C_ss_tomo_limber_nointerp(l, Z1NZ, Z2NZ, use_linear_ps_limber));
      }
    }
    #pragma omp parallel for collapse(2)
    for (int k=1; k<NSIZE; k++)
    {
      for (int i=0; i<nell; i++)
      {
        const int Z1NZ = Z1(k);
        const int Z2NZ = Z2(k); 
        const double lnl = lnlmin + i*dlnl;
        const double l = exp(lnl);
        table[k][i]= log(C_ss_tomo_limber_nointerp(l, Z1NZ, Z2NZ, use_linear_ps_limber));
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

  const int q = N_shear(ni, nj);
  if(q > NSIZE - 1)
  {
    log_fatal("error in selecting bin number");
    exit(1);
  }
  
  const double f1 = exp(interpol_fitslope(table[q], nell, lnlmin, lnlmax, dlnl, lnl, 1.));
  return isnan(f1) ? 0.0 : f1;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

// -----------------------------------------------------------------------------
// GS ANGULAR CORRELATION FUNCTION - TATT
// -----------------------------------------------------------------------------

double int_for_C_gs_tomo_limber_TATT(double a, void* params)
{
  if(!(a>0) || !(a<1)) 
  {
    log_fatal("a>0 and a<1 not true");
    exit(1);
  }

  double* ar = (double*) params;
  const double l = ar[2];
  const int ni = (int) ar[0];
  const int nj = (int) ar[1];

  const double growfac_a = growfac(a);
  struct chis chidchi = chi_all(a);
  const double hoverh0 = hoverh0v2(a, chidchi.dchida);

  const double ell = l + 0.5;
  const double fK = f_K(chidchi.chi);
  const double k = ell / fK;

  // radial n_z weight for source bin (for use with IA term)
  const double ws = W_source(a, nj, hoverh0);
  // radial lens efficiency for source bin
  const double wk = W_kappa(a, fK, nj);
  // IA parameters for first source bin
  const double C1 = C1_TA(a, nj, growfac_a);
  const double b_ta = b_TA(a, nj);
  const double C2 = C2_TT(a, nj, growfac_a);
  // radial n_z weight for lens bin (for use with clustering term)
  const double w_density = W_HOD(a, ni, hoverh0);
  // lens efficiency *b_mag for lens bin (for lens magnification)
  const double w_mag = W_mag(a, fK, ni) * gbias.b_mag[ni];
  // galaxy bias parameters for lens bin
  const double b1 = gbias.b1_function(1. / a - 1., ni);
  const double b2 = gbias.b2[ni];
  const double bs2 = gbias.bs2[ni];

  const double g4 = growfac_a*growfac_a*growfac_a*growfac_a;
  const double Pnl = Pdelta(k, a);
  double P_1loop = b1 * Pnl;
  if (w_density * b2 != 0)
  {
    P_1loop += g4 * (0.5 * b2 * PT_d1d2(k) + 0.5 * bs2 * PT_d1s2(k) +
                     0.5 * b3nl_from_b1(b1) * PT_d1d3(k));
  }

  // 1-loop P_gm ggl terms
  double res = w_density * wk * P_1loop;
  // lens magnification x G term
  res += w_mag * wk * Pnl;
  // (linear bias lens density + lens magnification) with TATT_GI terms

  //COCOA: Took these evaluations of the parenthesis - to force them to update
  //COCOA: the static variables in the first call that is done outside OpenMP loop
  const double tmp1 = TATT_GI_E(k, a, C1, C2, b_ta, growfac_a, Pnl);
  if (C1 || C2)
  {
    res += (b1 * w_density + w_mag) * ws * tmp1;
  }
  return res * chidchi.dchida / fK / fK;
}

// -----------------------------------------------------------------------------
// GS ANGULAR CORRELATION FUNCTION - NON TATT
// -----------------------------------------------------------------------------

double int_for_C_gs_tomo_limber(double a, void* params)
{
  if(!(a>0) || !(a<1)) 
  {
    log_fatal("a>0 and a<1 not true");
    exit(1);
  }
  double* ar = (double*) params;

  const int use_linear_ps = (int) ar[3];

  // prefactor correction (1812.05995 eqs 74-79)
  const double ell_prefactor1 = (ar[2])*(ar[2]+1.);
  double ell_prefactor2 = (ar[2]-1.)*ell_prefactor1*(ar[2]+2.);
  if (ell_prefactor2 <= 0.)
  {
    ell_prefactor2 = 0.;
  }
  else
  {
    ell_prefactor2 = sqrt(ell_prefactor2);
  }

  const double growfac_a = growfac(a);
  struct chis chidchi = chi_all(a);
  const double hoverh0 = hoverh0v2(a, chidchi.dchida);

  const double ell = ar[2]+0.5;
  const double fK = f_K(chidchi.chi);
  const double k = ell/fK;

  const double wgal = W_gal(a, ar[0], chidchi.chi, hoverh0) +
    W_mag(a, fK, ar[0])*(ell_prefactor1/ell/ell -1.)*gbias.b_mag[(int) ar[0]];

  double res = 0.0;
  switch(like.IA)
  {
    case 0:
    {
      res = W_kappa(a, fK, ar[1]);

      break;
    }
    case 1:
    {
      const double norm =
        A_IA_Joachimi(a)*cosmology.Omega_m*nuisance.c1rhocrit_ia/growfac_a;

      res = (W_kappa(a, fK, ar[1]) - W_source(a, ar[1], hoverh0)*norm);

      break;
    }
    case 3:
    {
      const double norm = nuisance.A_z[(int)ar[1]]*cosmology.Omega_m*
        nuisance.c1rhocrit_ia/growfac_a;

      res = (W_kappa(a, fK, ar[1]) - W_source(a, ar[1], hoverh0)*norm);

      break;
    }
    case 4:
    {
      const double norm = cosmology.Omega_m*nuisance.c1rhocrit_ia/growfac_a*
        nuisance.A_ia*pow(1./(a*nuisance.oneplusz0_ia),nuisance.eta_ia);

      res = (W_kappa(a, fK, ar[1]) - W_source(a, ar[1], hoverh0)*norm);

      break;
    }
    default:
      log_fatal("like.IA = %d not supported", like.IA);
      exit(1);
  }

  const double PK = use_linear_ps == 1 ? p_lin(k,a) : Pdelta(k,a);

  if(include_RSD_GS == 1)
  {
    const double chi_0 = f_K(ell/k);
    const double chi_1 = f_K((ell+1.)/k);
    const double a_0 = a_chi(chi_0);
    const double a_1 = a_chi(chi_1);

    res *= (wgal + W_RSD(ell, a_0, a_1, ar[0]))*;
  }
  else
  {
    res *= wgal;
  }
  return res*PK*(chidchi.dchida/(fK*fK))*(ell_prefactor2/(ell*ell));
}

double int_for_C_gs_tomo_limber_withb2(double a, void* params)
{
  if(!(a>0) || !(a<1)) 
  {
    log_fatal("a>0 and a<1 not true");
    exit(1);
  }
  double* ar = (double*) params;

  // prefactor correction (1812.05995 eqs 74-79)
  const double ell_prefactor1 = (ar[2])*(ar[2]+1.);
  double ell_prefactor2 = (ar[2]-1.)*ell_prefactor1*(ar[2]+2.);
  if (ell_prefactor2 <= 0.)
  {
    ell_prefactor2 = 0.;
  }
  else
  {
    ell_prefactor2 = sqrt(ell_prefactor2);
  }

  const double b1 = gbias.b1_function(1./a - 1., (int)ar[0]);
  const double b2 = gbias.b2[(int)ar[0]];
  const double bs2 = gbias.bs2[(int)ar[0]];

  const double growfac_a = growfac(a);
  struct chis chidchi = chi_all(a);
  const double hoverh0 = hoverh0v2(a, chidchi.dchida);
  const double g4 = growfac_a*growfac_a*growfac_a*growfac_a;

  const double ell = ar[2] + 0.5;
  const double fK = f_K(chidchi.chi);
  const double k = ell/fK;

  const double wgal = W_gal(a, ar[0], chidchi.chi, hoverh0) +
    W_mag(a, fK, ar[0])*(ell_prefactor1/ell/ell -1.)*gbias.b_mag[(int) ar[0]];

  const double PK = Pdelta(k,a);
  const double ws = W_source(a, ar[1], hoverh0);
  const double wk = W_kappa(a, fK, ar[1]);

  double linear_part = 0.0;
  double non_linear_part = 0.0;

  switch (like.IA)
  {
    case 0:
    {
      linear_part = wk;
      non_linear_part = wk;

      break;
    }
    case 1:
    {
      const double norm =
        A_IA_Joachimi(a)*cosmology.Omega_m*nuisance.c1rhocrit_ia/growfac_a;

      linear_part = (wk - ws*norm);
      non_linear_part = (wk - ws*norm);

      break;
    }
    case 3:
    {
      const double norm = nuisance.A_z[(int)ar[1]]*cosmology.Omega_m*
        nuisance.c1rhocrit_ia/growfac_a;

      linear_part = (wk - ws*norm);
      non_linear_part = (wk - ws*norm);

      break;
    }
    case 4:
    {
      const double norm = cosmology.Omega_m*nuisance.c1rhocrit_ia/growfac_a*
        nuisance.A_ia*pow(1./(a*nuisance.oneplusz0_ia),nuisance.eta_ia);

      linear_part = (wk - ws*norm);
      non_linear_part = (wk - ws*norm);

      break;
    }
    default:
      log_fatal("like.IA = %d not supported", like.IA);
      exit(1);
  }

  non_linear_part *= W_HOD(a, ar[0], hoverh0);
  non_linear_part *= g4*(0.5*b2*PT_d1d2(k) +
      0.5*bs2*PT_d1s2(k) + 0.5*b3nl_from_b1(b1)*PT_d1d3(k));

  if(include_RSD_GS == 1)
  {
    const double chi_0 = f_K(ell/k);
    const double chi_1 = f_K((ell+1.)/k);
    const double a_0 = a_chi(chi_0);
    const double a_1 = a_chi(chi_1);

    linear_part *= (wgal + W_RSD(ell, a_0, a_1, ar[0]))*PK;
  }
  else
  {
    linear_part *= wgal*PK;
  }

  return (linear_part + non_linear_part)*
    (chidchi.dchida/(fK*fK))*(ell_prefactor2/(ell*ell));
}

double C_gs_tomo_limber_nointerp(double l, int nl, int ns, int use_linear_ps, 
int init_static_vars_only)
{
  if(nl < -1 || nl > tomo.clustering_Nbin -1 || ns < -1 || ns > tomo.shear_Nbin -1)
  {
    log_fatal("invalid bin input (ni, nj) = (%d, %d)", nl, ns);
    exit(1);
  }

  double array[4] = {(double) nl, (double) ns, l, use_linear_ps};
  double result;

  if(like.IA == 0 || like.IA == 1 || like.IA == 3 || like.IA == 4)
  {
    if (gbias.b2[nl] && use_linear_ps == 0)
    {
      if (init_static_vars_only == 1)
      {
        result = int_for_C_gs_tomo_limber_withb2(amin_lens(nl), (void*) array);
      }
      else
      {
        result = int_gsl_integrate_low_precision(int_for_C_gs_tomo_limber_withb2, (void*) array,
          amin_lens(nl), 0.99999, NULL, GSL_WORKSPACE_SIZE);
      }
    }
    if (init_static_vars_only == 1)
    {
      result = int_for_C_gs_tomo_limber(amin_lens(nl), (void*) array);
    }
    else
    {
      result = int_gsl_integrate_medium_precision(int_for_C_gs_tomo_limber, (void*) array,
        amin_lens(nl), 0.99999, NULL, GSL_WORKSPACE_SIZE);
    }
  }
  else if (like.IA == 5 || like.IA == 6)
  {
    if(use_linear_ps == 1)
    {
      log_fatal("use linear power spectrum option not implemented with TATT");
      exit(1);
    }
    if (init_static_vars_only == 1)
    {
      result = int_for_C_gs_tomo_limber_TATT(amin_lens(nl), (void*) array);
    }
    else
    {
      result = int_gsl_integrate_low_precision(int_for_C_gs_tomo_limber_TATT,
        (void*) array, amin_lens(nl), 0.9999, NULL, GSL_WORKSPACE_SIZE);
    }
  }
  else
  {
    log_fatal("like.IA = %d not supported", like.IA);
    exit(1);
  }
  return result;
}

double C_gs_tomo_limber(double l, int ni, int nj)
{
  if(ni < -1 || ni > tomo.clustering_Nbin -1 || nj < -1 || nj > tomo.shear_Nbin -1)
  {
    log_fatal("invalid bin input (ni, nj) = (%d, %d)", ni, nj);
    exit(1);
  }

  static cosmopara C;
  static nuisancepara N;
  static galpara G;
  static double** table;
  static double* sig;
  static int osc[100];

  const int NSIZE = tomo.ggl_Npowerspectra;
  const int nell = (like.IA == 5 || like.IA == 6) ? Ntable.N_ell_TATT :  Ntable.N_ell;
  const double lnlmin = log(fmax(limits.LMIN_tab - 1., 1.0));
  const double lnlmax = log(limits.LMAX + 1);
  const double dlnl = (lnlmax - lnlmin) / (nell - 1.);

  if (table == 0) 
  {
    table = (double**) malloc(sizeof(double*)*NSIZE);
    for(int i=0; i<NSIZE; i++) 
    {
      table[i] = (double*) malloc(sizeof(double)*nell);
    }
    sig = (double*) malloc(sizeof(double)*NSIZE);
  }
  if (recompute_gs(C, G, N))
  {
    if (like.IA == 5 || like.IA == 6) // TATT MODELING
    { 
      C_gs_tomo_limber_nointerp(lnlmin, ZL(0), ZS(0), use_linear_ps_limber, 1); // init static vars
      #pragma omp parallel for collapse(2)
      for (int k=0; k<NSIZE; k++)
      {
        for (int i=0; i<nell; i++)
        {
          const int ZLNZ = ZL(k);
          const int ZSNZ = ZS(k);
          const double lnl = lnlmin + i*dlnl;
          const double l = exp(lnl);
          table[k][i] = C_gs_tomo_limber_nointerp(l, ZLNZ, ZSNZ, use_linear_ps_limber, 0);
        }
      }
      #pragma omp parallel for
      for (int k=0; k<NSIZE; k++)
      {
        const int ZLNZ = ZL(k);
        const int ZSNZ = ZS(k);
        sig[k] = 1.;
        osc[k] = 0;
        if (C_gs_tomo_limber_nointerp(500., ZLNZ, ZSNZ, use_linear_ps_limber, 0) < 0)
        {
          sig[k] = -1.;
        }
        for (int i=0; i<nell; i++)
        {
          if (table[k][i] * sig[k] < 0.)
          {
            osc[k] = 1;
          }
        }
        if (osc[k] == 0)
        {
          for (int i=0; i<nell; i++)
          {
            table[k][i] = log(sig[k] * table[k][i]);
          }
        }
      }
    }
    else
    {
      C_gs_tomo_limber_nointerp(lnlmin, ZL(0), ZS(0), use_linear_ps_limber, 1); // init static vars
      #pragma omp parallel for collapse(2)
      for (int k=0; k<NSIZE; k++)
      {
        for (int i=0; i<nell; i++)
        {
          const int ZLNZ = ZL(k);
          const int ZSNZ = ZS(k);
          const double lnl = lnlmin + i*dlnl;
          const double l = exp(lnl);
          table[k][i] = log(C_gs_tomo_limber_nointerp(l, ZLNZ, ZSNZ, use_linear_ps_limber, 0));
        }
      }
    }
    update_cosmopara(&C);
    update_nuisance(&N);
    update_galpara(&G);
  }
  
  const double lnl = log(l);
  if (lnl < lnlmin || lnl > lnlmax)
  {
    log_fatal("l = %e outside look-up table range [%e,%e]", l, exp(lnlmin), exp(lnlmax));
    exit(1);
  }

  const int k = N_ggl(ni, nj);
  if(k > NSIZE - 1)
  {
    log_fatal("error in selecting bin number");
    exit(1);
  }

  if (like.IA == 5 || like.IA == 6) // TATT MODELING
  {     
    double f1 = 0.;
    if (test_zoverlap(ni, nj) && osc[k] == 0)
    {
      f1 = sig[k] * exp(interpol(table[k], nell, lnlmin, lnlmax, dlnl, lnl, 0, 0));
    }
    else if (test_zoverlap(ni, nj) && osc[k] == 1)
    {
      f1 = interpol(table[k], nell, lnlmin, lnlmax, dlnl, lnl, 0, 0);
    }
    return isnan(f1) ? 0.0 : f1;
  }
  else
  {
    if(test_zoverlap(ni, nj))
    {
      const double f1 = exp(interpol_fitslope(table[k], nell, lnlmin, lnlmax, dlnl, lnl, 1));
      return isnan(f1) ? 0.0 : f1;
    }
    else
    {
      return 0.0;
    }
  }
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

double int_for_C_gg_tomo_limber(double a, void* params)
{
  if(!(a>0) || !(a<1)) 
  {
    log_fatal("a>0 and a<1 not true");
    exit(1);
  }
  double* ar = (double*) params;

  const int use_linear_ps = (int) ar[3];

  // prefactor correction (1812.05995 eqs 74-79)
  const double ell_prefactor = (ar[1])*(ar[1] + 1.);

  struct chis chidchi = chi_all(a);
  const double hoverh0 = hoverh0v2(a, chidchi.dchida);

  const double ell = ar[2] + 0.5;
  const double fK = f_K(chidchi.chi);
  const double k = ell / fK;

  const double wgal1 = W_gal(a, ar[0], chidchi.chi, hoverh0) +
    W_mag(a, fK, ar[0])*(ell_prefactor/(ell*ell) - 1.0)*gbias.b_mag[(int) ar[0]];

  const double wgal2 = W_gal(a, ar[1], chidchi.chi, hoverh0) +
    W_mag(a, fK, ar[1])*(ell_prefactor/(ell*ell) - 1.0)*gbias.b_mag[(int) ar[1]];

  const double PK = use_linear_ps == 1 ? p_lin : Pdelta;
  
  double tmp = 0.0;
  if(include_RSD_GG == 1)
  {
    const double chi_0 = f_K(ell/k);
    const double chi_1 = f_K((ell + 1.)/k);
    const double a_0 = a_chi(chi_0);
    const double a_1 = a_chi(chi_1);
    tmp = (wgal1 + W_RSD(ell, a_0, a_1, ar[0]))*(wgal2 + W_RSD(ell, a_0, a_1, ar[1]));
  }
  else
  {
    tmp = wgal1*wgal2;
  }
  return (tmp*PK*chidchi.dchida/(fK*fK));
}

double int_for_C_gg_tomo_limber_withb2(double a, void* params)
{
  if(!(a>0) || !(a<1)) 
  {
    log_fatal("a>0 and a<1 not true");
    exit(1);
  }
  double* ar = (double*) params;

  // prefactor correction (1812.05995 eqs 74-79)
  const double ell_prefactor = (ar[1])*(ar[1] + 1.);

  const double b1 = gbias.b1_function(1./a - 1., (int) ar[0]);
  const double b2 = gbias.b2[(int) ar[0]];
  const double bs2 = gbias.bs2[(int) ar[0]];

  const double growfac_a = growfac(a);
  struct chis chidchi = chi_all(a);
  const double hoverh0 = hoverh0v2(a, chidchi.dchida);
  const double g4 = growfac_a*growfac_a*growfac_a*growfac_a;

  const double ell = ar[2] + 0.5;
  const double fK = f_K(chidchi.chi);
  const double k = ell/fK;
  const double s4 = 0.; // PT_sigma4(k);

  const double wgal1 = W_gal(a, ar[0], chidchi.chi, hoverh0) +
    W_mag(a, fK, ar[0])*(ell_prefactor/(ell*ell) - 1.0)*gbias.b_mag[(int) ar[0]];

  const double wgal2 = W_gal(a, ar[1], chidchi.chi, hoverh0) +
    W_mag(a, fK, ar[1])*(ell_prefactor/(ell*ell) - 1.0)*gbias.b_mag[(int) ar[1]];

  const double non_linear_part = g4*W_HOD(a, ar[0], hoverh0)*
    W_HOD(a, ar[1], hoverh0)*
    (b1 * b2 * PT_d1d2(k) + 0.25 * b2 * b2 * (PT_d2d2(k) - 2. * s4) +
    b1 * bs2 * PT_d1s2(k) + 0.5 * b2 * bs2 * (PT_d2s2(k) - 4. / 3. * s4) +
    0.25 * bs2 * bs2 * (PT_s2s2(k) - 8. / 9. * s4) +
    b1 * b3nl_from_b1(b1) * PT_d1d3(k));

  const double PK = Pdelta(k, a);
  
  double linear_part = 0.0;
  if(include_RSD_GG == 1)
  {
    const double chi_0 = f_K(ell/k);
    const double chi_1 = f_K((ell+1.)/k);
    const double a_0 = a_chi(chi_0);
    const double a_1 = a_chi(chi_1);

    linear_part = (wgal1 + W_RSD(ell, a_0, a_1, ar[0]))*
      (wgal2 + W_RSD(ell, a_0, a_1, ar[1]))*PK;
  }
  else
  {
    linear_part = wgal1*wgal2*PK;
  }
  return (linear_part + non_linear_part)*(chidchi.dchida/(fK*fK));
}

 // WARNING: C_gg beyond linear bias for cross-tomography bins not yet supported
double C_gg_tomo_limber_nointerp(double l, int ni, int nj, int use_linear_ps,
int init_static_vars_only)
{
  if(ni < -1 || ni > tomo.clustering_Nbin -1 || nj < -1 || nj > tomo.clustering_Nbin -1)
  {
    log_fatal("invalid bin input (ni, nj) = (%d, %d)", ni, nj);
    exit(1);
  }

  double array[4] = {(double) ni, (double) nj, l, (double) use_linear_ps};

  if ((gbias.b2[ni] || gbias.b2[nj]) && use_linear_ps == 0)
  {
    if (ni != nj)
    {
      log_fatal(
        "cross-tomography (ni, nj) = (%d, %d) bins not supported beyond linear bias", ni, nj);
      exit(1);
    }
    if (init_static_vars_only == 1)
    {
      return int_for_C_gg_tomo_limber_withb2(amin_lens(ni), (void*) array);
    }
    else
    {
      return int_gsl_integrate_low_precision(int_for_C_gg_tomo_limber_withb2, (void*) array,
        amin_lens(ni), amax_lens(ni), NULL, GSL_WORKSPACE_SIZE);
    }
  }
  else
  {
    if (init_static_vars_only == 1)
    {
      return int_for_C_gg_tomo_limber(amin_lens(nj), (void*) array);
    }
    else
    {
      return int_gsl_integrate_low_precision(int_for_C_gg_tomo_limber, (void*) array,
        amin_lens(nj), 0.99999, NULL, GSL_WORKSPACE_SIZE);
    }
  }
}

double C_gg_tomo_limber(double l, int ni, int nj) // cross redshift bin not supported 
{ 
  if(ni < -1 || ni > tomo.clustering_Nbin -1 || nj < -1 || nj > tomo.clustering_Nbin -1)
  {
    log_fatal("invalid bin input (ni, nj) = (%d, %d)", ni, nj);
    exit(1);
  }

  static cosmopara C;
  static nuisancepara N;
  static galpara G;
  static double** table;

  const int nell = Ntable.N_ell;
  const double lnlmin = log(fmax(limits.LMIN_tab - 1., 1.0));
  const double lnlmax = log(limits.LMAX+1);
  const double dlnl = (lnlmax - lnlmin) / (nell);
  const int NSIZE = tomo.clustering_Nbin;

  if (table == 0)
  {
    table = (double**) malloc(sizeof(double*)*NSIZE);
    for(int i=0; i<NSIZE; i++) 
    {
      table[i] = (double*) malloc(sizeof(double)*nell);
    }
  }
  if (recompute_gg(C, G, N))
  {
    C_gg_tomo_limber_nointerp(lnlmin, 0, 0, use_linear_ps_limber, 1) // init static vars only
    #pragma omp parallel for collapse(2)
    for (int k=0; k<NSIZE; k++)  
    {
      for(int p=0; p<nell; p++)
      {
        const int ZCL1 = k; // cross redshift bin not supported so not using ZCL1(k)
        const int ZCL2 = k; // cross redshift bin not supported so not using ZCL2(k)
        const double lnl = lnlmin + p*dlnl;
        const double l = exp(lnl);
        const double result = C_gg_tomo_limber_nointerp(l, ZCL1, ZCL2, use_linear_ps_limber, 0);
        if (result <= 0)
        {
          table[k][p] = -100;
        }
        else
        {
          table[k][p] = log(result);
        }
      }
    }
    update_galpara(&G);
    update_cosmopara(&C);
    update_nuisance(&N);
  }

  if (ni != nj) 
  {
    log_fatal("cross-tomography not supported");
    exit(1);
  }

  const int q = ni; // cross redshift bin not supported so not using N_CL(ni, nj) instead of ni
  if (q > NSIZE - 1)
  {
    log_fatal("error in selecting bin number");
    exit(1);
  }

  const double lnl = log(l);
  if (lnl < lnlmin || lnl > lnlmax)
  {
    log_fatal("l = %e outside look-up table range [%e,%e]", l, exp(lnlmin), exp(lnlmax));
    exit(1);
  }

  const double f1 = exp(interpol_fitslope(table[q], nell, lnlmin, lnlmax, dlnl, lnl, 1));
  return isnan(f1) ? 0.0 : f1;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

double int_for_C_gk_limber(double a, void* params)
{
  if(!(a>0) || !(a<1)) 
  {
    log_fatal("a>0 and a<1 not true");
    exit(1);
  }
  double* ar = (double*) params;

  const int use_linear_ps = (int) ar[2];

  // prefactor correction (1812.05995 eqs 74-79)
  const double ell_prefactor = (ar[1])*(ar[1] + 1.);

  struct chis chidchi = chi_all(a);
  const double hoverh0 = hoverh0v2(a, chidchi.dchida);

  const double ell = ar[1] + 0.5;
  const double fK = f_K(chidchi.chi);
  const double k = ell/fK;

  const double wgal = W_gal(a, ar[0], chidchi.chi, hoverh0) +
    W_mag(a, fK, ar[0])*(ell_prefactor/(ell*ell) - 1.0)*gbias.b_mag[(int) ar[0]];

  const double WK = W_k(a, fK);
  const double PK = use_linear_ps == 1 ? p_lin(k,a) : Pdelta(k,a);

  double tmp;
  if(include_RSD_GK == 1)
  {
    const double chi_0 = f_K(ell/k);
    const double chi_1 = f_K((ell+1.)/k);
    const double a_0 = a_chi(chi_0);
    const double a_1 = a_chi(chi_1);
    tmp = (wgal + W_RSD(ell, a_0, a_1, ar[0]))*WK;
  }
  else
  {
    tmp = wgal*WK;
  }
  return (res*PK*chidchi.dchida/(fK*fK))*(ell_prefactor/(ell*ell));
}

double int_for_C_gk_limber_withb2(double a, void* params)
{
  if(!(a>0) || !(a<1)) 
  {
    log_fatal("a>0 and a<1 not true");
    exit(1);
  }
  double* ar = (double*) params;

  // prefactor correction (1812.05995 eqs 74-79)
  const double ell_prefactor = (ar[1])*(ar[1] + 1.);

  //const double b1 = gbias.b1_function(1./a - 1.,(int) ar[0]);
  const double b2 = gbias.b2[(int) ar[0]];
  const double bs2 = gbias.bs2[(int) ar[0]];

  const double growfac_a = growfac(a);
  const double g4 = growfac_a*growfac_a*growfac_a*growfac_a;

  struct chis chidchi = chi_all(a);
  const double hoverh0 = hoverh0v2(a, chidchi.dchida);
  const double ell = ar[1] + 0.5;
  const double fK = f_K(chidchi.chi);
  const double k = ell/fK;
  const double WK = W_k(a, fK);

  const double wgal = W_gal(a, ar[0], chidchi.chi, hoverh0) +
    W_mag(a, fK, ar[0])*(ell_prefactor/(ell*ell) - 1.0)*gbias.b_mag[(int) ar[0]];

  const double non_linear_part = g4*W_HOD(a, ar[0], hoverh0)*WK*(
    0.5*b2*PT_d1d2(k) + 0.5*bs2*PT_d1s2(k));

  const double PK = Pdelta(k, a);

  double linear_part = 0.0;
  if(include_RSD_GK == 1)
  {
    const double chi_0 = f_K(ell/k);
    const double chi_1 = f_K((ell+1.)/k);
    const double a_0 = a_chi(chi_0);
    const double a_1 = a_chi(chi_1);

    linear_part = (wgal + W_RSD(ell, a_0, a_1, ar[0]))*WK*PK;
  }
  {
    linear_part = wgal*WK*PK;
  }

  return (linear_part + non_linear_part)*(chidchi.dchida/(fK*fK))*
    (ell_prefactor/(ell*ell));
}

double C_gk_tomo_limber_nointerp(double l, int nl, int use_linear_ps)
{
  if(nl < -1 || nl > tomo.clustering_Nbin -1)
  {
    log_fatal("invalid bin input ni = %d", nl);
    exit(1);
  }

  double array[3] = {(double) nl, l, (double) use_linear_ps};

  if ((gbias.b2[nl] || gbias.b2[nl]) && use_linear_ps == 0)
  {
    return int_gsl_integrate_medium_precision(int_for_C_gk_limber_withb2,
      (void*) array, amin_lens(nl), amax_lens(nl), NULL, GSL_WORKSPACE_SIZE);
  }
  return int_gsl_integrate_medium_precision(int_for_C_gk_limber,
    (void*) array, amin_lens(nl), 0.99999, NULL, GSL_WORKSPACE_SIZE);
}

double C_gk_tomo_limber(double l, int ni)
{
  if(ni < -1 || ni > tomo.clustering_Nbin -1)
  {
    log_fatal("invalid bin input ni = %d", ni);
    exit(1);
  }

  static cosmopara C;
  static nuisancepara N;
  static galpara G;
  static double** table;

  const int NSIZE = tomo.clustering_Nbin;
  const int nell = Ntable.N_ell;
  const double lnlmin = log(fmax(limits.LMIN_tab - 1., 1.0));
  const double lnlmax = log(limits.LMAX + 1);
  const double dlnl = (lnlmax - lnlmin)/(nell);

  if (table == 0)
  {
    table = (double**) malloc(sizeof(double*)*NSIZE);
    for(int i=0; i<NSIZE; i++) 
    {
      table[i] = (double*) malloc(sizeof(double)*nell);
    }
  }
  if (recompute_gk(C, G, N))
  {
    {
      const int k = 0;
      {
        const int i = 0;
        const double lnl = lnlmin + i*dlnl;
        const double l = exp(lnl);
        table[k][i] = log(C_gk_tomo_limber_nointerp(l, k, use_linear_ps_limber));
      }
      #pragma omp parallel for
      for (int i=1; i<nell; i++)
      {
        const double lnl = lnlmin + i*dlnl;
        const double l = exp(lnl);
        table[k][i]= log(C_gk_tomo_limber_nointerp(l, k, use_linear_ps_limber));
      }
    }
    #pragma omp parallel for collapse(2)
    for (int k=1; k<NSIZE; k++)
    {
      for (int i=0; i<nell; i++)
      {
        const double lnl = lnlmin + i*dlnl;
        const double l = exp(lnl);
        table[k][i]= log(C_gk_tomo_limber_nointerp(l, k, use_linear_ps_limber));
      }
    }
    update_cosmopara(&C);
    update_nuisance(&N);
    update_galpara(&G);
  }

  const double lnl = log(l);
  if (lnl < lnlmin || lnl > lnlmax)
  {
    log_fatal("l = %e outside look-up table range [%e,%e]", l, exp(lnlmin), exp(lnlmax));
    exit(1);
  }
 
  const int q =  ni; 
  if(q > NSIZE - 1)
  {
    log_fatal("error in selecting bin number");
    exit(1);
  } 

  const double f1 = exp(interpol(table[q], nell, lnlmin, lnlmax, dlnl, lnl, 1, 1));
  return isnan(f1) ? 0.0 : f1;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

double int_for_C_ks_limber(double a, void* params)
{
  if(!(a>0) || !(a<1)) 
  {
    log_fatal("a>0 and a<1 not true");
    exit(1);
  }
  double* ar = (double*) params;

  const int use_linear_ps = (int) ar[2];

  const double growfac_a = growfac(a);
  struct chis chidchi = chi_all(a);
  const double hoverh0 = hoverh0v2(a, chidchi.dchida);

  // prefactor correction (1812.05995 eqs 74-79)
  const double ell_prefactor1 = (ar[1])*(ar[1] + 1.);
  double ell_prefactor2 = (ar[1] - 1.)*ell_prefactor1*(ar[1] + 2.);
  if(ell_prefactor2 <= 0.)
  {
    ell_prefactor2 = 0.;
  }
  else
  {
    ell_prefactor2 = sqrt(ell_prefactor2);
  }

  const double ell = ar[1] + 0.5;
  const double fK = f_K(chidchi.chi);
  const double k = ell/fK;
  const double ell4 = ell*ell*ell*ell;

  const double ws1 = W_source(a, ar[0], hoverh0);
  const double wk1 = W_kappa(a, fK, ar[0]);
  const double wk2 = W_k(a, fK);
  
  const double PK = use_linear_ps == 1 ? p_lin(k,a) : Pdelta(k,a);

  double res = 0;
  switch(like.IA)
  {
    case 0:
    {
      res = wk1*wk2;

      break;
    }
    case 1:
    {
      const double norm = A_IA_Joachimi(a)*cosmology.Omega_m*nuisance.c1rhocrit_ia/growfac_a;

      res = (-ws1*wk2*norm + wk1*wk2);

      break;
    }
    case 3:
    {
      const double norm = cosmology.Omega_m*nuisance.c1rhocrit_ia/growfac_a;

      res = (-ws1*wk2*norm + wk1*wk2);

      break;
    }
    case 4:
    {
      const double norm = (cosmology.Omega_m*nuisance.c1rhocrit_ia/growfac_a)*
        nuisance.A_ia*pow(1.0/(a*nuisance.oneplusz0_ia), nuisance.eta_ia);

      res = (-ws1*wk2*norm + wk1*wk2);

      break;
    }
    default:
    {
      log_fatal("like.IA = %d not supported", like.IA);
      exit(1);
    }
  }
  return (res*PK*chidchi.dchida/(fK*fK))*ell_prefactor1*ell_prefactor2/ell4;
}

double C_ks_tomo_limber_nointerp(double l, int nj, int use_linear_ps)
{
  if(nj < -1 || nj > tomo.shear_Nbin -1)
  {
    log_fatal("invalid bin input nj = %d", nj);
    exit(1);
  }

  double array[3] = {(double) nj, l, (double) use_linear_ps};

  switch(like.IA)
  { // different IA might require different integrator precision
    case 0:
    {
      return int_gsl_integrate_medium_precision(int_for_C_ks_limber, (void*) array,
        amin_source(nj), 0.99999, NULL, GSL_WORKSPACE_SIZE);

      break;
    }
    case 1:
    {
      return int_gsl_integrate_medium_precision(int_for_C_ks_limber,
        (void*) array, amin_source(nj), amax_source(nj), NULL, GSL_WORKSPACE_SIZE);

      break;
    }
    case 3:
    { 
      return int_gsl_integrate_medium_precision(int_for_C_ks_limber, (void*) array,
        amin_source(nj), amax_source(nj), NULL, GSL_WORKSPACE_SIZE);

      break;
    }
    case 4:
    {
      return int_gsl_integrate_medium_precision(int_for_C_ks_limber, (void*) array,
        amin_source(nj), amax_source(nj), NULL, GSL_WORKSPACE_SIZE);

      break;
    }
    default:
      log_fatal("like.IA = %d not supported", like.IA);
      exit(1);
  }
}

double C_ks_tomo_limber(double l, int ni)
{
  if(ni < -1 || ni > tomo.shear_Nbin -1)
  {
    log_fatal("invalid bin input ni = %d", ni);
    exit(1);
  }

  static cosmopara C;
  static nuisancepara N;
  static double** table;
  static double* sig;
  static int osc[100];

  const int NSIZE = tomo.shear_Nbin;
  const int nell = Ntable.N_ell;
  const double lnlmin = log(fmax(limits.LMIN_tab - 1., 1.0));
  const double lnlmax = log(limits.LMAX + 1);
  const double dlnl = (lnlmax - lnlmin)/(nell);

  if (table == 0)
  {
    table = (double**) malloc(sizeof(double*)*NSIZE);
    for(int i=0; i<NSIZE; i++) 
    {
      table[i] = (double*) malloc(sizeof(double)*nell);
    }
    sig = (double*) malloc(sizeof(double)*NSIZE);
  }
  if (recompute_shear(C, N))
  {
    {
      const int k = 0;
      sig[k] = 1.;
      osc[k] = 0;
      if (C_ks_tomo_limber_nointerp(500., k, use_linear_ps_limber) < 0)
      {
        sig[k] = -1.;
      }
      #pragma omp parallel for
      for (int i=0; i<nell; i++)
      {
        const double lnl = lnlmin + i*dlnl;
        table[k][i] = C_ks_tomo_limber_nointerp(exp(lnl), k, use_linear_ps_limber);
      }
      for (int i=0; i<nell; i++)
      {
        if (table[k][i]*sig[k] < 0.)
        {
          osc[k] = 1;
        }
      }
      if (osc[k] == 0)
      {
        #pragma omp parallel for
        for(int i=0; i<nell; i++)
        {
          table[k][i] = log(sig[k]*table[k][i]);
        }
      }
    }
    #pragma omp parallel for collapse(2)
    for (int k=1; k<NSIZE; k++)
    {
      for (int i=0; i<nell; i++)
      {
        const double lnl = lnlmin + i*dlnl;
        table[k][i] = C_ks_tomo_limber_nointerp(exp(lnl), k, use_linear_ps_limber);
      }
    }
    
    for (int k=1; k<NSIZE; k++)
    {
      sig[k] = 1.;
      osc[k] = 0;
      if (C_ks_tomo_limber_nointerp(500., k, use_linear_ps_limber) < 0)
      {
        sig[k] = -1.;
      }
      for (int i=0; i<nell; i++)
      {
        if (table[k][i]*sig[k] < 0.)
        {
          osc[k] = 1;
        }
      }
      if (osc[k] == 0)
      {
        #pragma omp parallel for
        for(int i=0; i<nell; i++)
        {
          table[k][i] = log(sig[k]*table[k][i]);
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

  const int q =  ni; 
  if(q > NSIZE - 1)
  {
    log_fatal("error in selecting bin number");
    exit(1);
  } 

  double f1 = 0.;
  if (osc[ni] == 0)
  {
    f1 = sig[ni]*exp(interpol_fitslope(table[q], nell, lnlmin, lnlmax, dlnl, lnl, 1));
  }
  if (osc[ni] == 1)
  {
    f1 = interpol_fitslope(table[q], nell, lnlmin, lnlmax, dlnl, lnl, 1);
  }
  if (isnan(f1))
  {
    f1 = 0.;
  }
  return f1;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

double int_for_C_kk_limber(double a, void* params)
{
  if(!(a>0) || !(a<1)) 
  {
    log_fatal("a>0 and a<1 not true");
    exit(1);
  }
  double* ar = (double*) params;

  const int use_linear_ps = (int) ar[1];

  struct chis chidchi = chi_all(a);

  // prefactor correction (1812.05995 eqs 74-79)
  const double ell_prefactor = (ar[1])*(ar[1] + 1.);

  const double ell = ar[0] + 0.5;
  const double fK = f_K(chidchi.chi);
  const double k = ell/fK;
  const double ell4 = ell*ell*ell*ell;
  const double WK = W_k(a, fK);

  const double PK = use_linear_ps == 1 ? p_lin(k,a) : Pdelta(k,a);
  return WK*WK*PK*(chidchi.dchida/(fK*fK))*ell_prefactor*ell_prefactor/ell4;
}

double C_kk_limber_nointerp(double l, int use_linear_ps)
{
  double ar[2] = {l, (double) use_linear_ps};
  return int_gsl_integrate_medium_precision(int_for_C_kk_limber, (void*) ar,
    limits.a_min*(1.+1.e-5), 1.-1.e-5, NULL, GSL_WORKSPACE_SIZE);
}

double C_kk_limber(double l)
{
  static cosmopara C;
  static double* table;
  
  const int nell = Ntable.N_ell;
  const double lnlmin = log(fmax(limits.LMIN_tab - 1., 1.0));
  const double lnlmax = log(limits.LMAX + 1);
  const double dlnl = (lnlmax - lnlmin)/(nell);

  if (table == 0)
  {
    table = (double*) malloc(sizeof(double)*nell);
  }
  if (recompute_cosmo3D(C))
  {
    {
      const int i = 0;
      const double lnl = lnlmin + i*dlnl;
       table[i] = log(C_kk_limber_nointerp(exp(lnl), use_linear_ps_limber));
    }
    #pragma omp parallel for
    for (int i=1; i<nell; i++)
    {
      const double lnl = lnlmin + i*dlnl;
       table[i] = log(C_kk_limber_nointerp(exp(lnl), use_linear_ps_limber));
    }
    update_cosmopara(&C);
  }

  const double lnl = log(l);
  if (lnl < lnlmin || lnl > lnlmax)
  {
    log_fatal("l = %e outside look-up table range [%e,%e]", l, exp(lnlmin), exp(lnlmax));
    exit(1);
  }
  
  double f1 = exp(interpol(table, nell, lnlmin, lnlmax, dlnl, lnl, 1., 1.));
  return isnan(f1) ? 0.0 : f1;
}

// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// Non Limber (Angular Power Spectrum)
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------

// --------------------------------------------------------------------------------------
// Galaxy Clustering
// --------------------------------------------------------------------------------------

void f_chi_for_Psi_cl(double* chi, int Nchi, double* f_chi, int ni, double zmin, double zmax)
{ // Integrand for galaxy density
  if(ni < -1 || ni > tomo.clustering_Nbin -1)
  {
    log_fatal("invalid bin input ni = %d", ni);
    exit(1);
  }
  const double real_coverH0 = cosmology.coverH0/cosmology.h0; // unit Mpc
  {
    const int i = 0;
    const double a = a_chi(chi[i]/real_coverH0 /* convert unit to c/H0 */);
    const double z = 1. / a - 1.;
    const double tmp1 = pf_photoz(z, ni);
    const double pf = (tmp1 < 0.) ? 0 : tmp1; // get rid of unphysical negatives
    f_chi[i] = gbias.b1_function(z, ni)*chi[i]*pf*growfac(a)*hoverh0(a)/real_coverH0;
    if ((z < zmin) || (z > zmax))
    {
      f_chi[i] = 0.;
    }
  }
  #pragma omp parallel for
  for (int i=1; i<Nchi; i++)
  {
    const double a = a_chi(chi[i]/real_coverH0 /* convert unit to c/H0 */);
    const double z = 1.0/a - 1.;
    const double tmp1 = pf_photoz(z, ni);
    const double pf = (tmp1 < 0.) ? 0 : tmp1; // get rid of unphysical negatives
    f_chi[i] = gbias.b1_function(z, ni)*chi[i]*pf*growfac(a)*hoverh0(a)/real_coverH0;
    if ((z < zmin) || (z > zmax))
    {
      f_chi[i] = 0.;
    }
  }
}

void f_chi_for_Psi_cl_RSD(double* chi, int Nchi, double* f_chi, int ni, double zmin, double zmax)
{ // Integrand for galaxy density RSD
  if(ni < -1 || ni > tomo.clustering_Nbin -1)
  {
    log_fatal("invalid bin input ni = %d", ni);
    exit(1);
  }
  const double real_coverH0 = cosmology.coverH0 / cosmology.h0;
  {
    const int i = 0;
    const double a = a_chi(chi[i]/real_coverH0 /* convert unit to c/H0 */);
    const double z = 1.0 / a - 1.0;
    const double tmp1 = pf_photoz(z, ni);
    const double pf = (tmp1 < 0.) ? 0 : tmp1; // get rid of unphysical negatives
    struct growths tmp2 = growfac_all(a);
    f_chi[i] = -chi[i]*pf*tmp2.D*tmp2.f*hoverh0(a)/real_coverH0;
    if ((z < zmin) || (z > zmax))
    {
      f_chi[i] = 0.; 
    } 
  }
  #pragma omp parallel for
  for (int i=1; i<Nchi; i++) 
  {
    const double a = a_chi(chi[i]/real_coverH0 /* convert unit to c/H0 */);
    const double z = 1.0 / a - 1.0;
    const double tmp1 = pf_photoz(z, ni);
    const double pf = (tmp1 < 0.) ? 0 : tmp1; // get rid of unphysical negatives
    struct growths tmp2 = growfac_all(a);
    f_chi[i] = -chi[i]*pf*tmp2.D*tmp2.f*hoverh0(a)/real_coverH0;
    if ((z < zmin) || (z > zmax))
    {
      f_chi[i] = 0.;
    } 
  }
}

void f_chi_for_Psi_cl_Mag(double* chi, int Nchi, double* f_chi, int ni, double zmax)
{ // Integrand for lensing magnification of galaxy density
  if(ni < -1 || ni > tomo.clustering_Nbin -1)
  {
    log_fatal("invalid bin input ni = %d", ni);
    exit(1);
  }
  const double real_coverH0 = cosmology.coverH0 / cosmology.h0;
  {
    const int i = 0;
    const double a = a_chi(chi[i]/real_coverH0 /* convert unit to c/H0 */);
    const double z = 1. / a - 1.;
    const double fK = f_K(chi[i]/real_coverH0 /* convert unit to c/H0 */);
    const double wmag = W_mag(a, fK, ni);
    const double window_M = wmag/fK/(real_coverH0*real_coverH0);
    f_chi[i] = window_M * growfac(a); // unit [Mpc^-2]
    if (z > zmax)
    {
      f_chi[i] = 0.;
    } 
  }
  #pragma omp parallel for
  for (int i=1; i<Nchi; i++)
  {
    const double a = a_chi(chi[i]/real_coverH0 /* convert unit to c/H0 */);
    const double z = 1. / a - 1.;
    const double fK = f_K(chi[i]/real_coverH0 /* convert unit to c/H0 */);
    const double wmag = W_mag(a, fK, ni);
    const double window_M = wmag/fK/(real_coverH0*real_coverH0);
    f_chi[i] = window_M * growfac(a); // unit [Mpc^-2]
    if (z > zmax)
    {
      f_chi[i] = 0.;
    } 
  }
}

void C_cl_tomo(int L, int ni, int nj, double* Cl, double dev, double tol) 
{
  if(ni < -1 || ni > tomo.clustering_Nbin -1 || nj < -1 || nj > tomo.clustering_Nbin -1)
  {
    log_fatal("invalid bin input (ni, nj) = (%d, %d)", ni, nj);
    exit(1);
  }
  if (ni != nj)
  {
    log_fatal("Cocoa disabled cross-spectrum w_gg");
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

  const double real_coverH0 = cosmology.coverH0/cosmology.h0;
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
    for (int i = 0; i < Nell_block; i++) 
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
  for (int i=0; i < Nchi; i++)
  { // chi_min and chi_max are cosmology dependent
    chi_ar[i] = chi_min * exp(dlnchi * i); 
  }
  #pragma omp parallel for collapse(2)
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
  
  {
    const double zmin = tomo.clustering_zmin[ni];
    const double zmax = tomo.clustering_zmax[ni];
    f_chi_for_Psi_cl(chi_ar, Nchi, f1_chi, ni, zmin, zmax);
    f_chi_for_Psi_cl_RSD(chi_ar, Nchi, f1_chi_RSD, ni, zmin, zmax);
    f_chi_for_Psi_cl_Mag(chi_ar, Nchi, f1_chi_Mag, ni, zmax);
  }
  if (ni != nj)
  { 
    const double zmin = tomo.clustering_zmin[nj];
    const double zmax = tomo.clustering_zmax[nj];
    f_chi_for_Psi_cl(chi_ar, Nchi, f2_chi, nj, zmin, zmax);
    f_chi_for_Psi_cl_RSD(chi_ar, Nchi, f2_chi_RSD, nj, zmin, zmax);
    f_chi_for_Psi_cl_Mag(chi_ar, Nchi, f2_chi_Mag, nj, zmax);
  }

  config cfg;
  cfg.nu = 1.;
  cfg.c_window_width = 0.25;
  cfg.derivative = 0;
  cfg.N_pad = 200;
  cfg.N_extrap_low = 0;
  cfg.N_extrap_high = 0;

  config cfg_RSD;
  cfg_RSD.nu = 1.01;
  cfg_RSD.c_window_width = 0.25;
  cfg_RSD.derivative = 2;
  cfg_RSD.N_pad = 500;
  cfg_RSD.N_extrap_low = 0;
  cfg_RSD.N_extrap_high = 0;

  config cfg_Mag;
  cfg_Mag.nu = 1.;
  cfg_Mag.c_window_width = 0.25;
  cfg_Mag.derivative = 0;
  cfg_Mag.N_pad = 500;
  cfg_Mag.N_extrap_low = 0;
  cfg_Mag.N_extrap_high = 0;

  int i_block = 0;

  while ((fabs(dev) > tol) && (L < limits.LMAX_NOLIMBER))
  {
    for (int i=0; i<Nell_block; i++)
    {
      ell_ar[i] = i + i_block * Nell_block; 
    }

    cfftlog_ells(chi_ar, f1_chi, Nchi, &cfg, ell_ar, Nell_block, k1, Fk1);
    cfftlog_ells_increment(chi_ar, f1_chi_RSD, Nchi, &cfg_RSD, ell_ar, Nell_block, k1, Fk1);   
    cfftlog_ells(chi_ar, f1_chi_Mag, Nchi, &cfg_Mag, ell_ar, Nell_block, k1, Fk1_Mag);    
    if (ni != nj)
    {
      cfftlog_ells(chi_ar, f2_chi, Nchi, &cfg, ell_ar, Nell_block, k2, Fk2);
      cfftlog_ells_increment(chi_ar, f2_chi_RSD, Nchi, &cfg_RSD, ell_ar, Nell_block, k2, Fk2);
      cfftlog_ells(chi_ar, f2_chi_Mag, Nchi, &cfg_Mag, ell_ar, Nell_block, k2, Fk2_Mag);
    }

    #pragma omp parallel for collapse(2)
    for (int i=0; i<Nell_block; i++)
    {
      for (int j=0; j<Nchi; j++)
      {
        const double ell_prefactor = ell_ar[i] * (ell_ar[i] + 1.);
        Fk1[i][j] += gbias.b_mag[ni]*ell_prefactor*Fk1_Mag[i][j]/(k1[i][j]*k1[i][j]);
        if (ni != nj)
        {
          Fk2[i][j] += gbias.b_mag[nj]*ell_prefactor*Fk2_Mag[i][j]/(k2[i][j]*k2[i][j]);
        }
      }
    }
    
    {
      const int i = 0;
      double cl_temp = 0.;
      for (int j=0; j<Nchi; j++)
      {
        const double k1cH0 = k1[i][j] * real_coverH0;
        const double PK = p_lin(k1cH0, 1.0);
        const double k1cH03 = k1cH0*k1cH0*k1cH0;
        cl_temp += (ni == nj) ? Fk1[i][j]*Fk1[i][j]*k1cH03*PK : Fk1[i][j]*Fk2[i][j]*k1cH03*PK;
      }
      Cl[ell_ar[i]] = cl_temp * dlnk * 2.0 / M_PI +
         C_gg_tomo_limber_nointerp((double) ell_ar[i], ni, nj, use_linear_ps_limber, 0)
        -C_gg_tomo_limber_nointerp((double) ell_ar[i], ni, nj, 1, 0);
    }
    #pragma omp parallel for
    for (int i=1; i<Nell_block; i++)
    {
      double cl_temp = 0.;
      for (int j=0; j<Nchi; j++)
      {
        const double k1cH0 = k1[i][j] * real_coverH0;
        const double PK = p_lin(k1cH0, 1.0);
        const double k1cH03 = k1cH0*k1cH0*k1cH0;
        cl_temp += (ni == nj) ? Fk1[i][j]*Fk1[i][j]*k1cH03*PK : Fk1[i][j]*Fk2[i][j]*k1cH03*PK;
      }
      Cl[ell_ar[i]] = cl_temp * dlnk * 2. / M_PI +
         C_gg_tomo_limber_nointerp((double) ell_ar[i], ni, nj, use_linear_ps_limber, 0)
        -C_gg_tomo_limber_nointerp((double) ell_ar[i], ni, nj, 1, 0);
    }

    i_block++;

    if(L >= limits.LMAX_NOLIMBER - Nell_block)
    { //Xiao: break before memory leak in next iteration
      break;
    }
    L = i_block * Nell_block - 1;
    dev = Cl[L]/C_gg_tomo_limber_nointerp((double) L, ni, nj, use_linear_ps_limber, 0) - 1;
  }
  L++;
  
  Cl[limits.LMAX_NOLIMBER+1] = C_gg_tomo_limber((double) limits.LMAX_NOLIMBER+1, ni, nj);
  #pragma omp parallel for
  for (int l=L; l<limits.LMAX_NOLIMBER; l++)
  {
    Cl[l] = (l<limits.LMIN_tab) ? 
      C_gg_tomo_limber_nointerp((double) l, ni, nj, use_linear_ps_limber, 0) :
      C_gg_tomo_limber((double) l, ni, nj);
  }
}

// --------------------------------------------------------------------------------------
// Galaxy-Galaxylensing
// --------------------------------------------------------------------------------------

void f_chi_for_Psi_sh(double* chi, int Nchi, double* fchi, int nj, double zmax)
{
  if(nj < -1 || nj > tomo.shear_Nbin -1)
  {
    log_fatal("invalid bin input nj = %d", nj);
    exit(1);
  }

  const double real_coverH0 = cosmology.coverH0/cosmology.h0;
  {
    const int i = 0;
    const double a = a_chi(chi[i] / real_coverH0) ;
    const double z = 1./a - 1.;
    const double fK = f_K(chi[i]/real_coverH0);
    const double wkappa = W_kappa(a, fK, nj);
    fchi[i] = (wkappa/fK/(real_coverH0*real_coverH0))*growfac(a); // unit [Mpc^-2]
    if (z > zmax)
    {
      fchi[i] = 0.;
    }
  }
  #pragma omp parallel for
  for (int i=1; i<Nchi; i++)
  {
    const double a = a_chi(chi[i] / real_coverH0) ;
    const double z = 1./a - 1.;
    const double fK = f_K(chi[i]/real_coverH0);
    const double wkappa = W_kappa(a, fK, nj);
    fchi[i] = (wkappa/fK/(real_coverH0*real_coverH0))*growfac(a); // unit [Mpc^-2]
    if (z > zmax)
    {
      fchi[i] = 0.;
    }
  }
}

// TODO: ADD ALL IA POSSIBILITIES
void f_chi_for_Psi_sh_IA(double* chi, int Nchi, double* fchi, int nj, double zmin, double zmax)
{
  if(nj < -1 || nj > tomo.shear_Nbin -1)
  {
    log_fatal("invalid bin input nj = %d", nj);
    exit(1);
  }
  const double real_coverH0 = cosmology.coverH0 / cosmology.h0;
  {
    const int i = 0;
    // first convert unit of chi from Mpc to c/H0
    const double a = a_chi(chi[i] / real_coverH0) ;
    const double z = 1./a - 1.;
    struct chis chidchi = chi_all(a);
    const double hoverh0 = hoverh0v2(a, chidchi.dchida);
    const double fK = f_K(chi[i]/real_coverH0);
    const double norm = cosmology.Omega_m*nuisance.c1rhocrit_ia/growfac(a)*
      nuisance.A_ia*pow(1./(a*nuisance.oneplusz0_ia),nuisance.eta_ia);
    const double tmp1 = W_source(a, (double) nj, hoverh0);
    const double wsource = (tmp1 > 0.) ? tmp1 : 0.;
    fchi[i] = (-wsource*norm/fK/(real_coverH0*real_coverH0))*growfac(a); // unit [Mpc^-2]
    if ((z<zmin) || (z>zmax))
    {
      fchi[i] = 0.;
    }
  }
  #pragma omp parallel for
  for(int i=1; i<Nchi; i++)
  {
    // first convert unit of chi from Mpc to c/H0
    const double a = a_chi(chi[i] / real_coverH0) ;
    const double z = 1./a - 1.;
    struct chis chidchi = chi_all(a);
    const double hoverh0 = hoverh0v2(a, chidchi.dchida);
    const double fK = f_K(chi[i]/real_coverH0);
    const double norm = cosmology.Omega_m*nuisance.c1rhocrit_ia/growfac(a)*
      nuisance.A_ia*pow(1./(a*nuisance.oneplusz0_ia),nuisance.eta_ia);
    const double tmp1 = W_source(a, (double) nj, hoverh0);
    const double wsource = (tmp1 > 0.) ? tmp1 : 0.;
    const double window_ia = -wsource*norm/fK/(real_coverH0*real_coverH0);
    fchi[i] = -wsource*norm/fK/(real_coverH0*real_coverH0)*growfac(a); // unit [Mpc^-2]
    if ((z<zmin) || (z>zmax))
    {
      fchi[i] = 0.;
    }
  }
}

void C_gl_tomo(int L, int nl, int ns, double* Cl, double dev, double tolerance)
{
  if(nl < -1 || nl > tomo.clustering_Nbin -1 || ns < -1 || ns > tomo.clustering_Nbin -1)
  {
    log_fatal("invalid bin input (ni, nj) = (%d, %d)", nl, ns);
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
  double f1_chi_RSD[Nchi];
  double f1_chi_Mag[Nchi];
  double f2_chi[Nchi];
  double f2_chi_IA_ar[Nchi];

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


  {
    const double zmin = tomo.clustering_zmin[nl];
    const double zmax = tomo.clustering_zmax[nl];
    f_chi_for_Psi_cl(chi_ar, Nchi, f1_chi, nl, zmin, zmax);
    f_chi_for_Psi_cl_RSD(chi_ar, Nchi, f1_chi_RSD, nl, zmin, zmax);
    f_chi_for_Psi_cl_Mag(chi_ar, Nchi, f1_chi_Mag, nl, zmax);
  }
  {
    const double zmin = tomo.shear_zmin[ns];
    const double zmax = tomo.shear_zmax[ns];
    f_chi_for_Psi_sh(chi_ar, Nchi, f2_chi, ns, zmax);
    f_chi_for_Psi_sh_IA(chi_ar, Nchi, f2_chi_IA_ar, ns, zmin, zmax);
  }

  for(int j=0; j<Nchi; j++)
  {
    f2_chi[j] += f2_chi_IA_ar[j];
  }

  int i_block = 0;

  config cfg;
  my_config.nu = 1.;
  my_config.c_window_width = 0.25;
  my_config.derivative = 0;
  my_config.N_pad = 200;
  my_config.N_extrap_low = 0;
  my_config.N_extrap_high = 0;

  config cfg_RSD;
  my_config_RSD.nu = 1.01;
  my_config_RSD.c_window_width = 0.25;
  my_config_RSD.derivative = 2;
  my_config_RSD.N_pad = 200;
  my_config_RSD.N_extrap_low = 0;
  my_config_RSD.N_extrap_high = 0;

  config cfg_Mag;
  my_config_Mag.nu = 1.;
  my_config_Mag.c_window_width = 0.25;
  my_config_Mag.derivative = 0;
  my_config_Mag.N_pad = 1000;
  my_config_Mag.N_extrap_low = 0;
  my_config_Mag.N_extrap_high = 0;

  config cfg_shear;
  my_config_L.nu = 1.;
  my_config_L.c_window_width = 0.25;
  my_config_L.derivative = 0;
  my_config_L.N_pad = 1000.;
  my_config_L.N_extrap_low = 0;
  my_config_L.N_extrap_high = 0;

  // COCOA: LMAX_NOLIMBER might avoid infinite loop in weird models
  while ((fabs(dev) > tolerance) && (L < limits.LMAX_NOLIMBER))
  {
    for(int i=0; i<Nell_block; i++)
    {
      ell_ar[i]=i+i_block*Nell_block;
    }

    cfftlog_ells(chi_ar, f1_chi, Nchi, &cfg, ell_ar, Nell_block, k1, Fk1);
    cfftlog_ells_increment(chi_ar, f1_chi_RSD, Nchi, &cfg_RSD, ell_ar, Nell_block, k1, Fk1);
    cfftlog_ells(chi_ar, f1_chi_Mag, Nchi, &cfg_Mag, ell_ar, Nell_block, k1, Fk1_Mag);
    cfftlog_ells(chi_ar, f2_chi, Nchi, &cfg_shear, ell_ar, Nell_block, k2, Fk2);
    #pragma omp parallel for collapse(2)
    for(int i=0; i<Nell_block; i++)
    {
      for(int j=0; j<Nchi; j++)
      {
        {
          const double ell_prefactor = ell_ar[i]*(ell_ar[i]+1);
          Fk1[i][j] += gbias.b_mag[nl]*(ell_prefactor/(k1[i][j]*k1[i][j])*Fk1_Mag[i][j]);
        }
        {
          double ell_prefactor2 = (ell_ar[i]-1.)*ell_ar[i]*(ell_ar[i]+1.)*(ell_ar[i]+2.);
          if(ell_prefactor2 <= 0.)
          {
            ell_prefactor2 = 0.;
          }
          else
          {
            ell_prefactor2 = sqrt(ell_prefactor2);
          }
          Fk2[i][j] *= (ell_prefactor2/(k1[i][j]*k1[i][j]));
        }
      }
    }

    {
      const int i=0;
      double cl_temp = 0.;
      for(int j=0; j<Nchi; j++)
      {
        const double k1_cH0 = k1[i][j] * real_coverH0;
        cl_temp += Fk1[i][j]*Fk2[i][j]*k1_cH0*k1_cH0*k1_cH0*p_lin(k1_cH0, 1.0);
      }
      Cl[ell_ar[i]] = cl_temp * dlnk * 2./M_PI +
        C_gs_tomo_limber_nointerp((double) ell_ar[i], nl, ns, use_linear_ps_limber, 0)
       -C_gs_tomo_limber_nointerp((double) ell_ar[i], nl, ns, 1, 0);
    }
    #pragma omp parallel for
    for(int i=1; i<Nell_block; i++)
    {
      double cl_temp = 0.;
      for(int j=0; j<Nchi; j++)
      {
        const double k1_cH0 = k1[i][j] * real_coverH0;
        cl_temp += Fk1[i][j]*Fk2[i][j]*k1_cH0*k1_cH0*k1_cH0*p_lin(k1_cH0, 1.0);
      }
      Cl[ell_ar[i]] = cl_temp * dlnk * 2./M_PI +
        C_gs_tomo_limber_nointerp((double) ell_ar[i], nl, ns, use_linear_ps_limber, 0)
       -C_gs_tomo_limber_nointerp((double) ell_ar[i], nl, ns, 1, 0);
    }

    i_block++;

    if(L >= limits.LMAX_NOLIMBER - Nell_block)
    { //Xiao: break before memory leak in next iteration
      break;
    }

    L = i_block*Nell_block -1 ;
    dev =
    Cl[L]/C_gs_tomo_limber_nointerp((double) L, nl, ns, use_linear_ps_limber, 0)-1;
  }
  L++;

  Cl[L] = C_gg_tomo_limber((double) L, nl, ns);
  
  #pragma omp parallel for
  for (int l=L+1; l<limits.LMAX_NOLIMBER+1; l++)
  {
    Cl[l] = C_gg_tomo_limber((double) l, nl, ns);
  }
}