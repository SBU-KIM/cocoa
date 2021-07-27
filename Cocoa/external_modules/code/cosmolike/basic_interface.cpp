#include <string>
#include <vector>
#include <numeric>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <cmath>
#include <stdexcept>
#include <array>
#include <random>

#include <spdlog/spdlog.h>
#include <spdlog/sinks/basic_file_sink.h>
#include <spdlog/cfg/env.h>

#include <boost/algorithm/string.hpp>

#include "cosmolike/basics.h"
#include "cosmolike/bias.h"
#include "cosmolike/baryons.h"
#include "cosmolike/cosmo2D.h"
#include "cosmolike/cosmo3D.h"
#include "cosmolike/halo.h"
#include "cosmolike/radial_weights.h"
#include "cosmolike/recompute.h"
#include "cosmolike/pt_cfastpt.h"
#include "cosmolike/redshift_spline.h"
#include "cosmolike/structs.h"

#include "basic_interface.hpp"


namespace basic_interface
{

typedef std::vector<double> dvec;
typedef arma::Col<double>  armavec;
typedef arma::Mat<double> armamat;

// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// INIT FUNCTIONS
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------

void cpp_initial_setup()
{
  spdlog::cfg::load_env_levels();
  spdlog::debug("\x1b[90m{}\x1b[0m: Begins", "initial_setup");

  // restart variables to 0 so error check can flag bad initialization
  tomo.shear_Nbin = 0;
  tomo.clustering_Nbin = 0;
  tomo.cluster_Nbin = 0;
  tomo.shear_Npowerspectra = 0;
  tomo.cg_clustering_Npowerspectra = 0;
  tomo.cc_clustering_Npowerspectra = 0;
  tomo.ggl_Npowerspectra = 0;
  tomo.clustering_Npowerspectra = 0;
  tomo.cgl_Npowerspectra = 0;
  tomo.external_selection_cg_clustering = NULL;
  
  like.shear_shear = 0;
  like.shear_pos = 0;
  like.pos_pos = 0;
  like.clusterN = 0;
  like.clusterWL = 0;
  like.clusterCG = 0;
  like.clusterCC = 0;

  // reset bias
  for (int i=0; i<MAX_SIZE_ARRAYS; i++)
  {
    gbias.b[i] = 0.0;
    gbias.b2[i] = 0.0;
    gbias.b_mag[i] = 0.0;
  }

  // reset IA
  for (int i=0; i<MAX_SIZE_ARRAYS; i++)
  {
    nuisance.A_z[i] = 0.0;
    nuisance.A2_z[i] = 0.0;
    nuisance.b_ta_z[i] = 0.0;
  }

  spdlog::debug("\x1b[90m{}\x1b[0m: Ends", "initial_setup");
}

void cpp_init_survey(std::string surveyname, double area, double sigma_e)
{
  spdlog::debug("\x1b[90m{}\x1b[0m: Begins", "init_survey");

  if (surveyname.size() > CHAR_MAX_SIZE-1)
  {
    spdlog::critical( "\x1b[90m{}\x1b[0m: insufficient pre-allocated char memory (max = {}) for"
      "the string: {}", "init_survey", CHAR_MAX_SIZE-1, surveyname);
    exit(1);
  }
  if (!(surveyname.size()>0))
  {
    spdlog::critical("{}: incompatible input", "init_survey");
    exit(1);
  }
  if(!(area > 0))
  {
    spdlog::critical("{}: invalid survey area = {} input", "init_survey", area);
    exit(1);    
  }

  memcpy(survey.name, surveyname.c_str(), surveyname.size() + 1);
  survey.area = area;
  survey.sigma_e = sigma_e;

  spdlog::debug("\x1b[90m{}\x1b[0m: Ends", "init_survey");
}

void cpp_init_probes(std::string possible_probes) 
{
  spdlog::debug("\x1b[90m{}\x1b[0m: Begins", "init_probes");

  if (possible_probes.compare("xi") == 0)
  { 
    like.shear_shear = 1;
    spdlog::debug("\x1b[90m{}\x1b[0m: {} selected", "init_probes", "xi");
  }
  else if (possible_probes.compare("wtheta") == 0) 
  {
    like.pos_pos = 1;
    spdlog::debug("\x1b[90m{}\x1b[0m: {} selected", "init_probes", "wtheta");
  } 
  else if (possible_probes.compare("gammat") == 0) 
  {
    like.shear_pos = 1;
    spdlog::debug("\x1b[90m{}\x1b[0m: {} selected", "init_probes", "gammat");
  } 
  else if (possible_probes.compare("2x2pt") == 0) 
  {
    like.shear_pos = 1;
    like.pos_pos = 1;
    spdlog::debug("\x1b[90m{}\x1b[0m: {} selected", "init_probes", "2x2pt");
  } 
  else if (possible_probes.compare("3x2pt") == 0) 
  {
    like.shear_shear = 1;
    like.shear_pos = 1;
    like.pos_pos = 1;
    spdlog::debug("\x1b[90m{}\x1b[0m: {} selected", "init_probes", "3x2pt");
  } 
  else if (possible_probes.compare("xi_ggl") == 0)
  {
    like.shear_shear = 1;
    like.shear_pos = 1;
    spdlog::debug("\x1b[90m{}\x1b[0m: {} selected", "init_probes",  "xi + gammat");
  }
  else if (possible_probes.compare("5x2pt") == 0) 
  {
    like.shear_shear = 1;
    like.shear_pos = 1;
    like.pos_pos = 1;
    like.gk = 1;
    like.ks = 1;
    spdlog::debug("\x1b[90m{}\x1b[0m: {} selected", "init_probes", "5x2pt");
  } 
  else if (possible_probes.compare("6x2pt") == 0) 
  {
    like.shear_shear = 1;
    like.shear_pos = 1;
    like.pos_pos = 1;
    like.gk = 1;
    like.ks = 1;
    like.kk = 1;
    spdlog::debug("\x1b[90m{}\x1b[0m: {} selected", "init_probes", "6x2pt");
  } 
  else if (possible_probes.compare("ss") == 0) 
  { 
    like.shear_shear = 1;
    spdlog::debug("\x1b[90m{}\x1b[0m: {} selected", "init_probes", "ss");
  } 
  else if (possible_probes.compare("gg") == 0) 
  {
    like.pos_pos = 1;
    spdlog::debug("\x1b[90m{}\x1b[0m: {} selected", "init_probes", "gg");
  } 
  else if (possible_probes.compare("gs") == 0) 
  {
    like.shear_pos = 1;
    spdlog::debug("\x1b[90m{}\x1b[0m: {} selected", "init_probes", "gs");
  }
  else if (possible_probes.compare("ss_gs") == 0) 
  {
    like.shear_shear = 1;
    like.shear_pos = 1;
    spdlog::debug("\x1b[90m{}\x1b[0m: {} selected", "init_probes", "ss + gs (2x2pt)");
  }
  else if (possible_probes.compare("wtheta_cc") == 0)
  { // cosmolike c interface
    like.clusterCC = 1;

    spdlog::debug("\x1b[90m{}\x1b[0m: {} selected", "init_probes", "wtheta_cc");
  }
  else if (possible_probes.compare("wtheta_cg") == 0)
  { // cosmolike c interface
    like.clusterCG = 1;

    spdlog::debug("\x1b[90m{}\x1b[0m: {} selected", "init_probes", "wtheta_cg");
  }
  else if (possible_probes.compare("cluster_gammat") == 0)
  { // cosmolike c interface
    like.clusterWL;

    spdlog::debug("\x1b[90m{}\x1b[0m: {} selected", "init_probes", "cluster_gammat");
  }
  else if (possible_probes.compare("N") == 0)
  { // cosmolike c interface
    like.clusterN = 1;

    spdlog::debug("\x1b[90m{}\x1b[0m: {} selected", "init_probes", "cluster n_counts");
  }
  else if (possible_probes.compare("4x2ptN") == 0)
  { // cosmolike c interface
    like.pos_pos = 1;
    like.clusterN = 1;
    like.clusterWL = 1;
    like.clusterCG = 1;
    like.clusterCC = 1;

    spdlog::debug("\x1b[90m{}\x1b[0m: {} selected", "init_probes", "4x2ptN");
  }
  else 
  {
    spdlog::critical("\x1b[90m{}\x1b[0m: {} probe not supported", "init_probes", possible_probes);
    exit(1);
  }

  spdlog::debug("\x1b[90m{}\x1b[0m: Ends", "init_probes");
}

void cpp_init_size_data_vector_real(const bool shear_shear, const bool shear_pos, 
const bool pos_pos, const bool ks, const bool gk, const bool kk, const bool clusterCC,
const bool clusterN, const bool clusterWL, const bool clusterCG)
{
  spdlog::debug("\x1b[90m{}\x1b[0m: Begins", "init_size_data_vector");

  like.Ndata = 0.0;
  if (like.Ntheta == 0)
  {
    spdlog::critical("{}: {} not set prior to this function call",
      "init_size_data_vector_real", "like.Ntheta");
    exit(1);
  }

  if(clusterCG == 1) 
  {
    if (tomo.cg_clustering_Npowerspectra == 0)
    {
      spdlog::critical("{}: {} not set prior to this function call",
        "init_size_data_vector_real", "tomo.cg_clustering_Npowerspectra");
      exit(1);
    }

    like.Ndata += like.Ntheta*tomo.cg_clustering_Npowerspectra;
  }
  if(clusterWL == 1) 
  {
    if (tomo.cg_clustering_Npowerspectra == 0)
    {
      spdlog::critical("{}: {} not set prior to this function call",
        "init_size_data_vector_real", "tomo.cg_clustering_Npowerspectra");
      exit(1);
    }

    like.Ndata += like.Ntheta*tomo.cg_clustering_Npowerspectra;
  }
  if(clusterN == 1) 
  {
    if (Cluster.N200_Nbin == 0)
    {
      spdlog::critical("{}: {} not set prior to this function call",
        "init_size_data_vector_real", "Cluster.N200_Nbin");
      exit(1);
    }

    like.Ndata += like.Ntheta*Cluster.N200_Nbin;
  }
  if(clusterCC == 1) 
  {
    if (tomo.shear_Npowerspectra == 0)
    {
      spdlog::critical("{}: {} not set prior to this function call",
        "init_size_data_vector_real", "tomo.cc_clustering_Npowerspectra");
      exit(1);
    }

    like.Ndata += like.Ntheta*tomo.cc_clustering_Npowerspectra;
  }
  if(shear_shear == 1) 
  {
    if (tomo.shear_Npowerspectra == 0)
    {
      spdlog::critical("{}: {} not set prior to this function call",
        "init_size_data_vector_real", "tomo.shear_Npowerspectra");
      exit(1);
    }

    like.Ndata += like.Ntheta*2*tomo.shear_Npowerspectra;
  }
  if(shear_pos == 1) 
  {
    if (tomo.ggl_Npowerspectra == 0)
    {
      spdlog::critical("{}: {} not set prior to this function call",
        "init_size_data_vector_real", "tomo.ggl_Npowerspectra");
      exit(1);
    }

    like.Ndata += like.Ntheta*tomo.ggl_Npowerspectra;
  }
  if(pos_pos == 1) 
  {
    if (tomo.clustering_Npowerspectra == 0)
    {
      spdlog::critical("{}: {} not set prior to this function call",
        "init_size_data_vector_real", "tomo.clustering_Npowerspectra");
      exit(1);
    }

    like.Ndata += like.Ntheta*tomo.clustering_Npowerspectra;
  }
  if(ks == 1) 
  {
    if (tomo.shear_Nbin == 0)
    {
      spdlog::critical("{}: {} not set prior to this function call",
        "init_size_data_vector_real", "tomo.shear_Nbin");
      exit(1);
    }

    like.Ndata += like.Ntheta*tomo.shear_Nbin;
  }
  if(gk  == 1) 
  {
    if (tomo.clustering_Nbin == 0)
    {
      spdlog::critical("{}: {} not set prior to this function call",
        "init_size_data_vector_real", "tomo.clustering_Nbin");
      exit(1);
    }

    like.Ndata += like.Ntheta*tomo.clustering_Nbin;
  }
  if(kk == 1) // kk always in fourier space
  {
    if (like.Ncl == 0)
    {
      spdlog::critical("{}: {} not set prior to this function call",
        "init_size_data_vector_real", "like.Ncl");
      exit(1);
    }
    like.Ndata +=  like.Ncl;
  }

  spdlog::debug("\x1b[90m{}\x1b[0m: {} = {} selected.", "init_size_data_vector_real", "Ndata", 
    like.Ndata);

  spdlog::debug("\x1b[90m{}\x1b[0m: Ends", "init_size_data_vector_real");
}

void cpp_init_size_data_vector_fourier(const bool shear_shear, const bool shear_pos, 
const bool pos_pos, const bool ks, const bool gk, const bool kk)
{
  spdlog::debug("\x1b[90m{}\x1b[0m: Begins", "init_size_data_vector");

  like.Ndata = 0.0;
  if(shear_shear == 1) 
  {
    if (like.Ntheta == 0)
    {
      spdlog::critical("{}: {} not set prior to this function call",
        "init_size_data_vector_fourier", "like.Ntheta");
      exit(1);
    }
    if (tomo.shear_Npowerspectra == 0)
    {
      spdlog::critical("{}: {} not set prior to this function call",
        "init_size_data_vector_fourier", "tomo.shear_Npowerspectra");
      exit(1);
    }

    like.Ndata += like.Ncl*tomo.shear_Npowerspectra;
  }
  if(shear_pos == 1) 
  {
    if (like.Ntheta == 0)
    {
      spdlog::critical("{}: {} not set prior to this function call",
        "init_size_data_vector_fourier", "like.Ntheta");
      exit(1);
    }
    if (tomo.ggl_Npowerspectra == 0)
    {
      spdlog::critical("{}: {} not set prior to this function call",
        "init_size_data_vector_fourier", "tomo.ggl_Npowerspectra");
      exit(1);
    }

    like.Ndata += like.Ncl*tomo.ggl_Npowerspectra;
  }
  if(pos_pos == 1) 
  {
    if (like.Ntheta == 0)
    {
      spdlog::critical("{}: {} not set prior to this function call",
        "init_size_data_vector_fourier", "like.Ntheta");
      exit(1);
    }
    if (tomo.clustering_Npowerspectra == 0)
    {
      spdlog::critical("{}: {} not set prior to this function call",
        "init_size_data_vector_fourier", "tomo.clustering_Npowerspectra");
      exit(1);
    }

    like.Ndata += like.Ncl*tomo.clustering_Npowerspectra;
  }
  if(ks == 1) 
  {
    if (like.Ntheta == 0)
    {
      spdlog::critical("{}: {} not set prior to this function call",
        "init_size_data_vector_fourier", "like.Ntheta");
      exit(1);
    }
    if (tomo.shear_Nbin == 0)
    {
      spdlog::critical("{}: {} not set prior to this function call",
        "init_size_data_vector_fourier", "tomo.shear_Nbin");
      exit(1);
    }

    like.Ndata += like.Ncl*tomo.shear_Nbin;
  }
  if(gk  == 1) 
  {
    if (like.Ntheta == 0)
    {
      spdlog::critical("{}: {} not set prior to this function call",
        "init_size_data_vector_fourier", "like.Ntheta");
      exit(1);
    }
    if (tomo.clustering_Nbin == 0)
    {
      spdlog::critical("{}: {} not set prior to this function call",
        "init_size_data_vector_fourier", "tomo.clustering_Nbin");
      exit(1);
    }

    like.Ndata += like.Ncl*tomo.clustering_Nbin;
  }
  if(kk == 1) // kk always in fourier space
  {
    if (like.Ncl == 0)
    {
      spdlog::critical("{}: {} not set prior to this function call",
        "init_size_data_vector_fourier", "like.Ncl");
      exit(1);
    }
    like.Ndata += like.Ncl;
  }

  spdlog::debug("\x1b[90m{}\x1b[0m: {} = {} selected.", "init_size_data_vector_fourier", "Ndata", 
    like.Ndata);
  spdlog::debug("\x1b[90m{}\x1b[0m: Ends", "init_size_data_vector_fourier");
}

void cpp_init_IA(const int N)
{
  spdlog::debug("\x1b[90m{}\x1b[0m: Begins", "init_IA");

  spdlog::debug("\x1b[90m{}\x1b[0m: {} = {} selected.", "init_IA", "IA", N);

  if (N == 3 || N == 4 || N == 5 || N == 6)
  {
    like.IA = N;
  }
  else
  {
    spdlog::critical("{}: IA choice {} = {} not supported", "init_IA", "like.IA", N);
    exit(1);
  }

  spdlog::debug("\x1b[90m{}\x1b[0m: Ends", "init_IA");
}

void cpp_init_cosmo_runmode(const bool is_linear)
{
  spdlog::debug("\x1b[90m{}\x1b[0m: Begins", "init_cosmo_runmode");

  std::string mode = is_linear ? "linear" : "Halofit";
  const size_t size = mode.size();
  memcpy(pdeltaparams.runmode, mode.c_str(), size + 1);

  spdlog::debug("\x1b[90m{}\x1b[0m: {} = {} selected", "init_cosmo_runmode", "runmode", mode);

  spdlog::debug("\x1b[90m{}\x1b[0m: Ends", "init_cosmo_runmode");
}

void cpp_init_binning_real(const int Ntheta, const double thetamin_arcmin, 
const double thetamax_arcmin)
{
  spdlog::debug("\x1b[90m{}\x1b[0m: Begins", "init_binning");

  if (!(Ntheta > 0))
  {
    spdlog::critical("\x1b[90m{}\x1b[0m: {} = {} not supported", "init_binning",
      "like.Ntheta", Ntheta);
    exit(1);
  }

  spdlog::debug("\x1b[90m{}\x1b[0m: {} = {} selected.", "init_binning", "Ntheta", Ntheta);
  spdlog::debug("\x1b[90m{}\x1b[0m: {} = {} selected.", "init_binning", "min", thetamin_arcmin);
  spdlog::debug("\x1b[90m{}\x1b[0m: {} = {} selected.", "init_binning", "max", thetamax_arcmin);

  like.Ntheta = Ntheta;
  like.vtmin = thetamin_arcmin * 2.90888208665721580e-4;
  like.vtmax = thetamax_arcmin * 2.90888208665721580e-4;
  const double logdt = (std::log(like.vtmax) - std::log(like.vtmin))/like.Ntheta;
  like.theta = (double*) calloc(like.Ntheta, sizeof(double));

  constexpr double x = 2./ 3.;

  for (int i=0; i<like.Ntheta; i++)
  {
    const double thetamin = std::exp(log(like.vtmin) + (i + 0.0) * logdt);
    const double thetamax = std::exp(log(like.vtmin) + (i + 1.0) * logdt);
    like.theta[i] = x * (thetamax*thetamax*thetamax - thetamin*thetamin*thetamin)/
      (thetamax*thetamax - thetamin*thetamin);

    spdlog::debug("\x1b[90m{}\x1b[0m: Bin {:d} - {} = {:.4e}, {} = {:.4e} and {} = {:.4e}",
      "init_binning", i, "theta_min [rad]", thetamin, "theta [rad]", like.theta[i], 
      "theta_max [rad]", thetamax);
  }

  spdlog::debug("\x1b[90m{}\x1b[0m: Ends", "init_binning");
}

void cpp_init_binning_fourier(const int Ncl, const double lmin, const double lmax) 
{
  spdlog::debug("\x1b[90m{}\x1b[0m: Begins", "init_binning");
  if (!(Ncl > 0)) 
  {
    spdlog::critical("\x1b[90m{}\x1b[0m: {} = {} not supported",
      "init_binning", "like.Ncl", Ncl);
    exit(1);
  }

  spdlog::debug("\x1b[90m{}\x1b[0m: {} = {} selected.", "init_binning", "Ncl", Ncl);
  spdlog::debug("\x1b[90m{}\x1b[0m: {} = {} selected.", "init_binning", "l_min", lmin);
  spdlog::debug("\x1b[90m{}\x1b[0m: {} = {} selected.", "init_binning", "l_max", lmax);

  like.Ncl = Ncl;
  like.lmin = lmin;
  like.lmax = lmax;
  const double logdl = (std::log(like.lmax) - std::log(like.lmin))/like.Ncl;
  like.ell = (double*) calloc(Ncl, sizeof(double)); 
  
  for (int i = 0; i < like.Ncl; i++) 
  {
    like.ell[i] = std::exp(std::log(like.lmin)+(i+0.5)*logdl);
    spdlog::debug("\x1b[90m{}\x1b[0m: Bin {:d} - {} = {:.4e}, {} = {:.4e} and {} = {:.4e}",
     "init_binning", i, "lmin", lmin, "ell", like.ell[i], "lmax", lmax);
  }

  spdlog::debug("\x1b[90m{}\x1b[0m: Ends", "init_binning");
}

void cpp_init_lens_sample(std::string multihisto_file, const int Ntomo, const double ggl_cut)
{
  spdlog::debug("\x1b[90m{}\x1b[0m: Begins", "init_lens_sample");

  if (tomo.shear_Nbin == 0)
  {
    spdlog::critical("{}: {} not set prior to this function call", "init_lens_sample", 
      "tomo.shear_Nbin");
    exit(1);
  }
  if (multihisto_file.size()>CHAR_MAX_SIZE-1)
  {
    spdlog::critical( "\x1b[90m{}\x1b[0m: insufficient pre-allocated char memory (max = {}) for"
      "the string: {}", "init_lens_sample", CHAR_MAX_SIZE-1, multihisto_file);
    exit(1);
  }
  if (!(multihisto_file.size() > 0))
  {
    spdlog::critical("\x1b[90m{}\x1b[0m: empty {} string not supported",
      "init_lens_sample", "multihisto_file");
    exit(1);
  }
  if (!(Ntomo > 0) || Ntomo > MAX_SIZE_ARRAYS)
  {
    spdlog::critical("\x1b[90m{}\x1b[0m: {} = {} not supported (max = {})",
      "init_lens_sample", "Ntomo", Ntomo, MAX_SIZE_ARRAYS);
    exit(1);
  }

  memcpy(redshift.clustering_REDSHIFT_FILE, multihisto_file.c_str(), multihisto_file.size()+1);

  redshift.clustering_photoz = 4;
  tomo.clustering_Nbin = Ntomo;
  tomo.clustering_Npowerspectra = tomo.clustering_Nbin;

  spdlog::debug("\x1b[90m{}\x1b[0m: {} = {} selected.", "init_lens_sample",
    "clustering_REDSHIFT_FILE", multihisto_file);

  spdlog::debug("\x1b[90m{}\x1b[0m: {} = {} selected.", "init_lens_sample",
    "clustering_Nbin", Ntomo);

  if (ggl_cut > 0)
  {
    survey.ggl_overlap_cut = ggl_cut;
  }
  else
  {
    survey.ggl_overlap_cut = 0.0;
  }

  spdlog::debug("\x1b[90m{}\x1b[0m: {} = {} selected.", "init_lens_sample",
    "survey.ggl_overlap_cut", survey.ggl_overlap_cut);

  pf_photoz(0.1, 0);
  {
    int n = 0;
    for (int i = 0; i < tomo.clustering_Nbin; i++)
    {
      for (int j = 0; j < tomo.shear_Nbin; j++)
      {
        n += test_zoverlap(i, j);
      }
    }
    tomo.ggl_Npowerspectra = n;

    spdlog::debug("\x1b[90m{}\x1b[0m: tomo.ggl_Npowerspectra = {}",
      "init_lens_sample", tomo.ggl_Npowerspectra);
  }
  spdlog::debug("\x1b[90m{}\x1b[0m: Ends", "init_lens_sample");
}

void cpp_init_source_sample(std::string multihisto_file, const int Ntomo)
{
  spdlog::debug("\x1b[90m{}\x1b[0m: Begins", "init_source_sample");

  if (multihisto_file.size() > CHAR_MAX_SIZE - 1)
  {
    spdlog::critical(
      "\x1b[90m{}\x1b[0m: insufficient pre-allocated char memory (max = {}) for"
      "the string: {}", "init_source_sample", CHAR_MAX_SIZE-1, multihisto_file);
    exit(1);
  }
  if (!(multihisto_file.size() > 0))
  {
    spdlog::critical("\x1b[90m{}\x1b[0m: empty {} string not supported",
      "init_source_sample", "multihisto_file");
    exit(1);
  }
  if (!(Ntomo > 0) || Ntomo > MAX_SIZE_ARRAYS - 1)
  {
    spdlog::critical("\x1b[90m{}\x1b[0m: {} = {} not supported (max = {})",
      "init_source_sample", "Ntomo", Ntomo, MAX_SIZE_ARRAYS - 1);
    exit(1);
  }

  // convert std::string to char*
  memcpy(redshift.shear_REDSHIFT_FILE, multihisto_file.c_str(), multihisto_file.size() + 1);

  redshift.shear_photoz = 4;
  tomo.shear_Nbin = Ntomo;
  tomo.shear_Npowerspectra = tomo.shear_Nbin * (tomo.shear_Nbin + 1) / 2;

  spdlog::debug("\x1b[90m{}\x1b[0m: tomo.shear_Npowerspectra = {}", 
    "init_source_sample", tomo.shear_Npowerspectra);

  for (int i=0; i<tomo.shear_Nbin; i++)
  {
    nuisance.bias_zphot_shear[i] = 0.0;
    const int zmean = zmean_source(i);
    spdlog::info("\x1b[90m{}\x1b[0m: bin {} - {} = {}.", "init_source_sample", i, "<z_s>", zmean);
  }

  spdlog::debug("\x1b[90m{}\x1b[0m: {} = {} selected.", "init_source_sample",
    "shear_REDSHIFT_FILE", multihisto_file);

  spdlog::debug("\x1b[90m{}\x1b[0m: {} = {} selected.", "init_source_sample", "shear_Nbin", Ntomo);

  spdlog::debug("\x1b[90m{}\x1b[0m: Ends", "init_source_sample");
}

void cpp_init_linear_power_spectrum(dvec io_log10k, dvec io_z, dvec io_lnP)
{
  spdlog::debug("\x1b[90m{}\x1b[0m: Begins", "init_linear_power_spectrum");

  {
    bool debug_fail = false;
    if (io_z.size()*io_log10k.size() != io_lnP.size())
    {
      debug_fail = true;
    }
    else
    {
      if (io_z.size() == 0 || io_log10k.size() == 0 || io_z.size() < 5 || io_log10k.size() < 5)
      {
        debug_fail = true;
      }
    }
    if (debug_fail)
    {
      spdlog::critical("\x1b[90m{}\x1b[0m: incompatible input w/ k.size = {}, z.size = {}, and"
        "lnP.size = {}", "init_linear_power_spectrum", io_log10k.size(), io_z.size(), io_lnP.size());
      exit(1);
    }
  }

  int nlog10k = static_cast<int>(io_log10k.size());
  int nz = static_cast<int>(io_z.size());
  double* log10k = io_log10k.data();
  double* z = io_z.data();
  double* lnP = io_lnP.data();
  setup_p_lin(&nlog10k, &nz, &log10k, &z, &lnP, 1);

  // force initialization - helps avoiding seg fault when openmp is on
  const double io_a = 1.0;
  const double io_k = 0.1*cosmology.coverH0;
  p_lin(io_k, io_a);

  spdlog::debug("\x1b[90m{}\x1b[0m: Ends", "init_linear_power_spectrum");
}

void cpp_init_non_linear_power_spectrum(dvec io_log10k, dvec io_z, dvec io_lnP)
{
  spdlog::debug("\x1b[90m{}\x1b[0m: Begins", "init_non_linear_power_spectrum");

  {
    bool debug_fail = false;
    if (io_z.size()*io_log10k.size() != io_lnP.size())
    {
      debug_fail = true;
    }
    else
    {
      if (io_z.size() == 0 || io_log10k.size() == 0 || io_z.size() < 5 || io_log10k.size() < 5)
      {
        debug_fail = true;
      }
    }
    if (debug_fail)
    {
      spdlog::critical("\x1b[90m{}\x1b[0m: incompatible input w/ k.size = {}, z.size = {}, "
        "and lnP.size = {}", "init_non_linear_power_spectrum", io_log10k.size(),
        io_z.size(), io_lnP.size());
      exit(1);
    }
  }

  int nlog10k = static_cast<int>(io_log10k.size());
  int nz = static_cast<int>(io_z.size());
  double* log10k = io_log10k.data();
  double* z = io_z.data();
  double* lnP = io_lnP.data();
  setup_p_nonlin(&nlog10k, &nz, &log10k, &z, &lnP, 1);

  // force initialization - helps avoiding seg fault when openmp is on
  const double io_a = 1.0;
  const double io_k = 0.1*cosmology.coverH0;
  p_nonlin(io_k, io_a);

  spdlog::debug("\x1b[90m{}\x1b[0m: Ends", "init_non_linear_power_spectrum");
}

// Growth: D = G * a
void cpp_init_growth(dvec io_z, dvec io_G)
{
  spdlog::debug("\x1b[90m{}\x1b[0m: Begins", "init_growth");

  {
    bool debug_fail = false;
    if (io_z.size() != io_G.size())
    {
      debug_fail = true;
    }
    else
    {
      if (io_z.size() == 0 || io_z.size() < 5)
      {
        debug_fail = true;
      }
    }
    if (debug_fail)
    {
      spdlog::critical("\x1b[90m{}\x1b[0m: incompatible input w/ z.size = {} and G.size = {}",
        "init_growth", io_z.size(), io_G.size());
      exit(1);
    }
  }

  int nz = static_cast<int>(io_z.size());
  double* z = io_z.data();
  double* G = io_G.data();
  setup_growth(&nz, &z, &G, 1);

  // force initialization - helps avoiding seg fault when openmp is on
  const double io_a = 1.0;
  const double zz = 0.0;
  f_growth(zz);
  growfac_all(io_a);
  growfac(io_a);

  spdlog::debug("\x1b[90m{}\x1b[0m: Ends", "init_growth");
}

void cpp_init_distances(dvec io_z, dvec io_chi)
{
  spdlog::debug("\x1b[90m{}\x1b[0m: Begins", "init_distances");

  {
    bool debug_fail = false;
    if (io_z.size() != io_chi.size())
    {
      debug_fail = true;
    }
    else
    {
      if (io_z.size() == 0)
      {
        debug_fail = true;
      }
    }
    if (debug_fail)
    {
      spdlog::critical(
        "\x1b[90m{}\x1b[0m: incompatible input w/ z.size = {} and G.size = {}",
        "init_distances",
        io_z.size(),
        io_chi.size()
      );
      exit(1);
    }
  }

  int nz = static_cast<int>(io_z.size());
  double* vz = io_z.data();
  double* vchi = io_chi.data();
  setup_chi(&nz, &vz, &vchi, 1);

  // force initialization - imp to avoid seg fault when openmp is on
  const double io_a = 1.0;
  chi(io_a);

  spdlog::debug("\x1b[90m{}\x1b[0m: Ends", "init_distances");

  return;
}

void cpp_init_cmb(std::string cmb_noise_file) 
{
  spdlog::debug("\x1b[90m{}\x1b[0m: Begins", "init_cmb");

  if (cmb_noise_file.size() > CHAR_MAX_SIZE-1) 
  {
    spdlog::critical("\x1b[90m{}\x1b[0m: insufficient pre-allocated char memory (max = {}) for"
      "the string: {}", "init_cmb", CHAR_MAX_SIZE-1, cmb_noise_file);
    exit(1);
  }
  if (!(cmb_noise_file.size() > 0)) 
  {
    spdlog::critical("\x1b[90m{}\x1b[0m: empty {} string not supported",
      "init_cmb", "cmb_noise_file");
    exit(1);
  }

  memcpy(cmb.pathLensRecNoise, cmb_noise_file.c_str(), cmb_noise_file.size()+1);

  like.lmax_kappacmb = 2999.;

  spdlog::info("\x1b[90m{}\x1b[0m: {} = {} selected.", "init_cmb", "cmb_noise_file", cmb_noise_file);
  spdlog::debug("\x1b[90m{}\x1b[0m: Ends", "init_cmb");
}

void cpp_init_baryons_contamination(const bool use_baryonic_simulations_contamination,
const std::string which_baryonic_simulations_contamination)
{
  spdlog::debug("\x1b[90m{}\x1b[0m: Begins", "init_baryons_contamination");

  spdlog::info("\x1b[90m{}\x1b[0m: {} = {} selected", "init_baryons_contamination", 
    "use_baryonic_simulations", use_baryonic_simulations_contamination);

  if (use_baryonic_simulations_contamination)
  {
    init_baryons(which_baryonic_simulations_contamination.c_str());

    spdlog::info("\x1b[90m{}\x1b[0m: {} = {} selected",
      "init_baryons_contamination", "which_baryonic_simulations_contamination",
      which_baryonic_simulations_contamination);
  }
  else
  {
    reset_bary_struct();
  }

  spdlog::debug("\x1b[90m{}\x1b[0m: Ends", "init_baryons_contamination");
}

// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// SET PARAM FUNCTIONS
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------

void cpp_set_cosmological_parameters(const double OM, const double H0, const bool is_cached_cosmo)
{
  spdlog::debug("\x1b[90m{}\x1b[0m: Begins", "set_cosmological_parameters");

  if(!is_cached_cosmo)
  {
    // Cosmolike should not need parameters from inflation or dark energy.
    // because Cobaya provides P(k,z), H(z), D(z), Chi(z)...
    // It may require H0 to set scales and \Omega_M to set the halo model

    // cosmolike c interface
    cosmology.Omega_m = OM;
    cosmology.Omega_v = 1.0-OM;
    // Cosmolike only needs to know that there are massive neutrinos (>0)
    cosmology.Omega_nu = 0.1;
    cosmology.h0 = H0/100.0; // assuming H0 in km/s/Mpc
    cosmology.MGSigma = 0.0;
    cosmology.MGmu = 0.0;

    // Technical Problem: we want Cosmolike to calculate the data vector when
    // Cobaya request (no cache). To avoid cache in Cosmolike, we use a
    // random number generators to set cosmology.random
    cosmology.random = RandomNumber::get_instance().get();
    cosmology.is_cached = 0;
  }
  else
  {
    cosmology.is_cached = 1;
  }

  spdlog::debug("\x1b[90m{}\x1b[0m: Ends", "set_cosmological_parameters");
}

void cpp_set_nuisance_shear_calib(dvec M)
{
  spdlog::debug("\x1b[90m{}\x1b[0m: Begins", "set_nuisance_shear_calib");

  if (tomo.shear_Nbin == 0)
  {
    spdlog::critical("\x1b[90m{}\x1b[0m: {} = 0 is invalid", "set_nuisance_shear_calib",
      "shear_Nbin");
    exit(1);
  }
  if (tomo.shear_Nbin != static_cast<int>(M.size()))
  {
    spdlog::critical("\x1b[90m{}\x1b[0m: incompatible input w/ size = {} (!= {})",
      "set_nuisance_shear_calib", M.size(), tomo.shear_Nbin);
    exit(1);
  }

  for (int i=0; i<tomo.shear_Nbin; i++)
  {
    nuisance.shear_calibration_m[i] = M[i];
  }

  spdlog::debug("\x1b[90m{}\x1b[0m: Ends", "set_nuisance_shear_calib");
}

void cpp_set_nuisance_shear_photoz(dvec SP)
{
  spdlog::debug("\x1b[90m{}\x1b[0m: Begins", "set_nuisance_shear_photoz");

  if (tomo.shear_Nbin == 0)
  {
    spdlog::critical("\x1b[90m{}\x1b[0m: {} = 0 is invalid", "set_nuisance_shear_photoz",
      "shear_Nbin");
    exit(1);
  }
  if (tomo.shear_Nbin != static_cast<int>(SP.size()))
  {
    spdlog::critical("\x1b[90m{}\x1b[0m: incompatible input w/ size = {} (!= {})",
      "set_nuisance_shear_photoz", SP.size(), tomo.shear_Nbin);
    exit(1);
  }

  for (int i=0; i<tomo.shear_Nbin; i++)
  {
    nuisance.bias_zphot_shear[i] = SP.at(i);
  }

  spdlog::debug("\x1b[90m{}\x1b[0m: Ends", "set_nuisance_shear_photoz");
}

void cpp_set_nuisance_clustering_photoz(dvec CP)
{
  spdlog::debug("\x1b[90m{}\x1b[0m: Begins", "set_nuisance_clustering_photoz");

  if (tomo.clustering_Nbin == 0)
  {
    spdlog::critical("\x1b[90m{}\x1b[0m: {} = 0 is invalid", "set_nuisance_clustering_photoz",
      "clustering_Nbin");
    exit(1);
  }
  if (tomo.clustering_Nbin != static_cast<int>(CP.size()))
  {
    spdlog::critical("\x1b[90m{}\x1b[0m: incompatible input w/ size = {} (!= {})",
      "set_nuisance_clustering_photoz", CP.size(), tomo.clustering_Nbin);
    exit(1);
  }

  for (int i=0; i<tomo.clustering_Nbin; i++)
  {
    nuisance.bias_zphot_clustering[i] = CP.at(i);
  }

  spdlog::debug("\x1b[90m{}\x1b[0m: Ends", "set_nuisance_clustering_photoz");
}

void cpp_set_nuisance_linear_bias(dvec B1)
{
  spdlog::debug("\x1b[90m{}\x1b[0m: Begins", "set_nuisance_linear_bias");

  if (tomo.clustering_Nbin == 0)
  {
    spdlog::critical("\x1b[90m{}\x1b[0m: {} = 0 is invalid", "set_nuisance_linear_bias", 
      "clustering_Nbin");
    exit(1);
  }
  if (tomo.clustering_Nbin != static_cast<int>(B1.size()))
  {
    spdlog::critical("\x1b[90m{}\x1b[0m: incompatible input w/ size = {} (!= {})",
      "set_nuisance_linear_bias", B1.size(), tomo.clustering_Nbin);
    exit(1);
  }

  for (int i=0; i<tomo.clustering_Nbin; i++)
  {
    gbias.b[i] = B1.at(i);
  }

  spdlog::debug("\x1b[90m{}\x1b[0m: Ends", "set_nuisance_linear_bias");
}

void cpp_set_nuisance_nonlinear_bias(dvec B1, dvec B2)
{
  spdlog::debug("\x1b[90m{}\x1b[0m: Begins", "set_nuisance_nonlinear_bias");

  if (tomo.clustering_Nbin == 0)
  {
    spdlog::critical("\x1b[90m{}\x1b[0m: {} = 0 is invalid",
      "set_nuisance_nonlinear_bias", "clustering_Nbin"
    );
    exit(1);
  }
  if (tomo.clustering_Nbin != static_cast<int>(B1.size()) ||
      tomo.clustering_Nbin != static_cast<int>(B2.size()))
  {
    spdlog::critical("\x1b[90m{}\x1b[0m: incompatible input w/ sizes = {} and {} (!= {})",
      "set_nuisance_nonlinear_bias", B1.size(), B2.size(), tomo.clustering_Nbin
    );
    exit(1);
  }

  constexpr double tmp = -4./7.;
  for (int i=0; i<tomo.clustering_Nbin; i++)
  {
    gbias.b2[i] = B2.at(i);
    gbias.bs2[i] = almost_equal(B2.at(i), 0.) ? 0 : tmp*(B1.at(i)-1.0);
  }

  spdlog::debug("\x1b[90m{}\x1b[0m: Ends", "set_nuisance_nonlinear_bias");
}

void cpp_set_nuisance_magnification_bias(dvec B_MAG)
{
  spdlog::debug("\x1b[90m{}\x1b[0m: Begins", "set_nuisance_magnification_bias");

  if (tomo.clustering_Nbin == 0)
  {
    spdlog::critical("\x1b[90m{}\x1b[0m: {} = 0 is invalid", "set_nuisance_magnification_bias",
      "clustering_Nbin");
    exit(1);
  }
  if (tomo.clustering_Nbin != static_cast<int>(B_MAG.size()))
  {
    spdlog::critical("\x1b[90m{}\x1b[0m: incompatible input w/ size = {} (!= {})",
      "set_nuisance_magnification_bias", B_MAG.size(), tomo.clustering_Nbin);
    exit(1);
  }

  for (int i=0; i<tomo.clustering_Nbin; i++)
  {
    gbias.b_mag[i] = B_MAG.at(i);
  }

  spdlog::debug("\x1b[90m{}\x1b[0m: Ends", "set_nuisance_magnification_bias");
}

void cpp_set_nuisance_bias(std::vector<double> B1, std::vector<double> B2,
std::vector<double> B_MAG)
{
  cpp_set_nuisance_linear_bias(B1);
  cpp_set_nuisance_nonlinear_bias(B1, B2);
  cpp_set_nuisance_magnification_bias(B_MAG);
}

void cpp_set_nuisance_ia(std::vector<double> A1, std::vector<double> A2,
std::vector<double> B_TA)
{
  spdlog::debug("\x1b[90m{}\x1b[0m: Begins", "set_nuisance_ia");

  if (tomo.shear_Nbin == 0)
  {
    spdlog::critical("\x1b[90m{}\x1b[0m: {} = 0 is invalid", "set_nuisance_ia", "shear_Nbin");
    exit(1);
  }
  if (tomo.shear_Nbin != static_cast<int>(A1.size()) ||
      tomo.shear_Nbin != static_cast<int>(A2.size()) ||
  		tomo.shear_Nbin != static_cast<int>(B_TA.size()))
  {
    spdlog::critical("\x1b[90m{}\x1b[0m: incompatible input w/ sizes = {}, {} and {} (!= {})",
      "set_nuisance_ia", A1.size(), A2.size(), B_TA.size(), tomo.shear_Nbin
    );
    exit(1);
  }

  nuisance.c1rhocrit_ia = 0.01389;
  if (like.IA == 3 || like.IA == 5)
  {
    for (int i=0; i<tomo.shear_Nbin; i++)
    {
      nuisance.A_z[i] = A1[i];
      nuisance.A2_z[i] = A2[i];
      nuisance.b_ta_z[i] = B_TA[i];
    }
  }
  else if (like.IA == 4 || like.IA == 6)
  {
    nuisance.A_ia = A1[0];
    nuisance.eta_ia = A1[1];
    nuisance.oneplusz0_ia = 1.62;

    nuisance.A2_ia = A2[0];
    nuisance.eta_ia_tt = A2[1];
    nuisance.b_ta_z[0] = B_TA[0];

    for (int i=2; i<tomo.shear_Nbin; i++)
    {
      if (!(almost_equal(A1[i], 0.)) || !(almost_equal(A2[i], 0.)) || !(almost_equal(B_TA[i], 0.)))
      {
        spdlog::critical("set_nuisance_ia: one of nuisance.A_z[{}]={}, nuisance.A2_z[{}]="
          "{}, nuisance.b_ta[{}]={} was specified w/ power-law evolution\n",
          i, nuisance.A_z[i], i, nuisance.A2_z[i], i, nuisance.b_ta_z[i]);
        exit(1);
      }
    }
  }

  spdlog::debug("\x1b[90m{}\x1b[0m: Ends", "set_nuisance_ia");
}

// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// COMPUTE FUNCTIONS
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------

double cpp_compute_baryon_ratio(double log10k, double a)
{
  const double KNL = pow(10.0,log10k)*cosmology.coverH0;
  return PkRatio_baryons(KNL, a);
}

// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// RESET FUNCTIONS
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------

void cpp_reset_baryionic_struct()
{
  reset_bary_struct();
  return;
}

// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// DATA CLASS MEMBER FUNCTIONS
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------

void Data::set_mask(std::string MASK) 
{
  if (!(like.Ndata>0))
  {
    spdlog::critical("\x1b[90m{}\x1b[0m: {} not set prior to this function call",
      "set_mask", "like.Ndata");
    exit(1);
  }

  this->ndata_ = like.Ndata;
  this->mask_.set_size(this->ndata_);
  
  arma::Mat<double> table = read_table(MASK);
  for (int i=0; i<this->ndata_; i++)
  {
    this->mask_(i) = static_cast<int>(table(i,1)+1e-13);
    if(!(this->mask_(i) == 0 || this->mask_(i) == 1))
    {
      spdlog::critical("\x1b[90m{}\x1b[0m: inconsistent mask", "set_mask");
      exit(1);
    }
  }
}

void Data::set_data(std::string DATA)
{
  if (!(this->is_mask_set_))
  {
    spdlog::critical("\x1b[90m{}\x1b[0m: {} not set prior to this function call", "set_data",
      "mask");
    exit(1);
  }

  this->data_masked_.set_size(this->ndata_);
  this->data_filename_ = DATA;

  armamat table = read_table(DATA);

  for(int i=0; i<like.Ndata; i++)
  {
    this->data_masked_(i) = table(i,1);
    this->data_masked_(i) *= this->get_mask(i);
  }

  this->is_data_set_ = true;
}

void Data::set_inv_cov(std::string COV)
{
  if (!(this->is_mask_set_))
  {
    spdlog::critical("\x1b[90m{}\x1b[0m: {} not set prior to this function call", "set_inv_cov",
      "mask");
    exit(1);
  }

  armamat table = read_table(COV); // this reads cov!

  this->cov_masked_.set_size(this->ndata_, this->ndata_);
  this->cov_masked_.zeros();

  this->inv_cov_masked_.set_size(this->ndata_, this->ndata_);
  this->inv_cov_masked_.zeros();

  switch (table.n_cols)
  {
    case 3:
    {
      for (int i=0; i<static_cast<int>(table.n_rows); i++)
      {
        const int j = static_cast<int>(table(i,0));
        const int k = static_cast<int>(table(i,1));

        this->cov_masked_(j,k) = table(i,2);
        this->inv_cov_masked_(j,k) = table(i,2);

        if (j!=k)
        {
          // apply mask to off-diagonal covariance elements
          this->cov_masked_(j,k) *= this->get_mask(j);
          this->cov_masked_(j,k) *= this->get_mask(k);

          this->inv_cov_masked_(j,k) *= this->get_mask(j);
          this->inv_cov_masked_(j,k) *= this->get_mask(k);

          // m(i,j) = m(j,i)
          this->cov_masked_(k,j) = this->cov_masked_(j,k);
          this->inv_cov_masked_(k,j) = this->inv_cov_masked_(j,k);
        }
      };
      break;
    }
    case 4:
    {
      for (int i=0; i<static_cast<int>(table.n_rows); i++)
      {
        const int j = static_cast<int>(table(i,0));
        const int k = static_cast<int>(table(i,1));

        this->cov_masked_(j,k) = table(i,2) + table(i,3);
        this->inv_cov_masked_(j,k) = table(i,2) + table(i,3);

        if (j!=k)
        {
          // apply mask to off-diagonal covariance elements
          this->cov_masked_(j,k) *= this->get_mask(j);
          this->cov_masked_(j,k) *= this->get_mask(k);

          this->inv_cov_masked_(j,k) *= this->get_mask(j);
          this->inv_cov_masked_(j,k) *= this->get_mask(k);

          // m(i,j) = m(j,i)
          this->cov_masked_(k,j) = this->cov_masked_(j,k);
          this->inv_cov_masked_(k,j) = this->inv_cov_masked_(j,k);
        }
      };
      break;
    }
    case 10:
    {
      for (int i=0; i<static_cast<int>(table.n_rows); i++)
      {
        const int j = static_cast<int>(table(i,0));
        const int k = static_cast<int>(table(i,1));

        this->cov_masked_(j,k) = table(i,8) + table(i,9);
        this->inv_cov_masked_(j,k) = table(i,8) + table(i,9);

        if (j!=k)
        {
          // apply mask to off-diagonal covariance elements
          this->cov_masked_(j,k) *= this->get_mask(j);
          this->cov_masked_(j,k) *= this->get_mask(k);

          this->inv_cov_masked_(j,k) *= this->get_mask(j);
          this->inv_cov_masked_(j,k) *= this->get_mask(k);

          // m(i,j) = m(j,i)
          this->cov_masked_(k,j) = this->cov_masked_(j,k);
          this->inv_cov_masked_(k,j) = this->inv_cov_masked_(j,k);
        }
      }
      break;
    }
    default:
    {
      spdlog::critical("{}: data format for covariance file = {} is invalid", "set_inv_cov", COV);
      exit(1);
    }
  }

  this->inv_cov_masked_ = arma::inv(this->inv_cov_masked_);

  // apply mask again, to make sure numerical errors in matrix
  // inversion don't cause problems...
  // also, set diagonal elements corresponding to datavector elements
  // outside mask to zero, so that these elements don't contribute to chi2
  for (int i=0; i<this->ndata_; i++)
  {
    this->inv_cov_masked_(i,i) *= this->get_mask(i)*this->get_mask(i);
    for (int j=0; j<i; j++)
    {
      this->inv_cov_masked_(i,j) *= this->get_mask(i)*this->get_mask(j);
      this->inv_cov_masked_(j,i) = this->inv_cov_masked_(i,j);
    }
  };

  this->cov_filename_ = COV;
  this->is_inv_cov_set_ = true;
}

arma::Col<int> Data::get_mask() const
{
  return this->mask_;
}

int Data::get_mask(const int ci) const
{
  if (ci > like.Ndata || ci < 0)
  {
    spdlog::critical(
      "\x1b[90m{}\x1b[0m: index i = {} is not valid (min = {}, max = {})",
      "get_mask", ci, 0.0, like.Ndata);
    exit(1);
  }

  return this->mask_(ci);
}

int Data::get_ndim() const
{
  return this->ndata_;
}

armavec Data::get_data_masked() const
{
  return this->data_masked_;
}

double Data::get_data_masked(const int ci) const
{
  if (ci > like.Ndata || ci < 0)
  {
    spdlog::critical("\x1b[90m{}\x1b[0m: index i = {} is not valid (min = {}, max = {})",
      "get_data_masked", ci, 0, like.Ndata);
    exit(1);
  }

  return this->data_masked_(ci);
}

armamat Data::get_inverse_covariance_masked() const
{
  return this->inv_cov_masked_;
}

double Data::get_inverse_covariance_masked(const int ci, const int cj) const
{
  if (ci > like.Ndata || ci < 0)
  {
    spdlog::critical("\x1b[90m{}\x1b[0m: index i = {} is not valid (min = {}, max = {})",
      "get_inverse_covariance_masked", ci, 0.0, like.Ndata);
    exit(1);
  }
  if (cj > like.Ndata || cj < 0)
  {
    spdlog::critical("\x1b[90m{}\x1b[0m: index j = {} is not valid (min = {}, max = {})",
      "get_inverse_covariance_masked", cj, 0.0, like.Ndata);
    exit(1);
  }

  return this->inv_cov_masked_(ci, cj);
}

armamat Data::get_covariance_masked() const
{
  return this->cov_masked_;
}

double Data::get_chi2(std::vector<double> datavector) const
{
  if (!(this->is_data_set_))
  {
    spdlog::critical("\x1b[90m{}\x1b[0m: {} not set prior to this function call",
      "get_chi2", "data_vector");
    exit(1);
  }
  if (!(this->is_mask_set_))
  {
    spdlog::critical("\x1b[90m{}\x1b[0m: {} not set prior to this function call", "get_chi2",
      "mask");
    exit(1);
  }
  if (!(this->is_inv_cov_set_))
  {
    spdlog::critical("\x1b[90m{}\x1b[0m: {} not set prior to this function call", "get_chi2", 
      "inv_cov");
    exit(1);
  }

  double sum = 0.0;
  for (int i=0; i<like.Ndata; i++)
  {
    double tmp = 0.0;
    if (this->get_mask(i))
    {
      const double x = datavector[i] - this->get_data_masked(i);
      for (int j=0; j<like.Ndata; j++)
      {
        if (this->get_mask(j))
        {
          const double y = datavector[j] - this->get_data_masked(j);
          tmp += x*this->get_inverse_covariance_masked(i,j)*y;
        }
      }
    }
    sum += tmp;
  }

  if (sum < 0.0)
  {
    spdlog::critical("\x1b[90m{}\x1b[0m: chi2 = {} (invalid)", "get_chi2", sum);
    exit(1);
  }
  return sum;
}

bool Data::is_mask_set() const
{
  return this->is_mask_set_;
}

bool Data::is_data_set() const
{
  return this->is_data_set_;
}

bool Data::is_inv_cov_set() const
{
  return this->is_inv_cov_set_;
}

// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// EXPANDED DATA CLASS MEMBER FUNCTIONS
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------

void ExpandedData::set_mask(std::string MASK) 
{	
	Data::set_mask(MASK);
}

void ExpandedData::set_data(std::string DATA)
{
	Data::set_data(DATA);

  this->data_masked_reduced_dim_.set_size(this->ndata_masked_);

  for(int i=0; i<like.Ndata; i++)
  {
    if(this->get_mask(i) == 1)
    {
      if(this->get_index_reduced_dim(i) < 0)
      {
        spdlog::critical("\x1b[90m{}\x1b[0m: logical error, internal inconsistent mask operation", 
        	"set_data");
        exit(1);
      }

      this->data_masked_reduced_dim_(this->get_index_reduced_dim(i)) = this->data_masked_(i);
    }
  }
}

void ExpandedData::set_inv_cov(std::string COV)
{
	Data::set_inv_cov(COV);

  this->cov_masked_reduced_dim_.set_size(this->ndata_masked_, this->ndata_masked_);
  this->inv_cov_masked_reduced_dim_.set_size(this->ndata_masked_, this->ndata_masked_);

  for(int i=0; i<this->ndata_; i++)
  {
    for(int j=0; j<this->ndata_; j++)
    {
      if((this->mask_(i)>0.99) && (this->mask_(j)>0.99))
      {
        if(this->get_index_reduced_dim(i) < 0)
        {
          spdlog::critical("\x1b[90m{}\x1b[0m: logical error, internal inconsistent mask operation", 
          	"set_inv_cov");
          exit(1);
        }
        if(this->get_index_reduced_dim(j) < 0)
        {
          spdlog::critical("\x1b[90m{}\x1b[0m: logical error, internal inconsistent mask operation", 
          	"set_inv_cov");
          exit(1);
        }

        this->cov_masked_reduced_dim_(this->get_index_reduced_dim(i),
          this->get_index_reduced_dim(j)) = this->cov_masked_(i,j);

        this->inv_cov_masked_reduced_dim_(this->get_index_reduced_dim(i),
          this->get_index_reduced_dim(j)) = this->inv_cov_masked_(i,j);
      }
    }
  }
}

void ExpandedData::set_reduced_dim()
{
	if (!(this->is_mask_set_))
  {
    spdlog::critical(
      "\x1b[90m{}\x1b[0m: {} not set prior to this function call", "set_data",
      "mask");
    exit(1);
  }

  this->index_reduced_dim_.set_size(this->ndata_);
 
  double j=0;
  for(int i=0; i<this->ndata_; i++)
  {
    if(this->get_mask(i) > 0)
    {
      this->index_reduced_dim_(i) = j;
      j++;
    }
    else
    {
      this->index_reduced_dim_(i) = -1;
    }
  }
  if(j != this->ndata_masked_)
  {
    spdlog::critical("\x1b[90m{}\x1b[0m: logical error, internal inconsistent mask operation",
     "set_mask");
    exit(1);
  }
}

int ExpandedData::get_nreduced_dim() const
{
  return this->ndata_masked_;
}

int ExpandedData::get_index_reduced_dim(const int ci) const
{
  if (ci > like.Ndata || ci < 0)
  {
    spdlog::critical("\x1b[90m{}\x1b[0m: index i = {} is not valid"
      " (min = {}, max = {})", "get_index_reduced_dim", ci, 0.0, like.Ndata);
    exit(1);
  }

  return this->index_reduced_dim_(ci);
}

armavec ExpandedData::get_data_masked_reduced_dim() const
{
  return this->data_masked_reduced_dim_;
}

double ExpandedData::get_data_masked_reduced_dim(const int ci) const
{
  if (ci > like.Ndata || ci < 0)
  {
    spdlog::critical("\x1b[90m{}\x1b[0m: index i = {} is not valid (min = {}, max = {})",
      "get_data_masked_reduced_dim", ci, 0, like.Ndata);
    exit(1);
  }

  return this->data_masked_reduced_dim_(ci);
}

armamat ExpandedData::get_inverse_covariance_masked_reduced_dim() const
{
  return this->inv_cov_masked_reduced_dim_;
}

double ExpandedData::get_inverse_covariance_masked_reduced_dim(const int ci,
const int cj) const
{
  if (ci > like.Ndata || ci < 0)
  {
    spdlog::critical("\x1b[90m{}\x1b[0m: index i = {} is not valid (min = {}, max = {})",
      "get_inverse_covariance_masked_reduced_dim", ci, 0.0, like.Ndata);
    exit(1);
  }
  if (cj > like.Ndata || cj < 0)
  {
    spdlog::critical("\x1b[90m{}\x1b[0m: index j = {} is not valid (min = {}, max = {})",
      "get_inverse_covariance_masked_reduced_dim", cj, 0.0, like.Ndata);
    exit(1);
  }

  return this->inv_cov_masked_reduced_dim_(ci, cj);
}

armamat ExpandedData::get_covariance_masked_reduced_dim() const
{
  return this->cov_masked_reduced_dim_;
}

armavec ExpandedData::get_expand_dim_from_masked_reduced_dim(armavec reduced_dim_vector) const
{
  if (this->ndata_masked_ != static_cast<int>(reduced_dim_vector.n_elem))
  {
    spdlog::critical("\x1b[90m{}\x1b[0m: {} invalid input vector",
      "get_expand_dim_from_masked_reduced_dim"
    );
    exit(1);
  }

  armavec vector;
  vector.set_size(this->ndata_);
  vector.zeros();

  for(int i=0; i<this->ndata_; i++)
  {
    if(this->mask_(i) > 0.99)
    {
      if(this->get_index_reduced_dim(i) < 0)
      {
        spdlog::critical("\x1b[90m{}\x1b[0m: logical error, internal inconsistent mask operation",
        	"get_expand_dim_from_masked_reduced_dim");
        exit(1);
      }
      vector(i) = reduced_dim_vector(this->get_index_reduced_dim(i));
    }
  }

  return vector;
}

// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// BaryonScenario CLASS MEMBER FUNCTIONS
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------

int BaryonScenario::nscenarios() const
{
  return this->nscenarios_;
}

void BaryonScenario::set_scenarios(std::string scenarios)
{
  std::vector<std::string> lines;
  lines.reserve(50);

  // Second: Split file into lines
  boost::trim_if(scenarios, boost::is_any_of("\t "));
  boost::trim_if(scenarios, boost::is_any_of("\n"));

  if (scenarios.empty())
  {
    spdlog::critical("\x1b[90m{}\x1b[0m: invalid string input (empty)","init_baryon_pca_scenarios");
    exit(1);
  }

  boost::split(lines, scenarios, boost::is_any_of("/"), boost::token_compress_on);

  this->nscenarios_ = lines.size();

  for(int i=0; i<this->nscenarios_; i++)
  {
  	//this->scenarios_.insert_or_assign(std::pair<int,std::string>(i,lines.at(i)));
    this->scenarios_[i] = lines.at(i);
  }

  return;
}

std::string BaryonScenario::get_scenario(const int i) const
{
  return this->scenarios_.at(i);
}

// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// PointMass CLASS MEMBER FUNCTIONS 
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------

void PointMass::set_pm_vector(std::vector<double> pm)
{
  this->pm_ = pm;
  return;
}

std::vector<double> PointMass::get_pm_vector() const
{
  return this->pm_;
}

double PointMass::get_pm(const int zl, const int zs, const double theta) const
{
  constexpr double G_over_c2 = 1.6e-23;
  const double a_lens = 1.0/(1.0 + zmean(zl));
  const double chi_lens = chi(a_lens);

  return 4*G_over_c2*this->pm_[zl]*1.e+13*g_tomo(a_lens, zs)/(theta*theta)/
    (chi_lens*a_lens);
}

// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// NON MEMBER FUNCIONS 
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------

armamat read_table(const std::string file_name)
{
  std::ifstream input_file(file_name);

  if (!input_file.is_open())
  {
    spdlog::critical("\x1b[90m{}\x1b[0m: file {} cannot be opened", "read_table", file_name);
    exit(1);
  }

  // Read the entire file into memory
  std::string tmp;
  input_file.seekg(0,std::ios::end);
  tmp.resize(static_cast<size_t>(input_file.tellg()));
  input_file.seekg(0,std::ios::beg);
  input_file.read(&tmp[0],tmp.size());
  input_file.close();
  if(tmp.empty())
  {
    spdlog::critical("\x1b[90m{}\x1b[0m: file {} is empty", "read_table", file_name);
    exit(1);
  }
  std::vector<std::string> lines;
  lines.reserve(50000);
  // Second: Split file into lines
  boost::trim_if(tmp,boost::is_any_of("\t "));
  boost::trim_if(tmp,boost::is_any_of("\n"));
  boost::split(lines, tmp,boost::is_any_of("\n"), boost::token_compress_on);
  // Erase comment/blank lines
  auto check = [](std::string mystr) -> bool
  {
    return boost::starts_with(mystr, "#");
  };
  lines.erase(std::remove_if(lines.begin(), lines.end(), check), lines.end());
  // Third: Split line into words
  armamat result;
  size_t ncols = 0;
  { // first line
    std::vector<std::string> words;
    words.reserve(100);
    boost::split(words,lines[0], boost::is_any_of(" \t"), boost::token_compress_on);
    ncols = words.size();
    result.set_size(lines.size(), ncols);
    for (size_t j=0; j<ncols; j++)
    {
      result(0,j) = std::stod(words[j]);
    }
  }
  #pragma omp parallel for
  for (size_t i=1; i<lines.size(); i++)
  {
    std::vector<std::string> words;
    boost::split(words, lines[i], boost::is_any_of(" \t"), boost::token_compress_on);
    if (words.size() != ncols)
    {
    	#pragma omp critical
      {
      	spdlog::critical("\x1b[90m{}\x1b[0m: file {} is not well formatted (regular table required)"
      	"read_table", file_name);
      	exit(1);
      }
    }
    for (size_t j=0; j<ncols; j++)
    {
      result(i,j) = std::stod(words[j]);
    }
  };
  return result;
}

dvec convert_arma_col_to_stl_vector(armavec in)
{
  dvec out(in.n_elem, 0.0);

  for(int i=0; i<static_cast<int>(in.n_elem); i++)
  {
    out[i] = in(i);
  }

  return out;
}

#ifdef ADDCLUSTER_BASICINTERFACE

void cpp_init_cluster_N200(int N200bins, std::vector<double> N200_edges)
{
  spdlog::debug("\x1b[90m{}\x1b[0m: Begins", "init_cluster_N200");

  if (!(N200bins > 0) || N200bins > MAX_SIZE_ARRAYS)
  {
    spdlog::critical("\x1b[90m{}\x1b[0m: {} = {} not supported (max = {})",
      "init_cluster_N200", "N200bins", N200bins, MAX_SIZE_ARRAYS);
    exit(1);
  }
  if (static_cast<int>(N200_edges.size()) != N200bins+1)
  {
    spdlog::critical("\x1b[90m{}\x1b[0m: incompatible {} array size ({} = {})",
      "init_cluster_N200", "N200_edges", "N200bins", N200bins);
    exit(1);
  }

  Cluster.N200_Nbin = N200bins;
  
  for (int i=0; i<Cluster.N200_Nbin; i++)
  {
    Cluster.N_min[i] = N200_edges[i];
    Cluster.N_max[i] = N200_edges[i+1];    

    spdlog::debug("\x1b[90m{}\x1b[0m: Cluster lambda-bin {}: [{}, {}]", "init_cluster_N200", i,
      Cluster.N_min[i], Cluster.N_max[i]);
  }

  spdlog::debug("\x1b[90m{}\x1b[0m: Cluster.N200_Nbin = {}", "init_cluster_N200", 
    Cluster.N200_Nbin);

  spdlog::debug("\x1b[90m{}\x1b[0m: Ends", "init_cluster_N200");
}

void cpp_init_cluster_tomo(int Ntomo, std::vector<double> zbin_edges) // it does not set w_cg! 
{
  spdlog::debug("\x1b[90m{}\x1b[0m: Begins", "init_tomo_cluster");
  
  if (!(Ntomo > 0) || Ntomo > MAX_SIZE_ARRAYS)
  {
    spdlog::critical("\x1b[90m{}\x1b[0m: {} = {} not supported (max = {})",
      "init_cluster_tomo", "Ntomo", Ntomo, MAX_SIZE_ARRAYS);
    exit(1);
  }
  if (static_cast<int>(zbin_edges.size()) != Ntomo + 1)
  {
    spdlog::critical("\x1b[90m{}\x1b[0m: incompatible {} array size ({} = {})",
      "init_cluster_tomo", "zbin_edges", "Ntomo", Ntomo);
    exit(1);
  }
  if (tomo.shear_Nbin == 0)
  {
    spdlog::critical("{}: {} not set prior to this function call", "init_cluster_tomo", 
      "tomo.shear_Nbin");
    exit(1);
  }

  tomo.cluster_Nbin = Ntomo;
  tomo.cc_clustering_Npowerspectra = Ntomo; // no support for cross-z bins

  for (int i=0; i<Ntomo; i++)
  {
    tomo.cluster_zmin[i] = zbin_edges[i];
    tomo.cluster_zmax[i] = zbin_edges[i+1];
    spdlog::debug("\x1b[90m{}\x1b[0m: Cluster z-bin {}: [{}, {}]", "init_tomo_cluster", i,
      Cluster.N_min[i], Cluster.N_max[i]);
  }
  
  int N =0;
  for (int j=0; j<tomo.cluster_Nbin; j++)
  {
    for (int i=0; i<tomo.shear_Nbin; i++)
    {
      if (test_zoverlap_c(j, i)) 
      {
        N++;
      }
    }
  }
  tomo.cgl_Npowerspectra = N;
  
  spdlog::debug("\x1b[90m{}\x1b[0m: Ends", "init_tomo_cluster");
}

// right now we are requiring external input on w_cg bin combination
void cpp_init_ext_choice_wcg_tomo(arma::Mat<int> bin_combination)
{
  const int nrows = static_cast<int>(bin_combination.n_rows);  // cluster bins
  const int ncols = static_cast<int>(bin_combination.n_cols); // galaxy bins
  if (ncols != 2)
  {
    spdlog::critical("\x1b[90m{}\x1b[0m: invalid matrix input", "init_ext_choice_wcg_tomo");
    exit(1);
  }

  tomo.cg_clustering_Npowerspectra = nrows; 
  tomo.external_selection_cg_clustering = (int**) malloc(ncols*sizeof(int*));
  tomo.external_selection_cg_clustering[0] = (int*) malloc(nrows*sizeof(int));
  tomo.external_selection_cg_clustering[1] = (int*) malloc(nrows*sizeof(int));
  
  for (int i=0; i<nrows; i++)
  {
    tomo.external_selection_cg_clustering[0][i] = bin_combination(i, 0); // cluster
    tomo.external_selection_cg_clustering[1][i] = bin_combination(i, 1); // galaxy
    if(bin_combination(i, 0) > tomo.cluster_Nbin - 1)
    {
      spdlog::critical("\x1b[90m{}\x1b[0m: invalid cluster redshift bin input ni = {} (max = {})",
        "init_ext_choice_wcg_tomo", bin_combination(i, 0), tomo.cluster_Nbin-1);
      exit(1);
    }
    if(bin_combination(i, 1) > tomo.clustering_Nbin - 1)
    {
      spdlog::critical("\x1b[90m{}\x1b[0m: invalid galaxy redshift bin input ni = {} (max = {})",
        "init_ext_choice_wcg_tomo", bin_combination(i, 1), tomo.clustering_Nbin-1);
      exit(1);
    }
  }
}

void set_mass_observable_relation(std::vector<double> MOR, std::string survey)
{
  spdlog::debug("\x1b[90m{}\x1b[0m: Begins", "set_mass_observable_relation");
  
  boost::to_upper(survey);

  if(survey.compare("SDSS") == 0) 
  {
      if(MOR.size() < 5) 
      {
        spdlog::critical("\x1b[90m{}\x1b[0m: MOR size {} = invalid input size",
          "set_mass_observable_relation", MOR.size());
        exit(1);
      }

      nuisance.cluster_MOR[0] = MOR[0]; // Mmin
      nuisance.cluster_MOR[1] = MOR[1]; // M1
      nuisance.cluster_MOR[2] = MOR[2]; // alpha
      nuisance.cluster_MOR[3] = MOR[3]; // sigma_int
      nuisance.cluster_MOR[4] = MOR[4]; // epsilon
      nuisance.N_cluster_MOR = 5;
  }
  else if(survey.compare("BUZZARD") == 0) 
  {
    if(MOR.size() < 4) 
    {
      spdlog::critical("\x1b[90m{}\x1b[0m: MOR size {} = invalid input size",
        "set_mass_observable_relation", MOR.size());
      exit(1);
    }

    nuisance.cluster_MOR[0] = MOR[0];
    nuisance.cluster_MOR[1] = MOR[1];
    nuisance.cluster_MOR[2] = MOR[2];
    nuisance.cluster_MOR[3] = MOR[3];
    nuisance.N_cluster_MOR = 4;
  }
  else
  {
    spdlog::critical("\x1b[90m{}\x1b[0m: survey = {} invalid survey option",
      "set_mass_observable_relation", survey);
    exit(1);
  }

  spdlog::debug("\x1b[90m{}\x1b[0m: Ends", "set_mass_observable_relation");
}

#endif

} // namespace basic_interface