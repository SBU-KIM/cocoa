#define ARMA_DONT_USE_WRAPPER
#include <armadillo>
#include <map>
#include <carma/carma.h>

#ifndef __COSMOLIKE_BASIC_INTERFACE_HPP
#define __COSMOLIKE_BASIC_INTERFACE_HPP

namespace basic_interface
{

// ------------------------------------------------------------------------
// ------------------------------------------------------------------------
// ------------------------------------------------------------------------
// ------------------------------------------------------------------------

class RandomNumber
{// Singleton Class that holds a random number generator
public:
  static RandomNumber& get_instance()
  {
 	static RandomNumber instance;
	return instance;
  }
  ~RandomNumber() = default;

  double get()
  {
    return dist_(mt_);
  }
protected:
  std::random_device rd_;
  std::mt19937 mt_;
  std::uniform_real_distribution<double> dist_;
private:
  RandomNumber() :
    rd_(),
    mt_(rd_()),
    dist_(0.0, 1.0) {
	};
  RandomNumber(RandomNumber const&) = delete;
};

// ------------------------------------------------------------------------
// ------------------------------------------------------------------------
// ------------------------------------------------------------------------
// ------------------------------------------------------------------------

class Data
{ 
public:
  virtual void set_data(std::string DATA);

  virtual void set_mask(std::string MASK);

  virtual void set_inv_cov(std::string COV);

  // ----- const -----------------------------------------------------

  int get_ndim() const;

  arma::Col<int> get_mask() const;

  int get_mask(const int ci) const;

  arma::Col<double> get_data_masked() const;

  double get_data_masked(const int ci) const;

  arma::Mat<double> get_covariance_masked() const;

  arma::Mat<double> get_inverse_covariance_masked() const;

  double get_inverse_covariance_masked(const int ci, const int cj) const;

  double get_chi2(std::vector<double> datavector) const;

  bool is_mask_set() const;

  bool is_data_set() const;

  bool is_inv_cov_set() const;
protected:
  bool is_mask_set_ = false;
  bool is_data_set_ = false;
  bool is_inv_cov_set_ = false;

  int ndata_;
  int ndata_masked_;

  std::string mask_filename_;
  std::string cov_filename_;
  std::string data_filename_;

  arma::Col<int> mask_;
  arma::Col<double> data_masked_;

  arma::Mat<double> cov_masked_;
  arma::Mat<double> inv_cov_masked_;
};

// ------------------------------------------------------------------------
// ------------------------------------------------------------------------
// ------------------------------------------------------------------------
// ------------------------------------------------------------------------

class ExpandedData : public Data
{
public:
  virtual void set_data(std::string DATA);

  virtual void set_mask(std::string MASK);

  virtual void set_inv_cov(std::string COV);

  void set_reduced_dim();

  // ----- const -----------------------------------------------------

  int get_nreduced_dim() const;

  int get_index_reduced_dim(const int ci) const;

  arma::Col<double> get_data_masked_reduced_dim() const;

  double get_data_masked_reduced_dim(const int ci) const;

  arma::Mat<double> get_covariance_masked_reduced_dim() const;

  arma::Mat<double> get_inverse_covariance_masked_reduced_dim() const;

  double get_inverse_covariance_masked_reduced_dim(const int ci, const int cj) const;

  arma::Col<double> get_expand_dim_from_masked_reduced_dim(
  	arma::Col<double> reduced_dim_vector) const;
protected:
  arma::Col<int> index_reduced_dim_;
  arma::Col<double> data_masked_reduced_dim_; // for baryon project, reduced dim
  arma::Mat<double> cov_masked_reduced_dim_;  // for baryon project, reduced dim
  arma::Mat<double> inv_cov_masked_reduced_dim_;
};

// ------------------------------------------------------------------------
// ------------------------------------------------------------------------
// ------------------------------------------------------------------------
// ------------------------------------------------------------------------

class PointMass
{
public:
  void set_pm_vector(std::vector<double> pm);

  std::vector<double> get_pm_vector() const;

  double get_pm(const int zl, const int zs, const double theta) const;
protected:
  std::vector<double> pm_;
};

// ------------------------------------------------------------------------
// ------------------------------------------------------------------------
// ------------------------------------------------------------------------
// ------------------------------------------------------------------------

class BaryonScenario
{
public:
  int nscenarios() const;

  void set_scenarios(std::string scenarios);

  std::string get_scenario(const int i) const;
protected:
  int nscenarios_;
  std::map<int, std::string> scenarios_;
};

// ------------------------------------------------------------------------
// ------------------------------------------------------------------------
// ------------------------------------------------------------------------
// ------------------------------------------------------------------------

arma::Mat<double> read_table(const std::string file_name);

std::vector<double> convert_arma_col_to_stl_vector(arma::Col<double> in);

// https://en.cppreference.com/w/cpp/types/numeric_limits/epsilon
template<class T>
typename std::enable_if<!std::numeric_limits<T>::is_integer, bool>::type
almost_equal(T x, T y, int ulp = 100)
{
    // the machine epsilon has to be scaled to the magnitude of the values used
    // and multiplied by the desired precision in ULPs (units in the last place)
    return std::fabs(x-y) <= std::numeric_limits<T>::epsilon() * std::fabs(x+y) * ulp
        // unless the result is subnormal
        || std::fabs(x-y) < std::numeric_limits<T>::min();
}

// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// INIT FUNCTIONS
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------

void cpp_initial_setup();

void cpp_init_survey(std::string surveyname, double area, double sigma_e);

void cpp_init_cosmo_runmode(const bool is_linear);

void cpp_init_binning_real(const int Ntheta, const double thetamin_arcmin, 
	const double thetamax_arcmin);

void cpp_init_binning_fourier(const int Ncl, const double lmin, const double lmax);

void cpp_init_lens_sample(std::string multihisto_file, const int Ntomo, const double ggl_cut);

void cpp_init_source_sample(std::string multihisto_file, const int Ntomo);

void cpp_init_linear_power_spectrum(std::vector<double> io_log10k, std::vector<double> io_z, 
	std::vector<double> io_lnP);

void cpp_init_non_linear_power_spectrum(std::vector<double> io_log10k, std::vector<double> io_z, 
	std::vector<double> io_lnP);

void cpp_init_IA(const int N);

void cpp_init_growth(std::vector<double> io_z, std::vector<double> io_G);

void cpp_init_distances(std::vector<double> io_z, std::vector<double> io_chi);

void cpp_init_probes(std::string possible_probes);

void cpp_init_size_data_vector_real(const bool shear_shear, const bool shear_pos, 
const bool pos_pos, const bool ks, const bool gk, const bool kk, const bool clusterCC,
const bool clusterN, const bool clusterWL, const bool clusterCG);

void cpp_init_size_data_vector_fourier(const bool shear_shear, const bool shear_pos, 
const bool pos_pos, const bool ks, const bool gk, const bool kk);

void cpp_reset_baryionic_struct();

void cpp_init_baryons_contamination(const bool use_baryonic_simulations_contamination,
const std::string which_baryonic_simulations_contamination);

// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// SET PARAM FUNCTIONS
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------

void cpp_set_cosmological_parameters(const double OM, const double H0, const bool is_cached_cosmo);

void cpp_set_nuisance_shear_calib(std::vector<double> M);

void cpp_set_nuisance_shear_photoz(std::vector<double> SP);

void cpp_set_nuisance_clustering_photoz(std::vector<double> CP);

void cpp_set_nuisance_linear_bias(std::vector<double> B1);

void cpp_set_nuisance_nonlinear_bias(std::vector<double> B1, std::vector<double> B2);

void cpp_set_nuisance_magnification_bias(std::vector<double> B_MAG);

void cpp_set_nuisance_ia(std::vector<double> A1, std::vector<double> A2,
	std::vector<double> B_TA);

void cpp_set_nuisance_bias(std::vector<double> B1, std::vector<double> B2,
	std::vector<double> B_MAG);

// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// COMPUTE FUNCTIONS
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------

double cpp_compute_baryon_ratio(double log10k, double a);

}  // namespace interface_mpp_aux
#endif // HEADER GUARD