#include "fastde/wmwtest.tpp"
#include "cpp11/doubles.hpp"
#include "cpp11/integers.hpp"
#include "cpp11/matrix.hpp"

template void sparse_sum_rank(std::pair<double, int> * temp, size_t const & nzc, 
  size_t const & count, double const & zero_val,
  std::unordered_map<int, size_t> const & z_cl_counts,
  std::unordered_map<int, size_t> & rank_sums,
  double & tie_sum
);

// ids points to start of the column's positive element row ids
// x  points to the start fo the columns positive element values
// count is hte number of positive elements in the column.
template void sparse_wmw_summary(double * in, int * ids,
  size_t const & nz_count, 
  int * labels, size_t const & count, double const & zero_val,
  std::vector<std::pair<int, size_t> > const & cl_counts,
  std::unordered_map<int, size_t> & rank_sums, double & tie_sum) ;

// rank_sums output:  map cluster to rank_sum.
// counts zeros and skip sorting of those.  this is faster than sorting all.
template void pseudosparse_wmw_summary(
  double * in, int * labels, size_t const & count, double const & zero_val,
  std::vector<std::pair<int, size_t> > const & cl_counts,
  std::unordered_map<int, size_t> & rank_sums, double & tie_sum);



// types:  0 = less, 1 = greater, 2 = twosided (default), 3 = U2
template void wmw(
  std::vector<std::pair<int, size_t> > const & cl_counts, 
  std::unordered_map<int, size_t> const & rank_sums, 
  double const & tie_sum,
  size_t const & count,
  double * out,
  int const & test_type, bool const & continuity) ;


template void omp_dense_wmw(
    double * mat, size_t const & nsamples, size_t const & nfeatures,
    int * lab, 
    int rtype, 
    bool continuity_correction, 
    std::vector<double> &pv,
    std::vector<std::pair<int, size_t> > &sorted_cluster_counts,
    int threads);


// =================================
template void omp_sparse_wmw(
    double * x, int * i, int * p, size_t nsamples, size_t nfeatures,
    int * lab,
    int rtype, 
    bool continuity_correction, 
    std::vector<double> &pv,
    std::vector<std::pair<int, size_t> > &sorted_cluster_counts,
    int threads);

template void omp_sparse_wmw(
    double * x, int * i, long * p, size_t nsamples, size_t nfeatures,
    int * lab,
    int rtype, 
    bool continuity_correction, 
    std::vector<double> &pv,
    std::vector<std::pair<int, size_t> > &sorted_cluster_counts,
    int threads);

// =========================================

template void sparse_wmw_summary_vec(
  cpp11::doubles const & in, cpp11::integers const & ids,
  size_t const & nz_offset, 
  size_t const & nz_count, 
  cpp11::integers const & labels, size_t const & count, double const & zero_val,
  std::vector<std::pair<int, size_t> > const & cl_counts,
  std::unordered_map<int, size_t> & rank_sums, double & tie_sum) ;

template void pseudosparse_wmw_summary_vec(
  cpp11::doubles_matrix<cpp11::by_column>::slice const & in, size_t const & offset, 
  cpp11::integers const & labels, size_t const & count, double const & zero_val,
  std::vector<std::pair<int, size_t> > const & cl_counts,
  std::unordered_map<int, size_t> & rank_sums, double & tie_sum);

template void pseudosparse_wmw_summary_vec(
  std::vector<double> const & in, size_t const & offset, 
  cpp11::integers const & labels, size_t const & count, double const & zero_val,
  std::vector<std::pair<int, size_t> > const & cl_counts,
  std::unordered_map<int, size_t> & rank_sums, double & tie_sum);


void wmw_vec(
  std::vector<std::pair<int, size_t> > const & cl_counts, 
  std::unordered_map<int, size_t> const & rank_sums, 
  double const & tie_sum,
  size_t const & count,
  std::vector<double> & out, size_t const & offset, size_t const & label_count,
  int const & test_type, bool const & continuity);

template void csc_dense_wmw_vec(
  std::vector<double> const & mat, size_t const & nsamples, size_t const & nfeatures,
  cpp11::integers const & lab, 
  int const & rtype, 
  bool const & continuity_correction, 
  std::vector<double> & pv,
  std::vector<std::pair<int, size_t> > const & sorted_cluster_counts,
  int const & threads);

template void csc_dense_wmw_mat(
  cpp11::doubles_matrix<cpp11::by_column> const & mat, size_t const & nsamples, 
  size_t const & nfeatures,
  cpp11::integers const & lab, 
  int const & rtype, 
  bool const & continuity_correction, 
  std::vector<double> & pv,
  std::vector<std::pair<int, size_t> > const & sorted_cluster_counts,
  int const & threads);


template void csc_sparse_wmw_vec(
  cpp11::doubles const & x, cpp11::integers const & i, cpp11::doubles const & p, 
  size_t nsamples, size_t nfeatures,
  cpp11::integers const & lab,
  int const & rtype, 
  bool const & continuity_correction, 
  std::vector<double> & pv,
  std::vector<std::pair<int, size_t> > const & sorted_cluster_counts,
  int const & threads);

template void csc_sparse_wmw_vec(
  cpp11::doubles const & x, cpp11::integers const & i, cpp11::integers const & p, 
  size_t nsamples, size_t nfeatures,
  cpp11::integers const & lab,
  int const & rtype, 
  bool const & continuity_correction, 
  std::vector<double> & pv,
  std::vector<std::pair<int, size_t> > const & sorted_cluster_counts,
  int const & threads);

  template void csc_sparse_wmw_vec(
  std::vector<double> const & x, std::vector<int> const & i, std::vector<int> const & p, 
  size_t nsamples, size_t nfeatures,
  cpp11::integers const & lab,
  int const & rtype, 
  bool const & continuity_correction, 
  std::vector<double> & pv,
  std::vector<std::pair<int, size_t> > const & sorted_cluster_counts,
  int const & threads);

    template void csc_sparse_wmw_vec(
  std::vector<double> const & x, std::vector<int> const & i, std::vector<long> const & p, 
  size_t nsamples, size_t nfeatures,
  cpp11::integers const & lab,
  int const & rtype, 
  bool const & continuity_correction, 
  std::vector<double> & pv,
  std::vector<std::pair<int, size_t> > const & sorted_cluster_counts,
  int const & threads);