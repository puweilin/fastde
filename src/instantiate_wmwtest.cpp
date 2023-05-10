#include "fastde/wmwtest.tpp"
#include "cpp11/doubles.hpp"
#include "cpp11/integers.hpp"
#include "cpp11/matrix.hpp"

template void sum_rank(std::vector<std::pair<double, int>> const & temp, 
  std::unordered_map<int, size_t> const & z_cl_counts,
  size_t const & count, double const & zero_val,
  std::unordered_map<int, size_t> & rank_sums,
  double & tie_sum
);

// ids points to start of the column's positive element row ids
// x  points to the start fo the columns positive element values
// count is hte number of positive elements in the column.
template void spmat_sort(double * in, int * ids,
  size_t const & nz_count, 
  int * labels, size_t const & count, double const & zero_val,
  std::vector<std::pair<int, size_t> > const & cl_counts,
  std::vector<std::pair<double, int> > & temp,
  std::unordered_map<int, size_t> & z_cl_counts) ;

// rank_sums output:  map cluster to rank_sum.
// counts zeros and skip sorting of those.  this is faster than sorting all.
template void pseudosparse_sort(
  double * in, int * labels, size_t const & count, double const & zero_val,
  std::vector<std::pair<int, size_t> > const & cl_counts,
  std::vector<std::pair<double, int> > & temp,
  std::unordered_map<int, size_t> & z_cl_counts);



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
//VEC
// direct write to doubles_matrix is not fast for cpp11:  proxy object creation and iterator creation....


template void spmat_sort_vec(
  cpp11::doubles const & in, cpp11::integers const & ids,
  size_t const & nz_offset, 
  size_t const & nz_count, 
  cpp11::integers const & labels, size_t const & count, double const & zero_val,
  std::vector<std::pair<int, size_t> > const & cl_counts,
  std::vector<std::pair<double, int> > & temp,
  std::unordered_map<int, size_t> & z_cl_counts);

template void pseudosparse_sort_vec(
  cpp11::doubles_matrix<cpp11::by_column>::slice const & in, size_t const & offset, 
  cpp11::integers const & labels, size_t const & count, double const & zero_val,
  std::vector<std::pair<int, size_t> > const & cl_counts,
  std::vector<std::pair<double, int> > & temp,
  std::unordered_map<int, size_t> & z_cl_counts);

template void pseudosparse_sort_vec(
  std::vector<double> const & in, size_t const & offset, 
  cpp11::integers const & labels, size_t const & count, double const & zero_val,
  std::vector<std::pair<int, size_t> > const & cl_counts,
  std::vector<std::pair<double, int> > & temp,
  std::unordered_map<int, size_t> & z_cl_counts);


void wmw_vec(
  std::vector<std::pair<int, size_t> > const & cl_counts, 
  std::unordered_map<int, size_t> const & rank_sums, 
  double const & tie_sum,
  size_t const & count,
  std::vector<double> & out, size_t const & offset,
  int const & test_type, bool const & continuity);

void wmw_vec(
  std::vector<std::pair<int, size_t> > const & cl_counts, 
  std::unordered_map<int, size_t> const & rank_sums, 
  double const & tie_sum,
  size_t const & count,
  cpp11::writable::doubles_matrix<cpp11::by_column>::slice & out, size_t const & offset,
  int const & test_type, bool const & continuity);


template void vecsc_wmw_vecsc(
  std::vector<double> const & mat, size_t const & nsamples, size_t const & nfeatures,
  cpp11::integers const & lab, 
  int const & rtype, 
  bool const & continuity_correction, 
  std::vector<double> & pv,
  std::vector<std::pair<int, size_t> > const & sorted_cluster_counts,
  int const & threads);

template void matc_wmw_vecsc(
  cpp11::doubles_matrix<cpp11::by_column> const & mat, size_t const & nsamples, 
  size_t const & nfeatures,
  cpp11::integers const & lab, 
  int const & rtype, 
  bool const & continuity_correction, 
  std::vector<double> & pv,
  std::vector<std::pair<int, size_t> > const & sorted_cluster_counts,
  int const & threads);

// template void matc_wmw_matc(
//   cpp11::doubles_matrix<cpp11::by_column> const & mat, size_t const & nsamples, 
//   size_t const & nfeatures,
//   cpp11::integers const & lab, 
//   int const & rtype, 
//   bool const & continuity_correction, 
//   cpp11::writable::doubles_matrix<cpp11::by_column> & pv,
//   std::vector<std::pair<int, size_t> > const & sorted_cluster_counts,
//   int const & threads);


template void csc_wmw_vecsc(
  cpp11::doubles const & x, cpp11::integers const & i, cpp11::doubles const & p, 
  size_t nsamples, size_t nfeatures,
  cpp11::integers const & lab,
  int const & rtype, 
  bool const & continuity_correction, 
  std::vector<double> & pv,
  std::vector<std::pair<int, size_t> > const & sorted_cluster_counts,
  int const & threads);

template void csc_wmw_vecsc(
  cpp11::doubles const & x, cpp11::integers const & i, cpp11::integers const & p, 
  size_t nsamples, size_t nfeatures,
  cpp11::integers const & lab,
  int const & rtype, 
  bool const & continuity_correction, 
  std::vector<double> & pv,
  std::vector<std::pair<int, size_t> > const & sorted_cluster_counts,
  int const & threads);

// template void csc_wmw_matc(
//   cpp11::doubles const & x, cpp11::integers const & i, cpp11::doubles const & p, 
//   size_t nsamples, size_t nfeatures,
//   cpp11::integers const & lab,
//   int const & rtype, 
//   bool const & continuity_correction, 
//   cpp11::writable::doubles_matrix<cpp11::by_column> & pv,
//   std::vector<std::pair<int, size_t> > const & sorted_cluster_counts,
//   int const & threads);

// template void csc_wmw_matc(
//   cpp11::doubles const & x, cpp11::integers const & i, cpp11::integers const & p, 
//   size_t nsamples, size_t nfeatures,
//   cpp11::integers const & lab,
//   int const & rtype, 
//   bool const & continuity_correction, 
//   cpp11::writable::doubles_matrix<cpp11::by_column> & pv,
//   std::vector<std::pair<int, size_t> > const & sorted_cluster_counts,
//   int const & threads);


template void csc_wmw_vecsc(
  std::vector<double> const & x, std::vector<int> const & i, std::vector<int> const & p, 
  size_t nsamples, size_t nfeatures,
  cpp11::integers const & lab,
  int const & rtype, 
  bool const & continuity_correction, 
  std::vector<double> & pv,
  std::vector<std::pair<int, size_t> > const & sorted_cluster_counts,
  int const & threads);

template void csc_wmw_vecsc(
  std::vector<double> const & x, std::vector<int> const & i, std::vector<long> const & p, 
  size_t nsamples, size_t nfeatures,
  cpp11::integers const & lab,
  int const & rtype, 
  bool const & continuity_correction, 
  std::vector<double> & pv,
  std::vector<std::pair<int, size_t> > const & sorted_cluster_counts,
  int const & threads);

// template void csc_wmw_matc(
//   std::vector<double> const & x, std::vector<int> const & i, std::vector<int> const & p, 
//   size_t nsamples, size_t nfeatures,
//   cpp11::integers const & lab,
//   int const & rtype, 
//   bool const & continuity_correction, 
//   cpp11::writable::doubles_matrix<cpp11::by_column> & pv,
//   std::vector<std::pair<int, size_t> > const & sorted_cluster_counts,
//   int const & threads);

// template void csc_wmw_matc(
//   std::vector<double> const & x, std::vector<int> const & i, std::vector<long> const & p, 
//   size_t nsamples, size_t nfeatures,
//   cpp11::integers const & lab,
//   int const & rtype, 
//   bool const & continuity_correction, 
//   cpp11::writable::doubles_matrix<cpp11::by_column> & pv,
//   std::vector<std::pair<int, size_t> > const & sorted_cluster_counts,
//   int const & threads);