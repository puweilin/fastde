#include <Rcpp.h>
#include <R.h>
#include <Rdefines.h>
using namespace Rcpp;


/*
 *  copied and modified from wave/include/transform/wilcox.hpp
 *
 *  Created on: Jun 8, 2021
 *  Author: Tony Pan
 *  Affiliation: Institute for Data Engineering and Science
 *  			Georgia Institute of Technology, Atlanta, GA 30332
 */

#include <cmath>  // sqrt
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <algorithm>

// Enable C++11
// [[Rcpp::plugins(cpp11)]]

// Enable OpenMP (excludes macOS)
// [[Rcpp::plugins(openmp)]]

#include <omp.h>

// #ifndef INV_SQRT2
// #define INV_SQRT2 0.70710678118
// #endif

#ifndef NROW
#define NROW(x) INTEGER(GET_DIM((x)))[0]
#endif

#ifndef NCOL
#define NCOL(x) INTEGER(GET_DIM((x)))[1]
#endif

template <typename IT, typename LABEL>
void dense_wmw_summary(IT const * in, LABEL const * labels, size_t const & count, 
  std::unordered_map<LABEL, std::pair<size_t, size_t>> & rank_sums, double & tie_sum) {

  // tuplize in and labels.
  if (count == 0) return;
  assert((count < 20) && "ERROR: count is too small (< 20) for a normal approximation\n");
  // if (count < 100) {
  //     printf("WARNING: count is small for a normal approximation: %ld\n", count);
  // }

  // ============== populate the sort structure and Sort
  std::vector<std::pair<IT, LABEL>> temp;
  temp.reserve(count);
  for (size_t i = 0; i < count; ++i) {
      temp.emplace_back(in[i], labels[i]);
  }

  // sort by label...
  std::sort(temp.begin(), temp.end(), [](std::pair<IT, LABEL> const & a, std::pair<IT, LABEL> const & b){
      return a.first < b.first;
  });

  // assign rank and accumulate count and rank
  rank_sums.clear();  // random access of this...

  // to handle ties, we need 0.5.  To make comparison easier and not float based, double the rank value here and divide later.
  // actually, we may need to have median of the tie ranks.
  // here is an approach:  keep undoubled rank, walk once from front, and once from back.  if ties, do not change rank until value changes.  each sum /2 would give correct rank
  //    at end of tie, reset to use array index as rank.
  // any tie detection would need 1x iteration.   each tie would need another pass to get median so may break prefetching, but this may not be common?
  // should minimize writes in the big array as well.

  size_t rank = 1;
  double tie;
  tie_sum = 0.0;
  LABEL key = temp[0].second;
  // insert first
  if (rank_sums.count(key) == 0) rank_sums[key] = std::make_pair(1, rank);
  else {
      ++rank_sums[key].first;
      rank_sums[key].second += rank;
  }
  // forward walk
  for (size_t i = 1; i < count; ++i) {
      // rank update
      if (temp[i-1].first < temp[i].first) {  // ties are dealt with only when exiting ties.
          tie = i + 1 - rank;  // rank is of the first element in the tie.
          tie_sum += tie * tie * tie - tie;
          rank = i + 1;  // not tie, so advance to new rank
      } // ELSE tie, so don't change the rank  (note :  this is the lower of the range for median)
      key = temp[i].second;
      // add rank to sum
      if (rank_sums.count(key) == 0) rank_sums[key] = std::make_pair(1, rank);
      else {
          ++rank_sums[key].first;
          rank_sums[key].second += rank;
      }
  }
    
  // reverse walk
  rank = count;
  for (size_t i = count - 1; i > 0; --i) {
      // add rank to sum
      rank_sums[temp[i].second].second += rank;

      // rank update
      if (temp[i].first > temp[i-1].first) rank = i;  // not tie, so advance to new rank
      // ELSE tie, so don't change the rank  (note :  this is the higher of the range for median)
  }
  // need to count ties...
  rank_sums[temp[0].second].second += rank;

}

// types:  0 = greater, 1 = less, 2 = twosided (default), 3 = U
template <typename LABEL, typename OT>
void wmw(std::unordered_map<LABEL, std::pair<size_t, size_t>> const & rank_sums, double const & tie_sum, size_t const & count, OT * out,
  int const & test_type, bool const & continuity) {
  
  // tuplize in and labels.
  if (count == 0) return;
  assert((count < 20) && "ERROR: count is too small (< 20) for a normal approximation\n");

  std::vector<LABEL> sorted_labels;
  LABEL key;
  for (auto item : rank_sums) {
    key = item.first;
    if (key == std::numeric_limits<LABEL>::max()) continue;

    sorted_labels.emplace_back(key);
  }
  std::sort(sorted_labels.begin(), sorted_labels.end());
  // dealing with ties amongst multiple classes?
  // since final differential comparisons are one vs all others, and the ties are split so 0.5 is added for all.

  // compute the 1 to other U stats, storing the minimum U
  double R, U1, U2, U, prod, mu, sigma, z;
  size_t n1;
  double c1 = static_cast<double>(count + 1);
  double tie_mean = tie_sum / static_cast<double>(count * (count - 1));
  constexpr double inv12 = 1.0 / 12.0;
  size_t i = 0;
  std::pair<size_t, size_t> val;
  for (auto key : sorted_labels) {
      val = rank_sums.at(key);
      R = static_cast<double>(val.second) * 0.5;
      n1 = val.first;

      // compute U stats
      prod = static_cast<double>(n1 * (count - n1));
      U1 = R - static_cast<double>((n1 * (n1 + 1)) >> 1);
      U2 = prod - U1; 
      if (test_type == 0)  // greater
          U = U1;
      else if (test_type == 1)   // less
          U = U2;
      else   // two sided
          U = std::max(U1, U2);

      // normal approximation
      mu = prod * 0.5;
      // sigma = sqrt((prod / 12.0) * ((count + 1) - tie_sum / (count * (count - 1))));
      sigma = sqrt(prod * inv12 * (c1 - tie_mean));
      z = U - mu - (continuity ? 0.5 : 0.0);
      z /= sigma;

      // convert to p-value
      // https://stackoverflow.com/questions/2328258/cumulative-normal-distribution-function-in-c-c
      // use erfc function  This is 2-tailed.
      if (test_type == 2)  // 2 sided
          out[i] = erfc( z * M_SQRT1_2 );   // this is correct.
      else if (test_type == 1)   // less
          out[i] = 1.0 - 0.5 * erfc( -z * M_SQRT1_2 );
      else if (test_type == 0)   // greater
          out[i] = 1.0 - 0.5 * erfc( -z * M_SQRT1_2 );  // this is 1-CDF = survival function
      else
          out[i] = std::min(U1, U2);
      ++i;

  }

}


//' Fast Wilcoxon-Mann-Whitney Test for dense matrix
//'
//' This implementation uses normal approximation, which works reasonably well if sample size is large (say N>=20)
//' 
//' @rdname densewmwfast
//' @param matrix an expression matrix, COLUMN-MAJOR, each col is a feature, each row a sample
//' @param labels an integer vector, each element indicating the group to which a sample belongs.
//' @param rtype 
//' \itemize{
//' \item{0} : p(greater)
//' \item{1} : p(less)
//' \item{2} : p(twoSided)
//' \item{3} : U
//' }
//' @param continuity_correction TRUE/FALSE for continuity correction
//' @param as_dataframe TRUE/FALSE - TRUE returns a dataframe, FALSE returns a matrix
//' @param threads  number of concurrent threads.
//' @return array or dataframe.  for each gene/feature, the rows for the clusters are ordered by id.
//' @name densewmwfast
//' @export
// [[Rcpp::export]]
extern SEXP densewmwfast(SEXP matrix, SEXP labels, SEXP rtype, 
    SEXP continuity_correction, 
    SEXP as_dataframe,
    SEXP threads) {
  // Rprintf("here 1\n");

  const int type=Rf_asInteger(rtype);
  if (type > 3) Rprintf("ERROR: unsupported type: %d. Supports only greater, less, twosided, and U\n", type);
  int nthreads=Rf_asInteger(threads);
  if (nthreads < 1) nthreads = 1;
  bool continuity = Rf_asLogical(continuity_correction);
  const size_t nsamples=NROW(matrix);  // n
  const size_t nfeatures=NCOL(matrix); // m
  bool _as_dataframe = Rf_asLogical(as_dataframe);
  // Rprintf("here 2\n");
  
  // get the number of unique labels.
  int *label_ptr=INTEGER(labels);
  std::unordered_set<int> info;
  for (size_t l = 0; l < nsamples; ++l, ++label_ptr) {
    info.insert(*label_ptr);
  }
  size_t label_count = info.size();
  std::vector<int> sorted_labels(info.begin(), info.end());
  std::sort(sorted_labels.begin(), sorted_labels.end());
  label_ptr = INTEGER(labels);
  // Rprintf("here 3\n");

  int intlen = snprintf(NULL, 0, "%lu", std::numeric_limits<size_t>::max());
  char * str = reinterpret_cast<char *>(malloc(intlen + 1));
  int proc_count = 0;

  // https://stackoverflow.com/questions/23547625/returning-a-data-frame-from-c-to-r

  // GET features.
  SEXP features = Rf_getAttrib(matrix, R_NamesSymbol);
  // check if features is null.  if so, make a new one.
  // https://stackoverflow.com/questions/25172419/how-can-i-get-the-sexptype-of-an-sexp-value
  if (TYPEOF(features) == NILSXP) {
    PROTECT(features = Rf_allocVector(STRSXP, nfeatures));
    ++proc_count;

    for (size_t j = 0; j < nfeatures; ++j) {
      // create string and set in clust.
      sprintf(str, "%lu", j);
      SET_STRING_ELT(features, j, Rf_mkChar(str));
      memset(str, 0, intlen + 1);
    }

  }
  // Rprintf("here 4\n");

  // double *matColPtr; // pointer to the current column of the matrix
  SEXP res, pv, clust, genenames, names, rownames, cls;
  int ncols = 1 + (_as_dataframe ? 2 : 0);
  int col_id = 0;

  if (_as_dataframe) {
    PROTECT(res = Rf_allocVector(VECSXP, ncols));   // cluster id and gene names are repeated.
    proc_count += 1;

    PROTECT(clust = Rf_allocVector(INTSXP, label_count * nfeatures));
    PROTECT(genenames = Rf_allocVector(STRSXP, label_count * nfeatures));
    PROTECT(pv = Rf_allocVector(REALSXP, label_count * nfeatures));
    proc_count += 3;

    PROTECT(cls = Rf_allocVector(STRSXP, 1)); // class attribute
    SET_STRING_ELT(cls, 0, Rf_mkChar("data.frame"));
    Rf_classgets(res, cls);
    ++proc_count;

  } else {
    // use clust for column names.
    PROTECT(clust = Rf_allocVector(STRSXP, label_count));
    ++proc_count;
    int key;
    size_t j = 0;
    for (auto key : sorted_labels) {
      // create string and set in clust.
      sprintf(str, "%d", key);
      SET_STRING_ELT(clust, j, Rf_mkChar(str));
      memset(str, 0, intlen + 1);
      ++j;
    }
    
    PROTECT(pv = Rf_allocMatrix(REALSXP, label_count, nfeatures));
    ++proc_count;

  }
  // Rprintf("here 5\n");

  // alloc output res (list?)
  if (_as_dataframe) {
    PROTECT(rownames = Rf_allocVector(STRSXP, label_count * nfeatures));  // dataframe column names.
    ++proc_count;

    // make the clusters vector.
    int * clust_ptr = INTEGER(clust);
    SEXP * features_ptr = STRING_PTR(features);
    size_t j = 0;
    for (size_t i = 0; i < nfeatures; ++i) {
      for (auto key : sorted_labels) {
        // rotate through cluster labels for this feature.        
        *clust_ptr = key;
        ++clust_ptr;

        // same feature name for the set of cluster
        SET_STRING_ELT(genenames, j, Rf_duplicate(*features_ptr));

        sprintf(str, "%lu", j);
        SET_STRING_ELT(rownames, j, Rf_mkChar(str));
        memset(str, 0, intlen + 1);
  
        ++j;
      }
      ++features_ptr;

    }

    PROTECT(names = Rf_allocVector(STRSXP, ncols));  // dataframe column names.
    ++proc_count;

    SET_VECTOR_ELT(res, col_id, clust);  
    SET_STRING_ELT(names, col_id, Rf_mkChar("cluster"));
    ++col_id;
    SET_VECTOR_ELT(res, col_id, genenames);
    SET_STRING_ELT(names, col_id, Rf_mkChar("gene"));
    ++col_id;
    SET_VECTOR_ELT(res, col_id, pv);
    SET_STRING_ELT(names, col_id, Rf_mkChar("p_val"));   
    ++col_id;
    Rf_namesgets(res, names);  // colnames

    // set row names - NEEDED to print the dataframe!!!
    Rf_setAttrib(res, R_RowNamesSymbol, rownames);


  } else {
    // set col and row names for pv.
    // https://stackoverflow.com/questions/5709940/r-extension-in-c-setting-matrix-row-column-names
    PROTECT(names = Rf_allocVector(VECSXP, 2));
    ++proc_count;
    SET_VECTOR_ELT(names, 0, clust);  // rows = clusters
    SET_VECTOR_ELT(names, 1, features);  // columns  = features (genes)
    Rf_setAttrib(pv, R_DimNamesSymbol, names);
  }
  // Rprintf("here 6\n");

  // Rprintf("inputsize  r %lu c %lu , output r %lu c %lu \n", nsamples, nfeatures, NROW(res), NCOL(res));

  omp_set_num_threads(nthreads);
  Rprintf("THREADING: using %d threads\n", nthreads);


#pragma omp parallel num_threads(nthreads)
  {
    int tid = omp_get_thread_num();
    size_t block = nfeatures / nthreads;
    size_t rem = nfeatures - nthreads * block;
    size_t offset = tid * block + (tid > rem ? rem : tid);
    int nid = tid + 1;
    size_t end = nid * block + (nid > rem ? rem : nid);

    double *matptr = REAL(matrix) + offset * nsamples;
    double *pvptr = REAL(pv) + offset * label_count;
    std::unordered_map<int, std::pair<size_t, size_t>> rank_sums;
    double tie_sum;
    for(size_t i=offset; i < end; ++i) {
      // directly compute matrix and res pointers.
      // Rprintf("thread %d processing feature %d\n", omp_get_thread_num(), i);
      dense_wmw_summary(matptr, label_ptr, nsamples, rank_sums, tie_sum);
      wmw(rank_sums, tie_sum, nsamples, pvptr, type, continuity);
      matptr += nsamples;
      pvptr += label_count;
    }
  }
  // Rprintf("here 7\n");

  UNPROTECT(proc_count);
  free(str);
  if (_as_dataframe) return(res);
  else return(pv);
}


// =================================
