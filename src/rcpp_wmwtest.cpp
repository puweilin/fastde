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

// types:  0 = greater, 1 = less, 2 = twosided (default), 3 = U
template <typename IT, typename LABEL, typename OT>
void wmw(IT const * in, LABEL const * labels, size_t const & count, OT * out, size_t const & label_count, int const & test_type, bool const & continuity) {
  
  // tuplize in and labels.
  if (count == 0) return;
  assert((count < 20) && "ERROR: count is too small (< 20) for a normal approximation\n");
  // if (count < 100) {
  //     printf("WARNING: count is small for a normal approximation: %ld\n", count);
  // }

  // populate the sort structure.
  std::vector<std::pair<IT, LABEL>> temp;
  temp.reserve(count);
  for (size_t i = 0; i < count; ++i) {
      temp.emplace_back(in[i], labels[i]);
  }

  // sort 
  std::sort(temp.begin(), temp.end(), [](std::pair<IT, LABEL> const & a, std::pair<IT, LABEL> const & b){
      return a.first < b.first;
  });

  // assign rank and accumulate count and rank
  std::unordered_map<LABEL, std::pair<size_t, size_t>> rank_sums;  // random access of this...

  // to handle ties, we need 0.5.  To make comparison easier and not float based, double the rank value here and divide later.
  // actually, we may need to have median of the tie ranks.
  // here is an approach:  keep undoubled rank, walk once from front, and once from back.  if ties, do not change rank until value changes.  each sum /2 would give correct rank
  //    at end of tie, reset to use array index as rank.
  // any tie detection would need 1x iteration.   each tie would need another pass to get median so may break prefetching, but this may not be common?
  // should minimize writes in the big array as well.

  size_t rank = 1;
  double tie, tie_sum = 0.0;
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

    
  // dealing with ties amongst multiple classes?
  // since final differential comparisons are one vs all others, and the ties are split so 0.5 is added for all.

  // compute the 1 to other U stats, storing the minimum U
  double R, U1, U2, U, prod, mu, sigma, z;
  size_t n1;
  double c1 = static_cast<double>(count + 1);
  double cc_1 = static_cast<double>(count * (count - 1));
  constexpr double inv12 = 1.0 / 12.0;
  for (size_t i = 0; i < label_count; ++i) {
      key = i+1;
      R = static_cast<double>(rank_sums[key].second) * 0.5;
      n1 = rank_sums[key].first;

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
      sigma = sqrt(prod * inv12 * (c1 - tie_sum / cc_1));
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

  }

}


/*! \brief Wilcoxon-Mann-Whitney Test
 *
 * \param matrix: an expression matrix, each row is a feature, each column corresponds to a samples
 * \param labels: an integer vector, each element indicating the group to which a sample belongs.
 * \param rtype:
 * \parblock
 * Define f(x)=abs(log10(x))
 * 0=p(greater), 1=p(less), 2=p(twoSided), 3=U,
 * 4=f(p(greater)),5=p(less),6=f(p(twoSided)), 7=p(greater)<p(less) ? f(p(twoSided)) : -f(p(twoSided))
 * \endparblock
 * 
 * This implementation uses normal approximation, which works reasonably well if sample size is large (say N>=20)
 */
// [[Rcpp::export]]
extern SEXP wmwtest(SEXP matrix, SEXP labels, SEXP rtype, SEXP continuity_correction, SEXP threads) {
  const int type=INTEGER(rtype)[0];
  if (type > 3) Rprintf("ERROR: unsupported type: %d. Supports only greater, less, twosided, and U\n", type);
  int nthreads=INTEGER(threads)[0];
  if (nthreads < 1) nthreads = 1;
  bool continuity = LOGICAL(continuity_correction)[0];
  const size_t nsamples=NROW(matrix);  // n
  const size_t nfeatures=NCOL(matrix); // m
  
  // get the number of unique labels.
  size_t label_count = 0;
  int *label_ptr=INTEGER(labels);
  {
    std::unordered_set<int> unique_labels;
    for (size_t l = 0; l < nsamples; ++l, ++label_ptr) {
      unique_labels.insert(*label_ptr);
    }
    label_count = unique_labels.size();   // k
  }
  label_ptr = INTEGER(labels);

  // double *matColPtr; // pointer to the current column of the matrix
  SEXP res;
  // double *resPtr;  // pointer to the current column of the output 
   

  res=PROTECT(Rf_allocMatrix(REALSXP, label_count, nfeatures));
  // Rprintf("inputsize  r %lu c %lu , output r %lu c %lu \n", nsamples, nfeatures, NROW(res), NCOL(res));

  omp_set_num_threads(nthreads);

#pragma omp parallel for num_threads(nthreads) 
  for(size_t i=0; i < nfeatures; ++i) {
    // directly compute matrix and res pointers.
    Rprintf("thread %d processing feature %d\n", omp_get_thread_num(), i);
    wmw(REAL(matrix) + i * nsamples, label_ptr, nsamples, REAL(res) + i * label_count, label_count, type, continuity);
  }

  UNPROTECT(1);
  return(res);
}


// =================================
