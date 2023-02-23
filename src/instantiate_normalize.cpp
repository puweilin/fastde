
// https://stackoverflow.com/questions/44795190/how-can-i-return-a-list-of-matrices-from-rcpp-to-r

/*
 *  similar structure as wilcox.hpp
 *
 *  Created on: Jun 8, 2021
 *  Author: Tony Pan
 *  Affiliation: Institute for Data Engineering and Science
 *  			Georgia Institute of Technology, Atlanta, GA 30332
 */

#include "fastde/normalize.tpp"

#include <vector>
#include <unordered_map>
#include <utility>

#include "cpp11/doubles.hpp"
#include "cpp11/integers.hpp"



template void csc_log_normalize_vec(cpp11::doubles const & x, cpp11::integers const & p, size_t const & cols, double const & scale_factor, cpp11::writable::doubles & out, int const & threads);

template void csc_log_normalize_vec(cpp11::doubles const & x, cpp11::doubles const & p, size_t const & cols, double const & scale_factor, cpp11::writable::doubles & out, int const & threads);



template void csc_clr_cols_vec(cpp11::doubles const & x, cpp11::integers const & i, cpp11::integers const & p, size_t const & rows, size_t const & cols, cpp11::writable::doubles & out, int const & threads);

template void csc_clr_cols_vec(cpp11::doubles const & x, cpp11::integers const & i, cpp11::doubles const & p, size_t const & rows, size_t const & cols, cpp11::writable::doubles & out, int const & threads);

template void csc_clr_rows_vec(cpp11::doubles const & x, cpp11::integers const & i, cpp11::integers const & p, size_t const & rows, size_t const & cols, cpp11::writable::doubles & out, int const & threads);

template void csc_clr_rows_vec(cpp11::doubles const & x, cpp11::integers const & i, cpp11::doubles const & p, size_t const & rows, size_t const & cols, cpp11::writable::doubles & out, int const & threads);


template void csc_relative_count_vec(cpp11::doubles const & x, cpp11::integers const & p, size_t const & cols, double const & scale_factor, cpp11::writable::doubles & out, int const & threads);

template void csc_relative_count_vec(cpp11::doubles const & x, cpp11::doubles const & p, size_t const & cols, double const & scale_factor, cpp11::writable::doubles & out, int const & threads);



