#include "fastde/normalize.hpp"

#include <chrono>

#include <cpp11/sexp.hpp>
#include <cpp11/matrix.hpp>
#include <cpp11/strings.hpp>
#include <cpp11/integers.hpp>
#include <cpp11/doubles.hpp>

#include "utils_data.hpp"
#include "fastde/benchmark_utils.hpp"


[[cpp11::register]]
extern cpp11::writable::list cpp11_sp_normalize(cpp11::doubles const & x,
    cpp11::integers const & i, cpp11::integers const & p, int const & nrow, int const & ncol,
    double const & scale_factor,
    int const & method, int const & threads) {

    size_t nz = p[ncol];

    cpp11::writable::doubles xv(nz);

    if (method == 0) {
      // log normal
      csc_log_normalize_vec(x, p, ncol, scale_factor, xv, threads);
    } else if (method == 1) {
      // clr
      csc_clr_cols_vec(x, i, p, nrow, ncol, xv, threads);
    } else if (method == 2) {
      // relative count.
      csc_relative_count_vec(x, p, ncol, scale_factor, xv, threads);
    }

    cpp11::writable::list out;
    // create output.
    out.push_back(xv);
    cpp11::writable::integers iv(i);
    out.push_back(iv);
    cpp11::writable::integers pv(p);
    out.push_back(pv);
    out.push_back(cpp11::as_sexp(ncol));
    out.push_back(cpp11::as_sexp(nrow));

    return out;

}



[[cpp11::register]]
extern cpp11::writable::list cpp11_sp64_normalize(cpp11::doubles const & x,
    cpp11::integers const & i, cpp11::doubles const & p, int const & nrow, int const & ncol, double const & scale_factor,
    int const & method, int const & threads) {
      
    size_t nz = p[ncol];

    cpp11::writable::doubles xv(nz);

    if (method == 0) {
      // log normal
      csc_log_normalize_vec(x, p, ncol, scale_factor, xv, threads);
    } else if (method == 1) {
      // clr
      csc_clr_cols_vec(x, i, p, nrow, ncol, xv, threads);
    } else if (method == 2) {
      // relative count.
      csc_relative_count_vec(x, p, ncol, scale_factor, xv, threads);
    }

    cpp11::writable::list out;
    // create output.
    out.push_back(xv);
    cpp11::writable::integers iv(i);
    out.push_back(iv);
    cpp11::writable::doubles pv(p);
    out.push_back(pv);
    out.push_back(cpp11::as_sexp(ncol));
    out.push_back(cpp11::as_sexp(nrow));

    return out;


}