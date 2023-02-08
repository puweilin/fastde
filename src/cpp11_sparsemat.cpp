#include "utils_sparsemat.hpp"
#include "cpp11/doubles.hpp"
#include "cpp11/integers.hpp"
#include "cpp11/list.hpp"
#include "cpp11/matrix.hpp"

#include "cpp11/r_vector.hpp"
#include "cpp11/list_of.hpp"
#include <R.h>
#include <vector>

#include "fastde/sparsemat.hpp"

// adapters for cpp11 r_vector to c++ implementation

[[cpp11::register]]
extern cpp11::writable::list cpp11_sp_transpose(cpp11::doubles const & x,
    cpp11::integers const & i, cpp11::integers const & p, int const & nrow, int const & ncol, int const & threads) {

    // size_t nz = p[ncol];

    // cpp11::writable::doubles xv(nz);
    // cpp11::writable::integers iv(nz);
    // cpp11::writable::integers pv(nrow + 1);

    // csc_transpose_csc(x.cbegin(), i.cbegin(), p.cbegin(), nrow, ncol, xv.begin(), iv.begin(), pv.begin(), threads);

    // cpp11::writable::list out;
    // // create output.
    // out.push_back(xv);
    // out.push_back(iv);
    // out.push_back(pv);
    // out.push_back(cpp11::as_sexp(ncol));
    // out.push_back(cpp11::as_sexp(nrow));
    
    // return out;

    if (threads == 1) return _sp_transpose(x, i, p, nrow, ncol, threads);
    else return _sp_transpose_par(x, i, p, nrow, ncol, threads);
}


[[cpp11::register]]
extern cpp11::writable::list cpp11_sp64_transpose(cpp11::doubles const & x,
    cpp11::integers const & i, cpp11::doubles const & p, int const & nrow, int const & ncol, int const & threads) {

    // size_t nz = p[ncol];

    // cpp11::writable::doubles xv(nz);
    // cpp11::writable::integers iv(nz);
    // cpp11::writable::doubles pv(nrow + 1);

    // csc_transpose_csc(x.cbegin(), i.cbegin(), p.cbegin(), nrow, ncol, xv.begin(), iv.begin(), pv.begin(), threads);

    // cpp11::writable::list out;
    // // create output.
    // out.push_back(xv);
    // out.push_back(iv);
    // out.push_back(pv);
    // out.push_back(cpp11::as_sexp(ncol));
    // out.push_back(cpp11::as_sexp(nrow));
    
    // return out;

    if (threads == 1) return _sp_transpose(x, i, p, nrow, ncol, threads);
    else return _sp_transpose_par(x, i, p, nrow, ncol, threads);
}



[[cpp11::register]]
extern cpp11::writable::doubles_matrix<cpp11::by_column> cpp11_sp_to_dense(cpp11::doubles const & x,
    cpp11::integers const & i, cpp11::integers const & p, int const & nrow, int const & ncol, int const & threads) {

    // std::vector<double> vec(nrow * ncol, 0);
    // csc_to_dense_c(x.cbegin(), i.cbegin(), p.cbegin(), nrow, ncol, vec.begin(), threads);

    // cpp11::writable::doubles_matrix<cpp11::by_column> out(nrow, ncol);

    // auto slice_start = out.begin();
    // auto slice_end = out.end();
    // auto viter = vec.cbegin();
    // for (; slice_start != slice_end; ++slice_start) {

    //     auto el_start = (*slice_start).begin();
    //     std::copy(viter, viter + nrow, el_start);
    // }

    // return out;
    return _sp_to_dense<cpp11::writable::doubles_matrix<cpp11::by_column>>(x, i, p, nrow, ncol, threads);
}


[[cpp11::register]]
extern cpp11::writable::doubles_matrix<cpp11::by_column> cpp11_sp64_to_dense(cpp11::doubles const & x,
    cpp11::integers const & i, cpp11::doubles const & p, int const & nrow, int const & ncol, int const & threads) {

    // std::vector<double> vec(nrow * ncol, 0);
    // csc_to_dense_c(x.cbegin(), i.cbegin(), p.cbegin(), nrow, ncol, vec.begin(), threads);

    // cpp11::writable::doubles_matrix<cpp11::by_column> out(nrow, ncol);

    // auto slice_start = out.begin();
    // auto slice_end = out.end();
    // auto viter = vec.cbegin();
    // for (; slice_start != slice_end; ++slice_start) {

    //     auto el_start = (*slice_start).begin();
    //     std::copy(viter, viter + nrow, el_start);
    // }

    // return out;
    return _sp_to_dense<cpp11::writable::doubles_matrix<cpp11::by_column>>(x, i, p, nrow, ncol, threads);
}



[[cpp11::register]]
extern cpp11::writable::doubles_matrix<cpp11::by_column> cpp11_sp_to_dense_transposed(
    cpp11::doubles const & x,
    cpp11::integers const & i, 
    cpp11::integers const & p, 
    int const & nrow, int const & ncol,
    int const & threads
) {
    // std::vector<double> vec(nrow * ncol, 0);
    // csc_to_dense_transposed_c(x.cbegin(), i.cbegin(), p.cbegin(), nrow, ncol, vec.begin(), threads);

    // cpp11::writable::doubles_matrix<cpp11::by_column> out(ncol, nrow);

    // auto slice_start = out.begin();
    // auto slice_end = out.end();
    // auto viter = vec.cbegin();
    // for (; slice_start != slice_end; ++slice_start) {

    //     auto el_start = (*slice_start).begin();
    //     std::copy(viter, viter + ncol, el_start);
    // }

    // return out;
    return _sp_to_dense_transposed<cpp11::writable::doubles_matrix<cpp11::by_column>>(x, i, p, nrow, ncol, threads);
}

[[cpp11::register]]
extern cpp11::writable::doubles_matrix<cpp11::by_column> cpp11_sp64_to_dense_transposed(
    cpp11::doubles const & x,
    cpp11::integers const & i, 
    cpp11::doubles const & p, 
    int const & nrow, int const & ncol,
    int const & threads
) {
    // std::vector<double> vec(nrow * ncol, 0);
    // csc_to_dense_transposed_c(x.cbegin(), i.cbegin(), p.cbegin(), nrow, ncol, vec.begin(), threads);

    // cpp11::writable::doubles_matrix<cpp11::by_column> out(ncol, nrow);

    // auto slice_start = out.begin();
    // auto slice_end = out.end();
    // auto viter = vec.cbegin();
    // for (; slice_start != slice_end; ++slice_start) {

    //     auto el_start = (*slice_start).begin();
    //     std::copy(viter, viter + ncol, el_start);
    // }

    // return out;
    return _sp_to_dense_transposed<cpp11::writable::doubles_matrix<cpp11::by_column>>(x, i, p, nrow, ncol, threads);
}


// this is a test to see how to use SEXP directly.
// namespace cpp11 {

// template <typename T>
// typename std::enable_if<std::is_same<T, std::vector<cpp11::r_vector<double>>>::value, T>::type as_cpp(SEXP from) {
//     T out;
//     if (TYPEOF(from) == VECSXP) {
//         R_xlen_t nlists = XLENGTH(from);
//         out.reserve(nlists);

//         for (R_xlen_t i = 0; i < nlists; ++i) {
//             // get list items in "from", convert to r_vector, and insert into output
//             out.emplace_back(cpp11::r_vector<double>(VECTOR_ELT(from, i)));
//         }
//     }
//     return out;
// }

// template <typename T>
// typename std::enable_if<std::is_same<T, std::vector<cpp11::r_vector<int>>>::value, T>::type as_cpp(SEXP from) {
//     T out;
//     if (TYPEOF(from) == VECSXP) {
//         R_xlen_t nlists = XLENGTH(from);
//         out.reserve(nlists);

//         for (R_xlen_t i = 0; i < nlists; ++i) {
//             // get list items in "from", convert to r_vector, and insert into output
//             out.emplace_back(cpp11::r_vector<int>(VECTOR_ELT(from, i)));
//         }
//     }
//     return out;
// }

// }


// this is much easier way to map.
// question - do we want to merge THEN compute using existing APIs?   would be easier....

// THIS IS A GOOD COMPILE TEST.
// [[cpp11::register]]
// extern int merge_sparse_mats(
//     cpp11::list_of<cpp11::doubles> const & xvecs,
//     cpp11::list_of<cpp11::integers> const & ivecs,
//     cpp11::list_of<cpp11::doubles> const & pvecs,
//     cpp11::integers const & nrows,
//     cpp11::integers const & ncols
// ) {
//     // have to do this because list_of.hpp's operator[] are not set up to work on const, and would throw an error.
//     // return cpp11::as_integers(ivecs.at(0))[0];

//     // below requries const operator[]
//     return ivecs[0][0];
// }

// "list_of" does note react well to "const &" - reports SEXPREC declared but not defined.
// this is because list is an r_vector<SEXP> and SEXP is an opaque pointer to SEXPREC, so it's not know that it is doubles array or integers array.
// data is in column major
[[cpp11::register]]
extern cpp11::writable::list cpp11_sp_rbind(
    cpp11::list_of<cpp11::doubles> const & xvecs,
    cpp11::list_of<cpp11::integers> const & ivecs,
    cpp11::list_of<cpp11::integers> const & pvecs,
    cpp11::integers const & nrows,
    cpp11::integers const & ncols,
    int const & threads,
    int const & method = 1
    ) {

    cpp11::writable::list out;

    int n_vecs = nrows.size();
    if (n_vecs == 0) return out;

    // check to see if all inputs have the same number of rows.
    int nrow = 0;
    int ncol = ncols[0];
    size_t nz = 0;
    bool equal_cols = true;

    if (method == 0) {
        // this is just for testing.
        std::vector<std::vector<double>> xsv(n_vecs);
        std::vector<std::vector<int>> isv(n_vecs);
        std::vector<std::vector<int>> psv(n_vecs);
        std::vector<std::vector<double>::const_iterator> xs;
        std::vector<std::vector<int>::const_iterator> is;
        std::vector<std::vector<int>::const_iterator> ps;
        xs.reserve(n_vecs);
        is.reserve(n_vecs);
        ps.reserve(n_vecs);
        std::vector<int> nr(n_vecs);
        std::vector<int> nc(n_vecs);

        size_t test_nz = 0;
        int test_col = 0;
        for (int i = 0; i < n_vecs; ++i) {
            test_col = ncols[i];
            test_nz = pvecs[i][test_col];

            equal_cols &= (ncol == test_col); // bad input.
            nrow += nrows[i];
            nz += test_nz;

            xsv[i].resize(test_nz);         std::copy(xvecs[i].cbegin(), xvecs[i].cend(), xsv[i].begin());
            isv[i].resize(test_nz);         std::copy(ivecs[i].cbegin(), ivecs[i].cend(), isv[i].begin());
            psv[i].resize(test_col + 1);    std::copy(pvecs[i].cbegin(), pvecs[i].cend(), psv[i].begin());

            xs.emplace_back(xsv[i].cbegin());
            is.emplace_back(isv[i].cbegin());
            ps.emplace_back(psv[i].cbegin());
            
            nr[i] = nrows[i];
            nc[i] = ncols[i];
        }
        if (! equal_cols) return out;

        // Rprintf("allocating %d non-zeros, %d columns, and %d rows!\n", nz, ncol, nrow);

        cpp11::writable::doubles xv(nz);
        cpp11::writable::integers iv(nz);

        if (nz >= __INT_MAX__) {
            cpp11::writable::doubles pv(ncol + 1);
            csc_rbind(xs, is, ps, nr, nc, xv.begin(), iv.begin(), pv.begin(), threads);
            // create output.
            out.push_back(xv);
            out.push_back(iv);
            out.push_back(pv);
            out.push_back(cpp11::as_sexp(nrow));
            out.push_back(cpp11::as_sexp(ncol));
            
        } else {
            cpp11::writable::integers pv(ncol + 1);
            csc_rbind(xs, is, ps, nr, nc, xv.begin(), iv.begin(), pv.begin(), threads);
            // create output.
            out.push_back(xv);
            out.push_back(iv);
            out.push_back(pv);
            out.push_back(cpp11::as_sexp(nrow));
            out.push_back(cpp11::as_sexp(ncol));

        }
        return out;
    } else if (method == 1) {
// -----------------

        for (int i = 0; i < n_vecs; ++i) {
            equal_cols &= (ncol == ncols[i]); // bad input.
            nrow += nrows[i];
            nz += pvecs[i][ncols[i]];
        }
        if (! equal_cols) return out;

        cpp11::writable::doubles xv(nz);
        cpp11::writable::integers iv(nz);

        if (nz >= __INT_MAX__) {
            cpp11::writable::doubles pv(ncol + 1);
            csc_rbind_vec(xvecs, ivecs, pvecs, nrows, ncols,  xv, iv, pv, threads);
            // create output.
            out.push_back(xv);
            out.push_back(iv);
            out.push_back(pv);
            out.push_back(cpp11::as_sexp(nrow));
            out.push_back(cpp11::as_sexp(ncol));
            
        } else {
            cpp11::writable::integers pv(ncol + 1);
            csc_rbind_vec(xvecs, ivecs, pvecs, nrows, ncols,  xv, iv, pv, threads);
            // create output.
            out.push_back(xv);
            out.push_back(iv);
            out.push_back(pv);
            out.push_back(cpp11::as_sexp(nrow));
            out.push_back(cpp11::as_sexp(ncol));

        }
        return out;
    } else if (method == 2) {
// -----------------
        // copy input to cpp
        // std::vector<std::vector<double>> xsv(n_vecs);
        // std::vector<std::vector<int>> isv(n_vecs);
        // std::vector<std::vector<int>> psv(n_vecs);

        size_t test_nz = 0;
        int test_col = 0;
        for (int i = 0; i < n_vecs; ++i) {
            test_col = ncols[i];
            test_nz = pvecs[i][test_col];

            equal_cols &= (ncol == test_col); // bad input.
            nrow += nrows[i];
            nz += test_nz;

            // xsv[i].resize(test_nz);         std::copy(xvecs[i].cbegin(), xvecs[i].cend(), xsv[i].begin());
            // isv[i].resize(test_nz);         std::copy(ivecs[i].cbegin(), ivecs[i].cend(), isv[i].begin());
            // psv[i].resize(test_col + 1);    std::copy(pvecs[i].cbegin(), pvecs[i].cend(), psv[i].begin());

        }
        if (! equal_cols) return out;


        if (nz >= __INT_MAX__) {

            std::vector<double> oxv(nz);
            std::vector<int> oiv(nz);
            std::vector<double> opv(ncol + 1);

            csc_rbind_vec(xvecs, ivecs, pvecs, nrows, ncols,  oxv, oiv, opv, threads);

            cpp11::writable::doubles xv(oxv.cbegin(), oxv.cend());
            cpp11::writable::integers iv(oiv.cbegin(), oiv.cend());
            cpp11::writable::doubles pv(opv.cbegin(), opv.cend());
            
            // create output.
            out.push_back(xv);
            out.push_back(iv);
            out.push_back(pv);
            out.push_back(cpp11::as_sexp(nrow));
            out.push_back(cpp11::as_sexp(ncol));
            
        } else {

            std::vector<double> oxv(nz);
            std::vector<int> oiv(nz);
            std::vector<int> opv(ncol + 1);

            csc_rbind_vec(xvecs, ivecs, pvecs, nrows, ncols, oxv, oiv, opv, threads);

            // check
            Rprintf("i0: %d, i1 %d, i2 %d\n", oiv[0], oiv[1], oiv[2]);
            Rprintf("p0: %d, p1: %d, p2: %d\n", opv[0], opv[1], opv[2]);

            cpp11::writable::doubles xv(oxv.cbegin(), oxv.cend());
            cpp11::writable::integers iv(oiv.cbegin(), oiv.cend());
            cpp11::writable::integers pv(opv.cbegin(), opv.cend());

            // check
            Rprintf("AFTER i0: %d, i1 %d, i2 %d\n", iv[0], iv[1], iv[2]);
            Rprintf("AFTER p0: %d, p1: %d, p2: %d\n", pv[0], pv[1], pv[2]);

            // create output.
            out.push_back(xv);
            out.push_back(iv);
            out.push_back(pv);
            out.push_back(cpp11::as_sexp(nrow));
            out.push_back(cpp11::as_sexp(ncol));

        }
        return out;
    } else {
        return _sp_rbind(xvecs, ivecs, pvecs, nrows, ncols, threads);
    }
}


[[cpp11::register]]
extern cpp11::writable::list cpp11_sp64_rbind(
    cpp11::list_of<cpp11::doubles> const & xvecs,
    cpp11::list_of<cpp11::integers> const & ivecs,
    cpp11::list_of<cpp11::doubles> const & pvecs,
    cpp11::integers const & nrows,
    cpp11::integers const & ncols,
    int const & threads,
    int const & method = 1
) {

    cpp11::writable::list out;

    int n_vecs = nrows.size();
    if (n_vecs == 0) return out;

    // check to see if all inputs have the same number of rows.
    int nrow = 0;
    int ncol = ncols[0];
    size_t nz = 0;
    bool equal_cols = true;

//------------------
    if (method == 0) {
        // this is just for testing.
        std::vector<std::vector<double>> xsv(n_vecs);
        std::vector<std::vector<int>> isv(n_vecs);
        std::vector<std::vector<double>> psv(n_vecs);
        std::vector<std::vector<double>::const_iterator> xs;
        std::vector<std::vector<int>::const_iterator> is;
        std::vector<std::vector<double>::const_iterator> ps;
        xs.reserve(n_vecs);
        is.reserve(n_vecs);
        ps.reserve(n_vecs);
        std::vector<int> nr(n_vecs);
        std::vector<int> nc(n_vecs);

        size_t test_nz = 0;
        int test_col = 0;
        for (int i = 0; i < n_vecs; ++i) {
            test_col = ncols[i];
            test_nz = pvecs[i][test_col];

            equal_cols &= (ncol == test_col); // bad input.
            nrow += nrows[i];
            nz += test_nz;

            xsv[i].resize(test_nz);         std::copy(xvecs[i].cbegin(), xvecs[i].cend(), xsv[i].begin());
            isv[i].resize(test_nz);         std::copy(ivecs[i].cbegin(), ivecs[i].cend(), isv[i].begin());
            psv[i].resize(test_col + 1);    std::copy(pvecs[i].cbegin(), pvecs[i].cend(), psv[i].begin());

            xs.emplace_back(xsv[i].cbegin());
            is.emplace_back(isv[i].cbegin());
            ps.emplace_back(psv[i].cbegin());
            
            nr[i] = nrows[i];
            nc[i] = ncols[i];
        }
        if (! equal_cols) return out;

        cpp11::writable::doubles xv(nz);
        cpp11::writable::integers iv(nz);
        cpp11::writable::doubles pv(ncol + 1);


        csc_rbind(xs, is, ps, nr, nc, xv.begin(), iv.begin(), pv.begin(), threads);

        // create output.
        out.push_back(xv);
        out.push_back(iv);
        out.push_back(pv);
        out.push_back(cpp11::as_sexp(nrow));
        out.push_back(cpp11::as_sexp(ncol));
        
        return out;
    } else if (method == 1) {
//------------------

        for (int i = 0; i < n_vecs; ++i) {
            equal_cols &= (ncol == ncols[i]); // bad input.
            nrow += nrows[i];
            nz += pvecs[i][ncols[i]];

        }
        if (! equal_cols) return out;

        cpp11::writable::doubles xv(nz);
        cpp11::writable::integers iv(nz);

        cpp11::writable::doubles pv(ncol + 1);
        csc_rbind_vec(xvecs, ivecs, pvecs, nrows, ncols, xv, iv, pv, threads);
        // create output.
        out.push_back(xv);
        out.push_back(iv);
        out.push_back(pv);
        out.push_back(cpp11::as_sexp(nrow));
        out.push_back(cpp11::as_sexp(ncol));
            
        return out;
    } else {
        return _sp_rbind(xvecs, ivecs, pvecs, nrows, ncols, threads);
    }
}

[[cpp11::register]]
extern cpp11::writable::list cpp11_sp_cbind(
    cpp11::list_of<cpp11::doubles> const & xvecs,
    cpp11::list_of<cpp11::integers> const & ivecs,
    cpp11::list_of<cpp11::integers> const & pvecs,
    cpp11::integers const & nrows,
    cpp11::integers const & ncols,
    int const & threads,
    int const & method = 1
) {

    cpp11::writable::list out;

    int n_vecs = nrows.size();
    if (n_vecs == 0) return out;

    // check to see if all inputs have the same number of rows.
    int nrow = nrows[0];
    int ncol = 0;
    size_t nz = 0;
    bool equal_rows = true;

// --------------------
    if (method == 0) {
        std::vector<std::vector<double>> xsv(n_vecs);
        std::vector<std::vector<int>> isv(n_vecs);
        std::vector<std::vector<int>> psv(n_vecs);
        std::vector<std::vector<double>::const_iterator> xs;
        std::vector<std::vector<int>::const_iterator> is;
        std::vector<std::vector<int>::const_iterator> ps;
        xs.reserve(n_vecs);
        is.reserve(n_vecs);
        ps.reserve(n_vecs);
        std::vector<int> nr(n_vecs);
        std::vector<int> nc(n_vecs);

        size_t test_nz = 0;
        int test_col = 0;
        for (int i = 0; i < n_vecs; ++i) {
            test_col = ncols[i];
            test_nz = pvecs[i][test_col];

            equal_rows &= (nrow == nrows[i]); // bad input.
            ncol += test_col;
            nz += test_nz;

            xsv[i].resize(test_nz);         std::copy(xvecs[i].cbegin(), xvecs[i].cend(), xsv[i].begin());
            isv[i].resize(test_nz);         std::copy(ivecs[i].cbegin(), ivecs[i].cend(), isv[i].begin());
            psv[i].resize(test_col + 1);    std::copy(pvecs[i].cbegin(), pvecs[i].cend(), psv[i].begin());

            xs.emplace_back(xsv[i].cbegin());
            is.emplace_back(isv[i].cbegin());
            ps.emplace_back(psv[i].cbegin());

            nr[i] = nrows[i];
            nc[i] = ncols[i];
        }
        if (! equal_rows) return out;

        cpp11::writable::doubles xv(nz);
        cpp11::writable::integers iv(nz);

        if (nz >= __INT_MAX__) {
            cpp11::writable::doubles pv(ncol + 1);
            csc_cbind(xs, is, ps, nr, nc, xv.begin(), iv.begin(), pv.begin(), threads);
            // create output.
            out.push_back(xv);
            out.push_back(iv);
            out.push_back(pv);
            out.push_back(cpp11::as_sexp(nrow));
            out.push_back(cpp11::as_sexp(ncol));
            
        } else {
            cpp11::writable::integers pv(ncol + 1);
            csc_cbind(xs, is, ps, nr, nc, xv.begin(), iv.begin(), pv.begin(), threads);
            // create output.
            out.push_back(xv);
            out.push_back(iv);
            out.push_back(pv);
            out.push_back(cpp11::as_sexp(nrow));
            out.push_back(cpp11::as_sexp(ncol));

        }

        return out;
    } else if (method == 1) {
//--------------
        for (int i = 0; i < n_vecs; ++i) {
            equal_rows &= (nrow == nrows[i]); // bad input.
            ncol += ncols[i];
            nz += pvecs[i][ncols[i]];

        }
        if (! equal_rows) return out;

        cpp11::writable::doubles xv(nz);
        cpp11::writable::integers iv(nz);

        if (nz >= __INT_MAX__) {
            cpp11::writable::doubles pv(ncol + 1);
            csc_cbind_vec(xvecs, ivecs, pvecs, nrows, ncols, xv, iv, pv, threads);
            // create output.
            out.push_back(xv);
            out.push_back(iv);
            out.push_back(pv);
            out.push_back(cpp11::as_sexp(nrow));
            out.push_back(cpp11::as_sexp(ncol));
            
        } else {
            cpp11::writable::integers pv(ncol + 1);
            csc_cbind_vec(xvecs, ivecs, pvecs, nrows, ncols, xv, iv, pv, threads);
            // create output.
            out.push_back(xv);
            out.push_back(iv);
            out.push_back(pv);
            out.push_back(cpp11::as_sexp(nrow));
            out.push_back(cpp11::as_sexp(ncol));

        }
        return out;
    } else {
        return _sp_cbind(xvecs, ivecs, pvecs, nrows, ncols, threads);
    }
}

[[cpp11::register]]
extern cpp11::writable::list cpp11_sp64_cbind(
    cpp11::list_of<cpp11::doubles> const & xvecs,
    cpp11::list_of<cpp11::integers> const & ivecs,
    cpp11::list_of<cpp11::doubles> const & pvecs,
    cpp11::integers const & nrows,
    cpp11::integers const & ncols,
    int const & threads,
    int const & method = 1
) {
    cpp11::writable::list out;

    int n_vecs = nrows.size();
    if (n_vecs == 0) return out;

    // check to see if all inputs have the same number of rows.
    int nrow = nrows[0];
    int ncol = 0;
    size_t nz = 0;
    bool equal_rows = true;

    // ---------------------
    if (method == 0) {
        std::vector<std::vector<double>> xsv(n_vecs);
        std::vector<std::vector<int>> isv(n_vecs);
        std::vector<std::vector<double>> psv(n_vecs);
        std::vector<std::vector<double>::const_iterator> xs;
        std::vector<std::vector<int>::const_iterator> is;
        std::vector<std::vector<double>::const_iterator> ps;
        xs.reserve(n_vecs);
        is.reserve(n_vecs);
        ps.reserve(n_vecs);
        std::vector<int> nr(n_vecs);
        std::vector<int> nc(n_vecs);

        size_t test_nz = 0;
        int test_col = 0;
        for (int i = 0; i < n_vecs; ++i) {
            test_col = ncols[i];
            test_nz = pvecs[i][test_col];

            equal_rows &= (nrow == nrows[i]); // bad input.
            ncol += test_col;
            nz += test_nz;

            xsv[i].resize(test_nz);         std::copy(xvecs[i].cbegin(), xvecs[i].cend(), xsv[i].begin());
            isv[i].resize(test_nz);         std::copy(ivecs[i].cbegin(), ivecs[i].cend(), isv[i].begin());
            psv[i].resize(test_col + 1);    std::copy(pvecs[i].cbegin(), pvecs[i].cend(), psv[i].begin());

            xs.emplace_back(xsv[i].cbegin());
            is.emplace_back(isv[i].cbegin());
            ps.emplace_back(psv[i].cbegin());

            nr[i] = nrows[i];
            nc[i] = ncols[i];
        }
        if (! equal_rows) return out;

        cpp11::writable::doubles xv(nz);
        cpp11::writable::integers iv(nz);
        cpp11::writable::doubles pv(ncol + 1);


        csc_cbind(xs, is, ps, nr, nc, xv.begin(), iv.begin(), pv.begin(), threads);

        // create output.
        out.push_back(xv);
        out.push_back(iv);
        out.push_back(pv);
        out.push_back(cpp11::as_sexp(nrow));
        out.push_back(cpp11::as_sexp(ncol));

        return out;
    } else if (method == 1) {

        for (int i = 0; i < n_vecs; ++i) {
            equal_rows &= (nrow == nrows[i]); // bad input.
            ncol += ncols[i];
            nz += pvecs[i][ncols[i]];
        }
        if (! equal_rows) return out;

        cpp11::writable::doubles xv(nz);
        cpp11::writable::integers iv(nz);

        cpp11::writable::doubles pv(ncol + 1);
        csc_cbind_vec(xvecs, ivecs, pvecs, nrows, ncols, xv, iv, pv, threads);
        // create output.
        out.push_back(xv);
        out.push_back(iv);
        out.push_back(pv);
        out.push_back(cpp11::as_sexp(nrow));
        out.push_back(cpp11::as_sexp(ncol));
    
        return out;
    } else {
        return _sp_cbind(xvecs, ivecs, pvecs, nrows, ncols, threads);
    }
}




[[cpp11::register]]
extern cpp11::writable::doubles cpp11_sp_colSums(
    cpp11::doubles const & x,
    cpp11::integers const & p,
    int const & threads,
    int const & method = 1
) {
    if (method == 0) {
        cpp11::writable::doubles out(p.size() - 1);
        csc_colsums_iter(x.cbegin(), p.cbegin(), static_cast<int>(p.size() - 1), out.begin(), threads);
        return out;
    } else if (method == 1) {
        cpp11::writable::doubles out(static_cast<R_xlen_t>(p.size() - 1));
        csc_colsums_vec(x, p, static_cast<int>(p.size() - 1), out, threads);
        return out;
    } else {
        return _sp_colsums(x, p, static_cast<int>(p.size() - 1), threads );
    }
}

[[cpp11::register]]
extern cpp11::writable::doubles cpp11_sp64_colSums(
    cpp11::doubles const & x,
    cpp11::doubles const & p,
    int const & threads,
    int const & method = 1
) {
    if (method == 0) {
        cpp11::writable::doubles out(p.size() - 1);
        csc_colsums_iter(x.cbegin(), p.cbegin(), static_cast<int>(p.size() - 1), out.begin(), threads);
        return out;
    } else if (method == 1) {
        cpp11::writable::doubles out(static_cast<R_xlen_t>(p.size() - 1));
        csc_colsums_vec(x, p, static_cast<int>(p.size() - 1), out, threads);
        return out;
    } else {
        return _sp_colsums(x, p, static_cast<int>(p.size() - 1), threads );
    }
}


[[cpp11::register]]
extern cpp11::writable::doubles cpp11_sp_rowSums(
    cpp11::doubles const & x,
    cpp11::integers const & i,
    int const & nrow,
    int const & threads,
    int const & method = 1
) {
    if (method == 0) {
        cpp11::writable::doubles out(nrow);
        csc_rowsums_iter(x.cbegin(), i.cbegin(), nrow, static_cast<size_t>(x.size()), out.begin(), threads);
        return out;
    } else if (method == 1) {
        cpp11::writable::doubles out(nrow);
        csc_rowsums_vec(x, i, nrow, out, threads);
        return out;
    } else {
        return _sp_rowsums(x, i, nrow, static_cast<size_t>(x.size()), threads );
    }
}
