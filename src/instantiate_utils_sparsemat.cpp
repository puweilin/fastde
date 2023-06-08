#include "utils_sparsemat.tpp"
#include "cpp11/doubles.hpp"
#include "cpp11/integers.hpp"

#include "fastde/sparsemat.tpp"


// NOTe:  there is no formal definition of sparse matrix.
// input is column major, so i has the row ids, and p is per column.
// template cpp11::writable::list _sp_transpose_sort(
//     cpp11::doubles const & x, 
//     cpp11::integers const & i, 
//     cpp11::integers const & p, int const & nrow, int const & ncol, int const & threads);


// NOTe:  there is no formal definition of sparse matrix.
// input is column major, so i has the row ids, and p is per column.
// template cpp11::writable::list _sp_transpose_sort(
//     cpp11::doubles const & x, 
//     cpp11::integers const & i, 
//     cpp11::doubles const & p, int const & nrow, int const & ncol, int const & threads);


// NOTe:  there is no formal definition of sparse matrix.
// input is column major, so i has the row ids, and p is per column.
template cpp11::writable::list _sp_transpose_par(
    cpp11::doubles const & x, 
    cpp11::integers const & i, 
    cpp11::integers const & p, 
    int const & nrow, int const & ncol, int const & threads);

// NOTe:  there is no formal definition of sparse matrix.
// input is column major, so i has the row ids, and p is per column.
template cpp11::writable::list _sp_transpose_par(
    cpp11::doubles const & x, 
    cpp11::integers const & i, 
    cpp11::doubles const & p, 
    int const & nrow, int const & ncol, int const & threads);



// NOTe:  there is no formal definition of sparse matrix.
// input is column major, so i has the row ids, and p is per column.
template cpp11::writable::list _sp_transpose(
    cpp11::doubles const & x, 
    cpp11::integers const & i, 
    cpp11::integers const & p, int const & nrow, int const & ncol, int const & threads) ;

// NOTe:  there is no formal definition of sparse matrix.
// input is column major, so i has the row ids, and p is per column.
template cpp11::writable::list _sp_transpose(
    cpp11::doubles const & x, 
    cpp11::integers const & i, 
    cpp11::doubles const & p, int const & nrow, int const & ncol, int const & threads) ;

// NOTe:  there is no formal definition of sparse matrix.
// input is column major, so i has the row ids, and p is per column.
// direct to stl vectors
template void _sp_transpose_par(
    cpp11::doubles const & x, 
    cpp11::integers const & i, 
    cpp11::integers const & p, 
    int const & nrow, int const & ncol,
    std::vector<double> & tx, 
    std::vector<int> & ti, 
    std::vector<int> & tp, 
    int const & threads);

// NOTe:  there is no formal definition of sparse matrix.
// input is column major, so i has the row ids, and p is per column.
// direct to stl vectors
template void _sp_transpose_par(
    cpp11::doubles const & x, 
    cpp11::integers const & i, 
    cpp11::doubles const & p, 
    int const & nrow, int const & ncol,
    std::vector<double> & tx, 
    std::vector<int> & ti, 
    std::vector<long> & tp, 
    int const & threads);

// NOTe:  there is no formal definition of sparse matrix.
// input is column major, so i has the row ids, and p is per column.
// direct to stl vectors
template void _sp_transpose_par(
    cpp11::doubles const & x, 
    cpp11::integers const & i, 
    cpp11::integers const & p, 
    int const & nrow, int const & ncol,
    double * tx, 
    int * ti, 
    int * tp, 
    int const & threads);

// NOTe:  there is no formal definition of sparse matrix.
// input is column major, so i has the row ids, and p is per column.
// direct to stl vectors
template void _sp_transpose_par(
    cpp11::doubles const & x, 
    cpp11::integers const & i, 
    cpp11::doubles const & p, 
    int const & nrow, int const & ncol,
    double * tx, 
    int * ti, 
    long * tp, 
    int const & threads);


// NOTe:  there is no formal definition of sparse matrix.
// input is column major, so i has the row ids, and p is per column.
// direct to stl vectors
template void _sp_transpose(
    cpp11::doubles const & x, 
    cpp11::integers const & i, 
    cpp11::integers const & p, 
    int const & nrow, int const & ncol, 
    std::vector<double> & tx, 
    std::vector<int> & ti, 
    std::vector<int> & tp, 
    int const & threads);

// NOTe:  there is no formal definition of sparse matrix.
// input is column major, so i has the row ids, and p is per column.
// direct to stl vectors
template void _sp_transpose(
    cpp11::doubles const & x, 
    cpp11::integers const & i, 
    cpp11::doubles const & p, 
    int const & nrow, int const & ncol, 
    std::vector<double> & tx, 
    std::vector<int> & ti, 
    std::vector<long> & tp, 
    int const & threads);

// NOTe:  there is no formal definition of sparse matrix.
// input is column major, so i has the row ids, and p is per column.
// direct to stl vectors
template void _sp_transpose(
    cpp11::doubles const & x, 
    cpp11::integers const & i, 
    cpp11::integers const & p, 
    int const & nrow, int const & ncol, 
    double * tx, 
    int * ti, 
    int * tp, 
    int const & threads);

// NOTe:  there is no formal definition of sparse matrix.
// input is column major, so i has the row ids, and p is per column.
// direct to stl vectors
template void _sp_transpose(
    cpp11::doubles const & x, 
    cpp11::integers const & i, 
    cpp11::doubles const & p, 
    int const & nrow, int const & ncol, 
    double * tx, 
    int * ti, 
    long * tp, 
    int const & threads);

// no names.
template cpp11::writable::doubles_matrix<cpp11::by_column> _sp_to_dense(cpp11::doubles const & x, 
    cpp11::integers const & i, cpp11::integers const & p, 
    int const & nrow, int const & ncol, int const & threads);

template cpp11::writable::doubles_matrix<cpp11::by_column> _sp_to_dense(cpp11::doubles const & x, 
    cpp11::integers const & i, cpp11::doubles const & p, 
    int const & nrow, int const & ncol, int const & threads);

template cpp11::writable::doubles_matrix<cpp11::by_column> _sp_to_dense_transposed(
    cpp11::doubles const & x, 
    cpp11::integers const & i, 
    cpp11::integers const & p, int const & nrow, int const & ncol, int const & threads);

template cpp11::writable::doubles_matrix<cpp11::by_column> _sp_to_dense_transposed(
    cpp11::doubles const & x, 
    cpp11::integers const & i, 
    cpp11::doubles const & p, int const & nrow, int const & ncol, int const & threads);

template cpp11::writable::list _sp_rbind(
    cpp11::list_of<cpp11::doubles> const & xvecs, 
    cpp11::list_of<cpp11::integers> const & ivecs, 
    cpp11::list_of<cpp11::integers> const & pvecs, 
    cpp11::integers const & nrows, 
    cpp11::integers const & ncols, int const & threads);

template cpp11::writable::list _sp_rbind(
    cpp11::list_of<cpp11::doubles> const & xvecs, 
    cpp11::list_of<cpp11::integers> const & ivecs, 
    cpp11::list_of<cpp11::doubles> const & pvecs, 
    cpp11::integers const & nrows, 
    cpp11::integers const & ncols, int const & threads);

template cpp11::writable::list _sp_cbind(
    cpp11::list_of<cpp11::doubles> const & xvecs, 
    cpp11::list_of<cpp11::integers> const & ivecs, 
    cpp11::list_of<cpp11::integers> const & pvecs, 
    cpp11::integers const & nrows, 
    cpp11::integers const & ncols, int const & threads);

template cpp11::writable::list _sp_cbind(
    cpp11::list_of<cpp11::doubles> const & xvecs, 
    cpp11::list_of<cpp11::integers> const & ivecs, 
    cpp11::list_of<cpp11::doubles> const & pvecs, 
    cpp11::integers const & nrows, 
    cpp11::integers const & ncols, int const & threads);


// csc
template cpp11::writable::doubles _sp_colsums(
    cpp11::doubles const & x, 
    cpp11::integers const & p, 
    int const & ncol, 
    int const & threads);

template cpp11::writable::doubles _sp_colsums(
    cpp11::doubles const & x, 
    cpp11::doubles const & p, 
    int const & ncol, 
    int const & threads);

// csc
template cpp11::writable::doubles _sp_rowsums(
    cpp11::doubles const & x, 
    cpp11::integers const & i, 
    int const & nrow, size_t const & nzcount, 
    int const & threads);



// // NOTE:  there is no formal definition of sparse matrix.
// // input is column major, so i has the row ids, and p is per column.
// // direct to stl vectors
// template void csc_transpose_csc(
//     cpp11::doubles::const_iterator x, 
//     cpp11::integers::const_iterator i, 
//     cpp11::integers::const_iterator p, 
//     int const & nrow, int const & ncol,
//     cpp11::writable::doubles::iterator tx, 
//     cpp11::writable::integers::iterator ti, 
//     cpp11::writable::integers::iterator tp, 
//     int const & threads);

// template void csc_transpose_csc(
//     cpp11::doubles::const_iterator x, 
//     cpp11::integers::const_iterator i, 
//     cpp11::doubles::const_iterator p, 
//     int const & nrow, int const & ncol,
//     cpp11::writable::doubles::iterator tx, 
//     cpp11::writable::integers::iterator ti, 
//     cpp11::writable::doubles::iterator tp, 
//     int const & threads);

// // no names.
// template void csc_to_dense_c(
//     cpp11::doubles::const_iterator x, 
//     cpp11::integers::const_iterator i, 
//     cpp11::integers::const_iterator p, 
//     int const & nrow, int const & ncol, 
//     std::vector<double>::iterator out,
//     int const & threads);

// // no names.
// template void csc_to_dense_c(
//     cpp11::doubles::const_iterator x, 
//     cpp11::integers::const_iterator i, 
//     cpp11::doubles::const_iterator p, 
//     int const & nrow, int const & ncol, 
//     std::vector<double>::iterator out,
//     int const & threads);

// template void csc_to_dense_transposed_c(
//     cpp11::doubles::const_iterator x, 
//     cpp11::integers::const_iterator i, 
//     cpp11::integers::const_iterator p, 
//     int const & nrow, int const & ncol, 
//     std::vector<double>::iterator out,
//     int const & threads);

// template void csc_to_dense_transposed_c(
//     cpp11::doubles::const_iterator x, 
//     cpp11::integers::const_iterator i, 
//     cpp11::doubles::const_iterator p, 
//     int const & nrow, int const & ncol, 
//     std::vector<double>::iterator out,
//     int const & threads);

// return total number of rows.
template int csc_rbind(
    std::vector<cpp11::doubles::const_iterator> const & xvecs, 
    std::vector<cpp11::integers::const_iterator> const & ivecs, 
    std::vector<cpp11::integers::const_iterator> const & pvecs, 
    std::vector<int> const & nrows, 
    std::vector<int> const & ncols, 
    cpp11::writable::doubles::iterator ox,
    cpp11::writable::integers::iterator oi,
    cpp11::writable::integers::iterator op,
    int const & threads);

// return total number of rows.
template int csc_rbind(
    std::vector<cpp11::doubles::const_iterator> const & xvecs, 
    std::vector<cpp11::integers::const_iterator> const & ivecs, 
    std::vector<cpp11::integers::const_iterator> const & pvecs, 
    std::vector<int> const & nrows, 
    std::vector<int> const & ncols, 
    cpp11::writable::doubles::iterator ox,
    cpp11::writable::integers::iterator oi,
    cpp11::writable::doubles::iterator op,
    int const & threads);

// return total number of rows.
template int csc_rbind(
    std::vector<cpp11::doubles::const_iterator> const & xvecs, 
    std::vector<cpp11::integers::const_iterator> const & ivecs, 
    std::vector<cpp11::doubles::const_iterator> const & pvecs, 
    std::vector<int> const & nrows, 
    std::vector<int> const & ncols, 
    cpp11::writable::doubles::iterator ox,
    cpp11::writable::integers::iterator oi,
    cpp11::writable::doubles::iterator op,
    int const & threads);

// return total number of rows.
template int csc_rbind(
    std::vector<std::vector<double>::const_iterator> const & xvecs, 
    std::vector<std::vector<int>::const_iterator> const & ivecs, 
    std::vector<std::vector<int>::const_iterator> const & pvecs, 
    std::vector<int> const & nrows, 
    std::vector<int> const & ncols, 
    cpp11::writable::doubles::iterator ox,
    cpp11::writable::integers::iterator oi,
    cpp11::writable::integers::iterator op,
    int const & threads);

// return total number of rows.
template int csc_rbind(
    std::vector<std::vector<double>::const_iterator> const & xvecs, 
    std::vector<std::vector<int>::const_iterator> const & ivecs, 
    std::vector<std::vector<int>::const_iterator> const & pvecs, 
    std::vector<int> const & nrows, 
    std::vector<int> const & ncols, 
    cpp11::writable::doubles::iterator ox,
    cpp11::writable::integers::iterator oi,
    cpp11::writable::doubles::iterator op,
    int const & threads);

// return total number of rows.
template int csc_rbind(
    std::vector<std::vector<double>::const_iterator> const & xvecs, 
    std::vector<std::vector<int>::const_iterator> const & ivecs, 
    std::vector<std::vector<double>::const_iterator> const & pvecs, 
    std::vector<int> const & nrows, 
    std::vector<int> const & ncols, 
    cpp11::writable::doubles::iterator ox,
    cpp11::writable::integers::iterator oi,
    cpp11::writable::doubles::iterator op,
    int const & threads);



template int csc_rbind_vec(
    std::vector<cpp11::doubles> const & xvecs, 
    std::vector<cpp11::integers> const & ivecs, 
    std::vector<cpp11::integers> const & pvecs, 
    cpp11::integers const & nrows, 
    cpp11::integers const & ncols, 
    std::vector<double> & ox,
    std::vector<int> & oi,
    std::vector<int> & op,
    int const & threads);

// return total number of rows.
template int csc_rbind_vec(
    std::vector<cpp11::doubles> const & xvecs, 
    std::vector<cpp11::integers> const & ivecs, 
    std::vector<cpp11::integers> const & pvecs, 
    cpp11::integers const & nrows, 
    cpp11::integers const & ncols, 
    std::vector<double> & ox,
    std::vector<int> & oi,
    std::vector<double> & op,
    int const & threads);

// return total number of rows.
template int csc_rbind_vec(
    std::vector<cpp11::doubles> const & xvecs, 
    std::vector<cpp11::integers> const & ivecs, 
    std::vector<cpp11::doubles> const & pvecs, 
    cpp11::integers const & nrows, 
    cpp11::integers const & ncols, 
    std::vector<double> & ox,
    std::vector<int> & oi,
    std::vector<double> & op,
    int const & threads);

template int csc_rbind_vec(
    std::vector<cpp11::doubles> const & xvecs, 
    std::vector<cpp11::integers> const & ivecs, 
    std::vector<cpp11::integers> const & pvecs, 
    cpp11::integers const & nrows, 
    cpp11::integers const & ncols, 
    cpp11::writable::doubles & ox,
    cpp11::writable::integers & oi,
    cpp11::writable::integers &  op,
    int const & threads);

// return total number of rows.
template int csc_rbind_vec(
    std::vector<cpp11::doubles> const & xvecs, 
    std::vector<cpp11::integers> const & ivecs, 
    std::vector<cpp11::integers> const & pvecs, 
    cpp11::integers const & nrows, 
    cpp11::integers const & ncols, 
    cpp11::writable::doubles & ox,
    cpp11::writable::integers & oi,
    cpp11::writable::doubles & op,
    int const & threads);

// return total number of rows.
template int csc_rbind_vec(
    std::vector<cpp11::doubles> const & xvecs, 
    std::vector<cpp11::integers> const & ivecs, 
    std::vector<cpp11::doubles> const & pvecs, 
    cpp11::integers const & nrows, 
    cpp11::integers const & ncols, 
    cpp11::writable::doubles & ox,
    cpp11::writable::integers & oi,
    cpp11::writable::doubles & op,
    int const & threads);


// template int csc_rbind_vec(
//     cpp11::list_of<cpp11::doubles> const & xvecs, 
//     cpp11::list_of<cpp11::integers> const & ivecs, 
//     cpp11::list_of<cpp11::integers> const & pvecs, 
//     cpp11::integers const & nrows, 
//     cpp11::integers const & ncols, 
//     std::vector<double> & ox,
//     std::vector<int> & oi,
//     std::vector<int> & op,
//     int const & threads);

// // return total number of rows.
// template int csc_rbind_vec(
//     cpp11::list_of<cpp11::doubles> const & xvecs, 
//     cpp11::list_of<cpp11::integers> const & ivecs, 
//     cpp11::list_of<cpp11::integers> const & pvecs, 
//     cpp11::integers const & nrows, 
//     cpp11::integers const & ncols, 
//     std::vector<double> & ox,
//     std::vector<int> & oi,
//     std::vector<double> & op,
//     int const & threads);

// // return total number of rows.
// template int csc_rbind_vec(
//     cpp11::list_of<cpp11::doubles> const & xvecs, 
//     cpp11::list_of<cpp11::integers> const & ivecs, 
//     cpp11::list_of<cpp11::doubles> const & pvecs, 
//     cpp11::integers const & nrows, 
//     cpp11::integers const & ncols, 
//     std::vector<double> & ox,
//     std::vector<int> & oi,
//     std::vector<double> & op,
//     int const & threads);


// template int csc_rbind_vec(
//     cpp11::list_of<cpp11::doubles> const & xvecs, 
//     cpp11::list_of<cpp11::integers> const & ivecs, 
//     cpp11::list_of<cpp11::integers> const & pvecs, 
//     cpp11::integers const & nrows, 
//     cpp11::integers const & ncols, 
//     cpp11::writable::doubles & ox,
//     cpp11::writable::integers & oi,
//     cpp11::writable::integers &  op,
//     int const & threads);

// // return total number of rows.
// template int csc_rbind_vec(
//     cpp11::list_of<cpp11::doubles> const & xvecs, 
//     cpp11::list_of<cpp11::integers> const & ivecs, 
//     cpp11::list_of<cpp11::integers> const & pvecs, 
//     cpp11::integers const & nrows, 
//     cpp11::integers const & ncols, 
//     cpp11::writable::doubles & ox,
//     cpp11::writable::integers & oi,
//     cpp11::writable::doubles & op,
//     int const & threads);

// // return total number of rows.
// template int csc_rbind_vec(
//     cpp11::list_of<cpp11::doubles> const & xvecs, 
//     cpp11::list_of<cpp11::integers> const & ivecs, 
//     cpp11::list_of<cpp11::doubles> const & pvecs, 
//     cpp11::integers const & nrows, 
//     cpp11::integers const & ncols, 
//     cpp11::writable::doubles & ox,
//     cpp11::writable::integers & oi,
//     cpp11::writable::doubles & op,
//     int const & threads);

// return total number of columns.
template int csc_cbind(
    std::vector<cpp11::doubles::const_iterator> const & xvecs, 
    std::vector<cpp11::integers::const_iterator> const & ivecs, 
    std::vector<cpp11::integers::const_iterator> const & pvecs, 
    std::vector<int> const & nrows, 
    std::vector<int> const & ncols, 
    cpp11::writable::doubles::iterator ox,
    cpp11::writable::integers::iterator oi,
    cpp11::writable::integers::iterator op,
    int const & threads);

// return total number of columns.
template int csc_cbind(
    std::vector<cpp11::doubles::const_iterator> const & xvecs, 
    std::vector<cpp11::integers::const_iterator> const & ivecs, 
    std::vector<cpp11::integers::const_iterator> const & pvecs, 
    std::vector<int> const & nrows, 
    std::vector<int> const & ncols, 
    cpp11::writable::doubles::iterator ox,
    cpp11::writable::integers::iterator oi,
    cpp11::writable::doubles::iterator op,
    int const & threads);

// return total number of columns.
template int csc_cbind(
    std::vector<cpp11::doubles::const_iterator> const & xvecs, 
    std::vector<cpp11::integers::const_iterator> const & ivecs, 
    std::vector<cpp11::doubles::const_iterator> const & pvecs, 
    std::vector<int> const & nrows, 
    std::vector<int> const & ncols, 
    cpp11::writable::doubles::iterator ox,
    cpp11::writable::integers::iterator oi,
    cpp11::writable::doubles::iterator op,
    int const & threads);

// return total number of columns.
template int csc_cbind(
    std::vector<std::vector<double>::const_iterator> const & xvecs, 
    std::vector<std::vector<int>::const_iterator> const & ivecs, 
    std::vector<std::vector<int>::const_iterator> const & pvecs, 
    std::vector<int> const & nrows, 
    std::vector<int> const & ncols, 
    cpp11::writable::doubles::iterator ox,
    cpp11::writable::integers::iterator oi,
    cpp11::writable::integers::iterator op,
    int const & threads);

// return total number of columns.
template int csc_cbind(
    std::vector<std::vector<double>::const_iterator> const & xvecs, 
    std::vector<std::vector<int>::const_iterator> const & ivecs, 
    std::vector<std::vector<int>::const_iterator> const & pvecs, 
    std::vector<int> const & nrows, 
    std::vector<int> const & ncols, 
    cpp11::writable::doubles::iterator ox,
    cpp11::writable::integers::iterator oi,
    cpp11::writable::doubles::iterator op,
    int const & threads);

// return total number of columns.
template int csc_cbind(
    std::vector<std::vector<double>::const_iterator> const & xvecs, 
    std::vector<std::vector<int>::const_iterator> const & ivecs, 
    std::vector<std::vector<double>::const_iterator> const & pvecs, 
    std::vector<int> const & nrows, 
    std::vector<int> const & ncols, 
    cpp11::writable::doubles::iterator ox,
    cpp11::writable::integers::iterator oi,
    cpp11::writable::doubles::iterator op,
    int const & threads);


// return total number of columns.
template int csc_cbind_vec(
    std::vector<cpp11::doubles> const & xvecs, 
    std::vector<cpp11::integers> const & ivecs, 
    std::vector<cpp11::integers> const & pvecs, 
    cpp11::integers const & nrows, 
    cpp11::integers const & ncols, 
    cpp11::writable::doubles & ox,
    cpp11::writable::integers & oi,
    cpp11::writable::integers & op,
    int const & threads);

// return total number of columns.
template int csc_cbind_vec(
    std::vector<cpp11::doubles> const & xvecs, 
    std::vector<cpp11::integers> const & ivecs, 
    std::vector<cpp11::integers> const & pvecs, 
    cpp11::integers const & nrows, 
    cpp11::integers const & ncols, 
    cpp11::writable::doubles & ox,
    cpp11::writable::integers & oi,
    cpp11::writable::doubles & op,
    int const & threads);

// return total number of columns.
template int csc_cbind_vec(
    std::vector<cpp11::doubles> const & xvecs, 
    std::vector<cpp11::integers> const & ivecs, 
    std::vector<cpp11::doubles> const & pvecs, 
    cpp11::integers const & nrows, 
    cpp11::integers const & ncols, 
    cpp11::writable::doubles & ox,
    cpp11::writable::integers & oi,
    cpp11::writable::doubles & op,
    int const & threads);

// // return total number of columns.
// template int csc_cbind_vec(
//     cpp11::list_of<cpp11::doubles> const & xvecs, 
//     cpp11::list_of<cpp11::integers> const & ivecs, 
//     cpp11::list_of<cpp11::integers> const & pvecs, 
//     cpp11::integers const & nrows, 
//     cpp11::integers const & ncols, 
//     cpp11::writable::doubles & ox,
//     cpp11::writable::integers & oi,
//     cpp11::writable::integers & op,
//     int const & threads);

// // return total number of columns.
// template int csc_cbind_vec(
//     cpp11::list_of<cpp11::doubles> const & xvecs, 
//     cpp11::list_of<cpp11::integers> const & ivecs, 
//     cpp11::list_of<cpp11::integers> const & pvecs, 
//     cpp11::integers const & nrows, 
//     cpp11::integers const & ncols, 
//     cpp11::writable::doubles & ox,
//     cpp11::writable::integers & oi,
//     cpp11::writable::doubles & op,
//     int const & threads);

// // return total number of columns.
// template int csc_cbind_vec(
//     cpp11::list_of<cpp11::doubles> const & xvecs, 
//     cpp11::list_of<cpp11::integers> const & ivecs, 
//     cpp11::list_of<cpp11::doubles> const & pvecs, 
//     cpp11::integers const & nrows, 
//     cpp11::integers const & ncols, 
//     cpp11::writable::doubles & ox,
//     cpp11::writable::integers & oi,
//     cpp11::writable::doubles & op,
//     int const & threads);

// csc
template void csc_colsums_iter(
    cpp11::doubles::const_iterator x, 
    cpp11::integers::const_iterator p, 
    int const & ncol, 
    cpp11::writable::doubles::iterator out,
    int const & threads);

template void csc_colsums_iter(
    cpp11::doubles::const_iterator x, 
    cpp11::doubles::const_iterator p, 
    int const & ncol, 
    cpp11::writable::doubles::iterator out,
    int const & threads);


template void csc_colsums_iter(
    cpp11::doubles::const_iterator x, 
    cpp11::integers::const_iterator p, 
    int const & ncol, 
    std::vector<double>::iterator out,
    int const & threads);

template void csc_colsums_iter(
    cpp11::doubles::const_iterator x, 
    cpp11::doubles::const_iterator p, 
    int const & ncol, 
    std::vector<double>::iterator out,
    int const & threads);


template void csc_colsums_vec(
    cpp11::doubles const & x, 
    cpp11::integers const & p, 
    int const & ncol, 
    cpp11::writable::doubles & out,
    int const & threads);


template void csc_colsums_vec(
    cpp11::doubles const & x, 
    cpp11::doubles const & p, 
    int const & ncol, 
    cpp11::writable::doubles & out,
    int const & threads);



// NOTe:  there is no formal definition of sparse matrix.
// input is column major, so i has the row ids, and p is per column.
template void csc_rowsums_iter(
    cpp11::doubles::const_iterator x, 
    cpp11::integers::const_iterator i, 
    int const & nrow, size_t const & nzcount, 
    cpp11::writable::doubles::iterator out,
    int const & threads);



// NOTe:  there is no formal definition of sparse matrix.
// input is column major, so i has the row ids, and p is per column.
template void csc_rowsums_vec(
    cpp11::doubles const & x, 
    cpp11::integers const & i, 
    int const & nrow,
    cpp11::writable::doubles & out,
    int const & threads);

