#include <iostream>
#include <stdlib.h>  // for rand()
#include <string>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <cmath>        //for std::abs
#include "Distributions.h"
#include "DoubleMatrix.h"
#include "DoubleRow.h"
#include "DoubleColumn.h"

using namespace std;

/**************************************************************************
 *
 * Constructors and factory functions
 *
 **************************************************************************/

// default constructor. Must be defined because we have other constructors.
DoubleMatrix::DoubleMatrix(){};  

/** Creates an n times m matrix filled with zeros.  */
DoubleMatrix::DoubleMatrix(int rows, int columns) {
    arma_mat = arma::mat(rows, columns, arma::fill::zeros);
}

/**
 * Creates an p_rows times p_columns matrix filled from the given source pointer. Data
 * is read sequentially from the source. If the boolean parameter is true the obtained
 * data are stored in the matrix row by row (used to read from C++ arrays). Otherwise,
 * the obtained data are stored in the matrix filling column by column (used to read
 * matrices stored in column-major format, e.g. Armadillo matrices, and also to transpose
 * another matrix).
 *
 * In the caller procedure use &a[0] if a is one-dimensional array,
 * &a[0][0] if a is two-dimensional array.
 */
DoubleMatrix::DoubleMatrix(int p_rows, int p_columns, const double* p_source, bool p_fill_row_by_row) {
    // arma_mat = arma::mat(p_rows, p_columns);
    if (p_fill_row_by_row) {
        // Copy row by row (row-major). Aramdillo fills column by columns
        // Therefore we fill the transposed matrix, and then transpose 
        // we use const_cast<double*>(p_source) because the constructor
        // expect double* when copy_aux_mem=false, instead of const double*.
        arma::mat tmp(const_cast<double*>(p_source), p_columns, p_rows, false,true);
        arma_mat = tmp.t();
    } else {
        // Create col-major view
        // Armadillo fills column-major, we don't need to transpose
        // we allow copy_aux_mem=true (default) to allow copying data to the new aram_mat
        // strict in this case (copy_aux_mem=true) is ignored. 
        arma_mat = arma::mat(p_source, p_rows, p_columns); 
    }
}

/** All the following constructors may be safely removed because they just copy arma_mat, 
 * which is handled automatically by the implicit copy/move constructors created by the 
 * compiler. They just exist because the older version had pointers, and required 
 * allocating memory for the internal data.   
*/

/** copy constructor: invoked by the compiler when 'other' is an lvalue
 * In this case the constructor should copy the internal components of (other) to (this) and keep the original (other) as it is*/
DoubleMatrix::DoubleMatrix(const DoubleMatrix& other)
: arma_mat(other.arma_mat)  {
	// std::cout << "matrix copy constructor\n";
}

/** move constructor: invoked by the compiler when 'other' is an rvalue.
 * In this case the constructor should move the internal components of (other) to (this)
 * DoubleMatrix&& is an rvalue reference to a DoubleMatrix object  */
DoubleMatrix::DoubleMatrix(DoubleMatrix&& other)
: arma_mat(other.arma_mat) {
    // std::cout << "matrix move constructor\n";
}


/** The assignment operator, used for assigning objects. 
 *  Works for both copy and move assignments.
 *  It is implemented to avoid assigning matrices/vectors of different sizes accedentally. 
 * Important:
 * 1. This assignment may involve implicit type conversion (row/column <-> matrix)
 * 2. Unlike other member functions, in this function we pass other by value:
 * The true genius of this method is what happens when the passed-by-value object to be assigned from
 * (the source object) is created as a temporary object in the assignment operator function. There are
 * three possibilities:
 *  i.  The source object was an lvalue: the source object’s copy constructor will be called to make a copy
 *      and the overload (=) will behave as a copy assignment operator.
 *  ii. The source object was an rvalue: nothing is called; the compiler elides the temporary copy and passes
 *      the object by value; the overload (=) will behave as a move assignment operator.
 *  iii.The source object was an rvalue: The source object’s move constructor will be called and the overload (=)
 *      will behave as a move assignment operator. */
DoubleMatrix& DoubleMatrix::operator=(DoubleMatrix other) {
	// std::cout << "matrix assignment operation\n";
	// supressed checking dimensions to allow assigning to matrices 
    // constructed by default constructors. Note that all matrix libraries don't restrict.
    // DoubleRow and DoubleColumns are handled differently using their conversion constructors
    // if (!(getRows() == other.getRows() && getColumns() == other.getColumns())) {
    //     std::cerr<< "matrix error: Unable to assign source matrix with "
    //     << other.getRows() <<" rows, " << other.getColumns()<< ", to matrix with " 
    //     << getRows() <<" rows, " << getColumns() << " columns\n";
    //     exit(1);
    // }
    arma_mat = other.arma_mat;
	return *this;
}


/** fill the matrix with random values uniformly in 0..1. */
DoubleMatrix& DoubleMatrix::fillRand() {
    arma_mat.randu();
    return *this;
}

/** fill the matrix with the given value */
DoubleMatrix& DoubleMatrix::fillWith(double v) {
    arma_mat.fill(v);
    return *this;
}

/** make this an identity (eye) matrix. It must be square matrix*/
DoubleMatrix& DoubleMatrix::makeItEye() {
    this->assertSquare();
    arma_mat.eye();
    return *this;
}

/**
 * Make it a diagonal matrix of arbitrary shape and set the diagonal according
 * to a passed vector d. length of d must be smaller than both rows and columns.
 */
DoubleMatrix& DoubleMatrix::makeItDiag(const DoubleMatrix&  d) {
    if (!(d.getRows()==1 || d.getColumns()==1)) {
        throw("matrix error: Unable to create diagonal matrix. The given vector is invalid");
    }
    if ( d.getLength()>getRows() || d.getLength()>getColumns()) {
        throw("matrix error: Unable to create diagonal matrix. The given vector is too large");
    }
    arma_mat.zeros();
    for (int i = 0; i < d.getLength(); i++) {
        arma_mat(i, i) = d.get(i);  
    }
    return *this;
}

/**
 * Concatenates this matrix to the given matrix b horizontally. 
 * Matrices must have identical numbers of rows.
 */
DoubleMatrix DoubleMatrix::concatHorizontally(const DoubleMatrix& b) const {
    if (getRows() != b.getRows()) {
        ostringstream oss;
        oss << "matrix error: The given matrix has a different number of rows (" << b.getRows() <<" != " << getRows()  <<".";
        string msg = oss.str();
        throw(msg);
    }
    DoubleMatrix result;
    result.arma_mat = arma::join_rows(arma_mat, b.arma_mat);
    return result;
}

/**
 * Concatenates this matrix to the given matrix b vertically. Matrices must have identical
 * numbers of columns.
 */
DoubleMatrix DoubleMatrix::concatVertically(const DoubleMatrix& b) const {
    if (getColumns() != b.getColumns()) {
        ostringstream oss;
        oss << "matrix error: The given matrix has a different number of columns (" << getColumns() <<" != " << b.getColumns()  <<".";
        string msg = oss.str();
        throw(msg);
    }
    DoubleMatrix result;
    result.arma_mat = arma::join_cols(arma_mat, b.arma_mat);
    return result;
}

/**************************************************************************
 * restricting the matrix to specified rows	and columns
 */

/**
 * Restrict the matrix in-place to the rows passed in the vector rowsIndexes
 * we assume that rowsIndexes elements are in ascending order
 * */
DoubleMatrix& DoubleMatrix::selectiRows(const std::vector<int>& rowsIndexes) {
    arma::uvec indices(rowsIndexes.size());
    for (size_t i = 0; i < rowsIndexes.size(); i++) {
        if (rowsIndexes[i] < 0 || rowsIndexes[i] >= getRows()) {
            std::cerr << "unable to select row: index is out of range";
            exit(0);
        }
        indices[i] = rowsIndexes[i];
    }
    arma_mat = arma_mat.rows(indices);
    return *this;
}

/** Restrict the matrix in-place to the columns passed in the vector columnIndexes */
DoubleMatrix& DoubleMatrix::selectiColumns(const std::vector<int>& columnsIndexes) {
    arma::uvec indices(columnsIndexes.size());
    for (size_t i = 0; i < columnsIndexes.size(); i++) {
        if (columnsIndexes[i] < 0 || columnsIndexes[i] >= getColumns()) {
            std::cerr << "unable to select column: index is out of range";
            exit(0);
        }
        indices[i] = columnsIndexes[i];
    }
    arma_mat = arma_mat.cols(indices);
    return *this;
}

/** Restrict the matrix to the rows passed in the vector rowsIndexes */
DoubleMatrix DoubleMatrix::selectRows(const std::vector<int>& rowsIndexes) const {
    DoubleMatrix result = *this;
    return result.selectiRows(rowsIndexes);
}

/** returns row r of this matrix. */
DoubleRow DoubleMatrix::selectRow(int r) const {
    if (r<0 | r >= getRows()) {
        std::cerr << "unable to select row: index is out of range";
        exit(0);
    }
    DoubleRow result;
    result.arma_mat = arma_mat.row(r);
    return result;
}

/** Restrict the matrix to the columns passed in the vector columnsIndexes  */
DoubleMatrix DoubleMatrix::selectColumns(const std::vector<int>& columnsIndexes) const {
    DoubleMatrix result = *this;
    return result.selectiColumns(columnsIndexes);
}

/** return column c of this matrix. */
DoubleColumn DoubleMatrix::selectColumn(int c) const {
    if (c<0 | c >= getColumns()) {
        std::cerr << "unable to select column: index is out of range";
        exit(0);
    }
    DoubleColumn result;
    result.arma_mat = arma_mat.col(c);
    return result;
}

/**************************************************************************
 * Basic operations (copying, transposing, resizing, element access)
 */
/** Return the transpose of this matrix. */
DoubleMatrix DoubleMatrix::transpose() const {
    DoubleMatrix result(getColumns(), getRows());
    result.arma_mat = arma_mat.t();
    return result;
}

/** Resize the matrix. All elements will be set to zero. */
DoubleMatrix& DoubleMatrix::resize(int newRows, int newColumns) {
    arma_mat.set_size(newRows, newColumns);
    arma_mat.zeros();
    return *this;
}

/** Reshape the matrix. Number of elements must not change. */
DoubleMatrix& DoubleMatrix::reshape(int newRows, int newColumns) {
    if (getLength() != newRows * newColumns) {
        throw ("matrix error: Cannot reshape matrix. Number of elements must not change.");
    }
    arma_mat.reshape(newRows, newColumns);
    return *this;
}

/** Put (copy) matrix a into this. The two matrices must have the same dimensions.
 *  Note that it is in-place operation
 */
DoubleMatrix& DoubleMatrix::put(const DoubleMatrix& a) {
    assertSameDimensions(a);
    arma_mat = a.arma_mat;
    return *this;
}

/** Swap two columns of a matrix. */
DoubleMatrix& DoubleMatrix::swapColumns(int i, int j) {
    arma_mat.swap_cols(i, j);
    return *this;
}

/** Swap two rows of a matrix. */
DoubleMatrix& DoubleMatrix::swapRows(int i, int j) {
    arma_mat.swap_rows(i, j);
    return *this;
}

/** Set matrix element */
DoubleMatrix& DoubleMatrix::put(int rowIndex, int columnIndex, double value) {
    arma_mat(rowIndex, columnIndex) = value;
    return *this;
}

/** Retrieve matrix element */
double DoubleMatrix::get(int rowIndex, int columnIndex) const {
    return arma_mat(rowIndex, columnIndex);
}

/** Get the linear index of an element */
int DoubleMatrix::index(int rowIndex, int columnIndex) const {
    return rowIndex*getColumns() + columnIndex;
}

/** Compute the row index of a linear index. */
int DoubleMatrix::rowIndex(int i) const {
    return i / getColumns();
}

/** Compute the column index of a linear index. */
int DoubleMatrix::columnIndex(int i) const {
    return i % getColumns();
}

/** Get a matrix element (using linear indexing). */
double DoubleMatrix::get(int i) const {
    return arma_mat(rowIndex(i), columnIndex(i));
}

/** Set a matrix element (using linear indexing). */
DoubleMatrix& DoubleMatrix::put(int i, double v) {
    arma_mat(rowIndex(i), columnIndex(i)) = v;
    return *this;
}

/** Get number of rows. */
int DoubleMatrix::getRows() const {
    return arma_mat.n_rows;
}

/** Get number of columns. */
int DoubleMatrix::getColumns() const {
    return arma_mat.n_cols;
}

/** Get total number of elements. */
int DoubleMatrix::getLength() const {
    return arma_mat.n_elem;
}

/** Get the diagonal of this matrix as a column vector */
DoubleColumn DoubleMatrix::getDiag() const {
    DoubleColumn result(getRows());
    result.arma_mat = arma_mat.diag();
    return result;
}

/**************************************************************************
 * Arithmetic Operations
 */

/** Add matrix a to this (in-place) and return this. */
DoubleMatrix& DoubleMatrix::addi(const DoubleMatrix& a) {
    assertSameLength(a);
    arma_mat += a.arma_mat;
    return *this;
}

/** Add scalar v (element-wise) to this (in-place) and return this. */
DoubleMatrix& DoubleMatrix::addi(double v) {
    arma_mat += v;
    return *this;
}


/** Subtract matrix a from this (in-place). */
DoubleMatrix& DoubleMatrix::subi(const DoubleMatrix& a) {
    assertSameLength(a);
    arma_mat -= a.arma_mat;
    return *this;
}

/** Subtract scalar v (element-wise) from this (in-place) and return this. */
DoubleMatrix& DoubleMatrix::subi(double v) {
    arma_mat -= v;
    return *this;
}

/** Elementwise multiplication (in-place): multiply this by a and return this*/
DoubleMatrix& DoubleMatrix::muli(const DoubleMatrix& a) {
    assertSameLength(a);
    arma_mat %= a.arma_mat;
    return *this;
}

/** Elementwise multiplication with a scalar (in-place): multiply this by a and return this. */
DoubleMatrix& DoubleMatrix::muli(double v) {
    arma_mat *= v;
    return *this;
}

/** Matrix-matrix multiplication (in-place): this = this.a 
 * Note that the multiplier matrix a must be square to allow in-place operation for this.
*/
DoubleMatrix& DoubleMatrix::mmuli(const DoubleMatrix& a) {
    assertMultipliesWith(a);
    if (!a.isSquare()) {
        throw "matrix error: Unable to multiply in place: the right-side matrix must be square";
    }
    arma_mat = arma_mat * a.arma_mat;
    return *this;
}

/** Matrix-multiply by a matrix: result = this * a.   */
DoubleMatrix DoubleMatrix::mmul(const DoubleMatrix& a) const {
    assertMultipliesWith(a);
    DoubleMatrix result;
    result.arma_mat = arma_mat * a.arma_mat;
    return result;
}

/** Elementwise division by a (in-place). */
DoubleMatrix& DoubleMatrix::divi(const DoubleMatrix& a) {
    assertSameLength(a);
    arma_mat /= a.arma_mat;
    return *this;
}

/** Elementwise division with a scalar (in-place). */
DoubleMatrix& DoubleMatrix::divi(double a) {
    arma_mat /= a;
    return *this;
}


/** Negate each element (in-place) */
DoubleMatrix& DoubleMatrix::negi() {
    arma_mat = -arma_mat;
    return *this;
}

/** Negate each element, and set the result to matrix result */
DoubleMatrix DoubleMatrix::neg() const {
    DoubleMatrix result(getRows(), getColumns());
    result.arma_mat = -arma_mat;
    return result;
}

/** Maps zero to 1.0 and all non-zero values to 0.0 (in-place). */
DoubleMatrix& DoubleMatrix::noti() {
	for (int i = 0; i < getLength(); i++) {
		put(i, get(i) == 0.0 ? 1.0 : 0.0);
	}
    return *this;
}

/** Maps zero to 1.0 and all non-zero values to 0.0. Set the result matrix*/
DoubleMatrix DoubleMatrix::nott() const {
	DoubleMatrix result(getRows(), getColumns());
	return result.put(*this).noti();
}

/** Maps zero to 0.0 and all non-zero values to 1.0 (in-place). */
DoubleMatrix& DoubleMatrix::truthi() {
    arma::umat cond = arma::conv_to<arma::umat>::from(arma_mat != 0);
    arma_mat = arma::conv_to<arma::mat>::from(cond);
    return *this;
}

/** Maps zero to 0.0 and all non-zero values to 1.0. Set the result matrix*/
DoubleMatrix DoubleMatrix::truth() const {
	DoubleMatrix result(getRows(),getColumns());
	return result.put(*this).truthi();
}

/****************************************************************
 * Logical operations
 */
/** Returns the minimal element of the matrix. */
double DoubleMatrix::min() const {
    if (isEmpty()) {
        return std::numeric_limits<double>::infinity();
    }
    return arma_mat.min();
}

/**
 * Returns the linear index of the minimal element. If there are
 * more than one elements with this value, the first one is returned.
 * if the matrix is empty or have all NaN numbers returns -1
 */
int DoubleMatrix::argmin() const {
	if (isEmpty()) {
		return -1;
	}
	double v = std::numeric_limits<double>::infinity();
	int a = -1;
	for (int i = 0; i < getLength(); i++) {
		if (!isnan(get(i)) && get(i) < v) {
			v = get(i);
			a = i;
		}
	}
	return a;
}

/**
 * Computes the minimum between this matrix and the matrix a. Returns the smaller of the
 * corresponding elements (in-place).
 */
DoubleMatrix& DoubleMatrix::mini(const DoubleMatrix& a) {
    assertSameDimensions(a);
    arma_mat = arma::min(arma_mat, a.arma_mat);
    return *this;
}

/**
 * Computes the minimum between 'this' matrix and matrix a. Set result to
 * the smaller of the corresponding elements in a new matrix.
 */
DoubleMatrix DoubleMatrix::min(const DoubleMatrix& a) const {
    DoubleMatrix result(getRows(), getColumns());
    result.arma_mat = arma::min(arma_mat, a.arma_mat);
    return result;
}

/**
 * Computes the minimum between every element of this matrix and the value v.
 * Returns the smaller of the corresponding elements (in-place).
 */
DoubleMatrix& DoubleMatrix::mini(double v) {
    arma_mat = arma::min(arma_mat, v);
    return *this;
}

/**
 * Computes the minimum between every element of this matrix and the value v.
 * Returns the smaller of the corresponding elements (in-place).
 */
DoubleMatrix DoubleMatrix::min(double v) const {
    DoubleMatrix result(getRows(), getColumns());
    result.arma_mat = arma::min(arma_mat, v);
    return result;
}

/** Returns the maximal element of the matrix. */
double DoubleMatrix::max() const {
    if (isEmpty()) {
        return -std::numeric_limits<double>::infinity();
    }
    return arma_mat.max();
}

/**
 * Returns the linear index of the maximal element. If there are
 * more than one elements with this value, the first one is returned.
 * if the matrix is empty or have all NaN numbers returns -1
 */
int DoubleMatrix::argmax() const {
	if (isEmpty()) {
		return -1;
	}
	double v = -std::numeric_limits<double>::infinity();
	int a = -1;
	for (int i = 0; i < getLength(); i++) {
		if (!isnan(get(i)) && get(i) > v) {
			v = get(i);
			a = i;
		}
	}
	return a;
}

/**
 * Computes the maximum between this matrix and the matrix a. Returns the smaller of the
 * corresponding elements (in-place).
 */
DoubleMatrix& DoubleMatrix::maxi(const DoubleMatrix& a) {
    arma_mat = arma::max(arma_mat, a.arma_mat);
    return *this;
}

/**
 * Computes the maximum between 'this' matrix and matrix a. Set the matrix result
 * containing the maximum of the corresponding elements in a new matrix.
 */
DoubleMatrix DoubleMatrix::max(const DoubleMatrix& a) const {
    DoubleMatrix result(getRows(), getColumns());
    result.arma_mat = arma::max(arma_mat, a.arma_mat);
    return result;
}

/**
 * Computes the maximum between every element of this matrix and the value v.
 * Returns the largest of the corresponding elements (in-place).
 */
DoubleMatrix& DoubleMatrix::maxi(double v) {
    arma_mat = arma::max(arma_mat, v);
    return *this;
}

/**
 * Computes the maximum between every element of this matrix and the value v.
 * Set the matrix results to the bigger of the corresponding elements.
 */
DoubleMatrix DoubleMatrix::max(double v) const {
    DoubleMatrix result(getRows(), getColumns());
    result.arma_mat = arma::max(arma_mat, v);
    return result;
}

/** Computes the sum of all elements of the matrix. */
double DoubleMatrix::sum() const {
    return arma::accu(arma_mat);
}

/** Computes the product of all elements of the matrix */
double DoubleMatrix::prod() const {
    return arma::prod(arma::vectorise(arma_mat));
}

/**
 * Computes the mean value of all elements in the matrix,
 * that is, x.sum() / x.length.
 */
double DoubleMatrix::mean() const {
    return arma::mean(arma::vectorise(arma_mat));
}

/**
 * Computes the cumulative sum, that is, the element at index i contains the sum of all elements
 * up to i in linear addressing (in-place).
 */
DoubleMatrix& DoubleMatrix::cumulativeSumi() {
	double s = 0.0;
	for (int i = 0; i < getLength(); i++) {
		s += get(i);
		put(i, s);
	}
	return *this;
}

/**
 * Computes the cumulative sum, that is, the element at index i contains the sum of all elements
 * up to i in linear addressing. Set the result matrix.
 */
DoubleMatrix DoubleMatrix::cumulativeSum() const {
	DoubleMatrix result(getRows(), getColumns());
	return result.put(*this).cumulativeSumi();
}

/** Returns scalar product of this with matrix a: i.e. the sum of element-wise products */
double DoubleMatrix::dot(const DoubleMatrix& a) const {
    assertSameDimensions(a);
    return arma::dot(arma::vectorise(arma_mat), arma::vectorise(a.arma_mat));
}

/**
 * Computes the projection coefficient of other on this.
 *
 * The returned scalar times this is the orthogonal projection
 * of other on this.
 */
double DoubleMatrix::project(const DoubleMatrix& a) const {
    a.checkLength(getLength());
    double norm = 0, dot = 0;
    for (int i = 0; i < getLength(); i++) {
        double x = get(i);
        norm += x * x;
        dot += x * a.get(i);
    }
    return dot / norm;
}

/**
 * The Euclidean norm of the matrix as vector, also the Frobenius
 * norm of the matrix.
 */
double DoubleMatrix::norm2() const {
    return arma::norm(arma_mat, "fro");
}

/**
 * The maximum norm of the matrix (maximal absolute value of the elements).
 */
double DoubleMatrix::normmax() const {
    return arma::norm(arma_mat, "inf");
}

/**
 * The 1-norm of the matrix as vector (sum of absolute values of elements).
 */
double DoubleMatrix::norm1() const {
    return arma::norm(arma::vectorise(arma_mat), 1);
}

/**
 * Returns the squared (Euclidean) distance.
 */
double DoubleMatrix::squaredDistance(const DoubleMatrix& a) const {
    a.checkLength(getLength());
    arma::mat diff = arma_mat - a.arma_mat;
    return arma::accu(diff % diff);
}

/**
 * Returns the Euclidean distance.
 */
double DoubleMatrix::distance2(const DoubleMatrix& a) const {
    return sqrt(squaredDistance(a));
}

/**
 * Returns the (1-norm) distance.
 */
double DoubleMatrix::distance1(const DoubleMatrix& a) const {
    a.checkLength(getLength());
    arma::mat diff = arma_mat - a.arma_mat;
    return arma::norm(arma::vectorise(diff), 1);
}

/** Return a row vector containing the sum of every column */
DoubleRow DoubleMatrix::columnSums() const {
    DoubleRow result(getColumns());
    result.arma_mat = arma::sum(arma_mat, 0); // sum along dimension 0 (columns)
    return result;
}

/** returns a row vector that contains the mean of every column. */
DoubleRow DoubleMatrix::columnMeans() const {
    DoubleRow result(getColumns());
    result.arma_mat = arma::mean(arma_mat, 0); // mean along dimension 0 (columns)
    return result;
}

/** returns a column vector that contains the sum of every row. */
DoubleColumn DoubleMatrix::rowSums() const {
    DoubleColumn result(getRows());
    result.arma_mat = arma::sum(arma_mat, 1);
    return result;
}

/** returns a column vector that contains the mean of every row. */
DoubleColumn DoubleMatrix::rowMeans() const {
    DoubleColumn result(getRows());
    result.arma_mat = arma::mean(arma_mat, 1);
    return result;
}

/************************************************************************
 * Column and rows access.
 * The implementation of these functions depends on the data ordering: column-major or row-major
 */

/** Set column c of this matrix to the value of column vector a */
DoubleMatrix& DoubleMatrix::putColumn(int c, const DoubleColumn& a) {
    a.checkRows(getRows());
    arma_mat.col(c) = a.arma_mat;
    return *this;
}

/** Set row r of this matrix to the row vector a */
DoubleMatrix& DoubleMatrix::putRow(int r, const DoubleRow& a) {
    a.checkRows(1);
    a.checkColumns(getColumns());
    arma_mat.row(r) = a.arma_mat;
    return *this;
}

/** returns a row vector that contains the column-wise minimums. */
DoubleRow DoubleMatrix::columnMins() const {
    DoubleRow result(getColumns());
    result.arma_mat = arma::min(arma_mat, 0); // compute the min for dimension 0 (across columns)
    return result;
}

/**
 * This function receives a pointer to an array, and populates it with the index of
 * minimal element per column. Use this with caution to create the array before passing it.
 * */
void DoubleMatrix::columnArgmins(int* argmins) const {
    arma::uword index;
    for (int c = 0; c < getColumns(); c++) {
        arma::colvec col = arma_mat.col(c);
        index = col.index_min();
        argmins[c] = static_cast<int>(index);
    }
}

/** return a row vector that contains the column-wise maximums. */
DoubleRow DoubleMatrix::columnMaxs() const {
    DoubleRow result(getColumns());
    result.arma_mat = arma::max(arma_mat, 0);
    return result;
}

/**
 * This function receives a pointer to an array, and populates it with the index of
 * maximal element per column. Use this with caution. Array must be created properly
 * before passing it.
 * */
void DoubleMatrix::columnArgmaxs(int* argmaxs) const {
    arma::uword index;
    for (int c = 0; c < getColumns(); c++) {
        arma::colvec col = arma_mat.col(c);
        index = col.index_max();
        argmaxs[c] = static_cast<int>(index);
    }
}

/** return a column vector that contains the row-wise minimums. */
DoubleColumn DoubleMatrix::rowMins() const {
    DoubleColumn result(getRows());
    result.arma_mat = arma::min(arma_mat, 1);
    return result;
}

/**
 * This function receives a pointer to an array, and populates it with the index of
 * minimal element per row. Use this with caution to create the array before passing it.
 * */
void DoubleMatrix::rowArgmins(int* argmins) const {
    arma::uword index;
    for (int r = 0; r < getRows(); r++) {
        arma::rowvec row = arma_mat.row(r);
        index = row.index_min();
        argmins[r] = static_cast<int>(index);
    }
}

/** Set column vector that contains the row-wise maximums. */
DoubleColumn DoubleMatrix::rowMaxs() const {
    DoubleColumn result(getRows());
    result.arma_mat = arma::max(arma_mat, 1);
    return result;
}

/**
 * This function receives a pointer to an array, and populates it with the index of
 * maximal element per row. Use this with caution to create the array before passing it.
 * */
void DoubleMatrix::rowArgmaxs(int* argmaxs) const {
    arma::uword index;
    for (int r = 0; r < getRows(); r++) {
        arma::rowvec row = arma_mat.row(r);
        index = row.index_max();
        argmaxs[r] = static_cast<int>(index);
    }
}

/**************************************************************************
 * Elementwise Functions
 * The following function depend on data ordering column major or row major
 */

/** Add a row vector to all rows of the matrix (in place). */
DoubleMatrix& DoubleMatrix::addiRowVector(const DoubleRow& a) {
    a.checkRows(1);
    a.checkColumns(getColumns());
    arma_mat.each_row() += a.arma_mat;
    return *this;
}

/** Add a row vector to all rows of the matrix, and return result */
DoubleMatrix DoubleMatrix::addRowVector(const DoubleRow& a) const {
    DoubleMatrix result(getRows(), getColumns());
    result.arma_mat = arma_mat;
    result.addiRowVector(a);
    return result;
}

/** Add a column vector to all columns of the matrix (in-place). */
DoubleMatrix& DoubleMatrix::addiColumnVector(const DoubleColumn& a) {
    a.checkRows(getRows());
    a.checkColumns(1);
    arma_mat.each_col() += a.arma_mat;
    return *this;
}

/** Add a column vector to all columns of the matrix, and return the result */
DoubleMatrix DoubleMatrix::addColumnVector(const DoubleColumn& a) const {
    DoubleMatrix result(getRows(), getColumns());
    result.arma_mat = arma_mat;
    result.addiColumnVector(a);
    return result;
}

/** Subtract a row vector from all rows of the matrix (in place). */
DoubleMatrix& DoubleMatrix::subiRowVector(const DoubleRow& a) {
    a.checkRows(1);
    a.checkColumns(getColumns());
    arma_mat.each_row() -= a.arma_mat;
    return *this;
}

/** Subtract a row vector from all rows of the matrix, and return the result matrix */
DoubleMatrix DoubleMatrix::subRowVector(const DoubleRow& a) const {
    DoubleMatrix result(getRows(), getColumns());
    result.arma_mat = arma_mat;
    result.subiRowVector(a);
    return result;
}

/** Subtract a column vector from all columns of the matrix (in-place). */
DoubleMatrix& DoubleMatrix::subiColumnVector(const DoubleColumn& a) {
    a.checkRows(getRows());
    a.checkColumns(1);
    arma_mat.each_col() -= a.arma_mat;
    return *this;
}

/** Subtract a column vector from all columns of the matrix, and return the result */
DoubleMatrix DoubleMatrix::subColumnVector(const DoubleColumn& a) const {
    DoubleMatrix result(getRows(), getColumns());
    result.arma_mat = arma_mat;
    result.subiColumnVector(a);
    return result;
}

/** Multiply row r by scalar scale (in-place). */
DoubleMatrix& DoubleMatrix::muliRow(int r, double scale) {
    arma_mat.row(r) *= scale;
    return *this;
}

/** Multiply column c by a scalar scale (in-place). */
DoubleMatrix& DoubleMatrix::muliColumn(int c, double scale) {
    arma_mat.col(c) *= scale;
    return *this;
}

/** Elementwise-multiply all columns with a column vector (in-place). */
DoubleMatrix& DoubleMatrix::muliColumn(const DoubleColumn& a) {
    a.checkRows(getRows());
    a.checkColumns(1);
    arma_mat.each_col() %= a.arma_mat;
    return *this;
}

/** Elementwise-multiply all columns with a column vector, return the result */
DoubleMatrix DoubleMatrix::mulColumn(const DoubleColumn& a) const {
    DoubleMatrix result(getRows(), getColumns());
    result.arma_mat = arma_mat;
    result.muliColumn(a);
    return result;
}

/** Elementwise-multiply all rows with a row vector (in-place). */
DoubleMatrix& DoubleMatrix::muliRow(const DoubleRow& a) {
    a.checkRows(1);
    a.checkColumns(getColumns());
    arma_mat.each_row() %= a.arma_mat;
    return *this;
}

/** Elementwise-multiply all rows with a row vector, and return result */
DoubleMatrix DoubleMatrix::mulRow(const DoubleRow& a) const {
    DoubleMatrix result(getRows(), getColumns());
    result.arma_mat = arma_mat;
    result.muliRow(a);
    return result;
}

/** Divide all rows by a row vector (in-place). */
DoubleMatrix& DoubleMatrix::diviRow(const DoubleRow& a) {
    a.checkRows(1);
    a.checkColumns(getColumns());
    arma_mat.each_row() /= a.arma_mat;
    return *this;
}

/** Divide all rows by a row vector, and set result */
DoubleMatrix DoubleMatrix::divRow(const DoubleRow& a) const {
    DoubleMatrix result(getRows(), getColumns());
    result.arma_mat = arma_mat;
    result.diviRow(a);
    return result;
}

/** Divide all columns by a column vector (in-place). */
DoubleMatrix& DoubleMatrix::diviColumn(const DoubleColumn& a) {
    a.checkRows(getRows());
    a.checkColumns(1);
    arma_mat.each_col() /= a.arma_mat;
    return *this;
}

/** Divide all columns by a column vector, and return the result */
DoubleMatrix DoubleMatrix::divColumn(const DoubleColumn& a) const {
    DoubleMatrix result(getRows(), getColumns());
    result.arma_mat = arma_mat;
    result.diviColumn(a);
    return result;
}


/** Add matrix a and set result. */
DoubleMatrix DoubleMatrix::add(const DoubleMatrix& a) const {
	DoubleMatrix result(getRows(),getColumns());
	return result.put(*this).addi(a);
}

/** Element-wise add a scalar v and return the result. */
DoubleMatrix DoubleMatrix::add(double v) const {
	DoubleMatrix result(getRows(),getColumns());
	return result.put(*this).addi(v);
}

/** Subtract matrix a and return the result. */
DoubleMatrix DoubleMatrix::sub(const DoubleMatrix& a) const {
	DoubleMatrix result(getRows(),getColumns());
	return result.put(*this).subi(a);
}

/** Element-wise subtract scalar v and return the result. */
DoubleMatrix DoubleMatrix::sub(double v) const {
	DoubleMatrix result(getRows(),getColumns());
	return result.put(*this).subi(v);
}

/** Elementwise divide by a matrix and return the result */
DoubleMatrix DoubleMatrix::div(const DoubleMatrix& a) const {
	DoubleMatrix result(getRows(),getColumns());
	return result.put(*this).divi(a);
}

/** Elementwise divide by a scalar and return the result */
DoubleMatrix DoubleMatrix::div(double v) const {
	DoubleMatrix result(getRows(),getColumns());
	return result.put(*this).divi(v);
}

/** Elementwise multiply by a matrix and return the result. */
DoubleMatrix DoubleMatrix::mul(const DoubleMatrix& a) const {
	DoubleMatrix result(getRows(),getColumns());
	return result.put(*this).muli(a);
}

/** Elementwise multiply by a scalar and return the result */
DoubleMatrix DoubleMatrix::mul(double v) const {
	DoubleMatrix result(getRows(),getColumns());
	return result.put(*this).muli(v);
}

/** Test for "less than" (in-place). */
DoubleMatrix& DoubleMatrix::lti(const DoubleMatrix& a) {
	assertSameLength(a);
	for (int i = 0; i < getLength(); i++)
		put(i, get(i) < a.get(i) ? 1.0 : 0.0);
	return *this;
}

/** Test for "less than", and return the result */
DoubleMatrix DoubleMatrix::lt(const DoubleMatrix& a) const {
	DoubleMatrix result(getRows(),getColumns());
	return result.put(*this).lti(a);
}

/** Test for "less than" against a scalar (in-place). */
DoubleMatrix& DoubleMatrix::lti(double value) {
	for (int i = 0; i < getLength(); i++)
		put(i, get(i) < value ? 1.0 : 0.0);
	return *this;
}

/** test for "less than" against a scalar, and return the result */
DoubleMatrix DoubleMatrix::lt(double value) const {
	DoubleMatrix result(getRows(),getColumns());
	return result.put(*this).lti(value);
}

/** Test for "greater than" (in-place). */
DoubleMatrix& DoubleMatrix::gti(const DoubleMatrix& a) {
	assertSameLength(a);
	for (int i = 0; i < getLength(); i++)
		put(i, get(i) > a.get(i) ? 1.0 : 0.0);
	return *this;
}

/** Test for "greater than", and return the result */
DoubleMatrix DoubleMatrix::gt(const DoubleMatrix& a) const {
	DoubleMatrix result(getRows(),getColumns());
	return result.put(*this).gti(a);
}

/** Test for "greater than" against a scalar (in-place). */
DoubleMatrix& DoubleMatrix::gti(double value) {
	for (int i = 0; i < getLength(); i++)
		put(i, get(i) > value ? 1.0 : 0.0);
	return *this;
}

/** test for "greater than" against a scalar, and return the result */
DoubleMatrix DoubleMatrix::gt(double value) const {
	DoubleMatrix result(getRows(),getColumns());
	return result.put(*this).gti(value);
}

/** Test for "less than or equal" (in-place). */
DoubleMatrix& DoubleMatrix::lei(const DoubleMatrix& a) {
	assertSameLength(a);
	for (int i = 0; i < getLength(); i++)
		put(i, get(i) <= a.get(i) ? 1.0 : 0.0);
	return *this;
}

/** Test for "less than or equal", and return the result */
DoubleMatrix DoubleMatrix::le(const DoubleMatrix& a) const {
	DoubleMatrix result(getRows(),getColumns());
	return result.put(*this).lei(a);
}

/** Test for "less than or equal" against a scalar (in-place). */
DoubleMatrix& DoubleMatrix::lei(double value) {
	for (int i = 0; i < getLength(); i++)
		put(i, get(i) <= value ? 1.0 : 0.0);
	return *this;
}

/** Test for "less than or equal" against a scalar, and return result */
DoubleMatrix DoubleMatrix::le(double value) const {
	DoubleMatrix result(getRows(),getColumns());
	return result.put(*this).lei(value);
}

/** Test for "greater than or equal" (in-place). */
DoubleMatrix& DoubleMatrix::gei(const DoubleMatrix& a) {
	assertSameLength(a);
	for (int i = 0; i < getLength(); i++)
		put(i, get(i) >= a.get(i) ? 1.0 : 0.0);
	return *this;
}

/** Test for "greater than or equal", and return result. */
DoubleMatrix DoubleMatrix::ge(const DoubleMatrix& a) const {
	DoubleMatrix result(getRows(),getColumns());
	return result.put(*this).gei(a);
}

/** Test for "greater than or equal" against a scalar (in-place). */
DoubleMatrix& DoubleMatrix::gei(double value) {
	for (int i = 0; i < getLength(); i++)
		put(i, get(i) >= value ? 1.0 : 0.0);
	return *this;
}

/** Test for "greater than or equal" against a scalar, and return the result. */
DoubleMatrix DoubleMatrix::ge(double value) const {
	DoubleMatrix result(getRows(),getColumns());
	return result.put(*this).gei(value);
}

/** Test for equality (in-place). */
DoubleMatrix& DoubleMatrix::eqi(const DoubleMatrix& a) {
	assertSameLength(a);
	for (int i = 0; i < getLength(); i++)
		put(i, get(i) == a.get(i) ? 1.0 : 0.0);
	return *this;
}

/** Test for equality, and return the result */
DoubleMatrix DoubleMatrix::eq(const DoubleMatrix& a) const {
	DoubleMatrix result(getRows(),getColumns());
	return result.put(*this).eqi(a);
}

/** Test for equality against a scalar (in-place). */
DoubleMatrix& DoubleMatrix::eqi(double value) {
	for (int i = 0; i < getLength(); i++)
		put(i, get(i) == value ? 1.0 : 0.0);
	return *this;
}

/** test for equality against a scalar, and return the result */
DoubleMatrix DoubleMatrix::eq(double value) const {
	DoubleMatrix result(getRows(),getColumns());
	return result.put(*this).eqi(value);
}

/** Test for inequality (in-place). */
DoubleMatrix& DoubleMatrix::nei(const DoubleMatrix& a) {
	assertSameLength(a);
	for (int i = 0; i < getLength(); i++)
		put(i, get(i) != a.get(i) ? 1.0 : 0.0);
	return *this;
}

/** Test for inequality, and return the result */
DoubleMatrix DoubleMatrix::ne(const DoubleMatrix& a) const {
	DoubleMatrix result(getRows(),getColumns());
	return result.put(*this).nei(a);
}

/** Test for inequality against a scalar (in-place). */
DoubleMatrix& DoubleMatrix::nei(double value) {
	for (int i = 0; i < getLength(); i++)
		put(i, get(i) != value ? 1.0 : 0.0);
	return *this;
}

/** test for inequality against a scalar, and return result */
DoubleMatrix DoubleMatrix::ne(double value) const {
	DoubleMatrix result(getRows(),getColumns());
	return result.put(*this).nei(value);
}

/** Compute elementwise logical and (in-place). */
DoubleMatrix& DoubleMatrix::andi(const DoubleMatrix& a) {
	assertSameLength(a);
	for (int i = 0; i < getLength(); i++)
		put(i, ((get(i) != 0.0) && (a.get(i) != 0.0)) ? 1.0 : 0.0);
	return *this;
}

/** Compute elementwise logical and, and set to result */
DoubleMatrix DoubleMatrix::_and(const DoubleMatrix& a) const {
	DoubleMatrix result(getRows(),getColumns());
	return result.put(*this).andi(a);
}

/** Compute elementwise logical and against a scalar (in-place). */
DoubleMatrix& DoubleMatrix::andi(double value) {
	bool val = (value != 0.0);
	for (int i = 0; i < getLength(); i++)
		put(i, ((get(i) != 0.0) && (val != 0.0)) ? 1.0 : 0.0);
	return *this;
}

/** Compute elementwise logical and against a scalar, and return result. */
DoubleMatrix DoubleMatrix::_and(double value) const {
	DoubleMatrix result(getRows(),getColumns());
	return result.put(*this).andi(value);
}

/** Compute elementwise logical or (in-place). */
DoubleMatrix& DoubleMatrix::ori(const DoubleMatrix& a) {
	assertSameLength(a);
	for (int i = 0; i < getLength(); i++)
		put(i, ((get(i) != 0.0) || (a.get(i) != 0.0)) ? 1.0 : 0.0);
	return *this;
}

/** Compute elementwise logical or, and set to result */
DoubleMatrix DoubleMatrix::_or(const DoubleMatrix& a) const {
	DoubleMatrix result(getRows(),getColumns());
	return result.put(*this).ori(a);
}

/** Compute elementwise logical or against a scalar (in-place). */
DoubleMatrix& DoubleMatrix::ori(double value) {
	bool val = (value != 0.0);
	for (int i = 0; i < getLength(); i++)
		put(i, ((get(i) != 0.0) || (val != 0.0)) ? 1.0 : 0.0);
	return *this;
}

/** Compute elementwise logical or against a scalar, and return result */
DoubleMatrix DoubleMatrix::_or(double value) const {
	DoubleMatrix result(getRows(),getColumns());
	return result.put(*this).ori(value);
}

/** Compute elementwise logical xor (in-place). */
DoubleMatrix& DoubleMatrix::xori(const DoubleMatrix& a) {
	assertSameLength(a);
	for (int i = 0; i < getLength(); i++)
		put(i, (get(i) != 0.0) ^ (a.get(i) != 0.0) ? 1.0 : 0.0);
	return *this;
}

/** Compute elementwise logical xor, and return result */
DoubleMatrix DoubleMatrix::_xor(const DoubleMatrix& a) const {
	DoubleMatrix result(getRows(),getColumns());
	return result.put(*this).xori(a);
}

/** Compute elementwise logical xor against a scalar (in-place). */
DoubleMatrix& DoubleMatrix::xori(double value) {
	bool val = (value != 0.0);
	for (int i = 0; i < getLength(); i++)
		put(i, (get(i) != 0.0) ^ val ? 1.0 : 0.0);
	return *this;
}

/** Compute elementwise logical xor against a scalar, and return result */
DoubleMatrix DoubleMatrix::_xor(double value) const {
	DoubleMatrix result(getRows(),getColumns());
	return result.put(*this).xori(value);
}


//*************** LAPACK functions (solving linear equations) *******************//

/** Solves the linear equation (*this)*X = B.
 * X can be one or multiple columns and B has the same dimensions as X */
DoubleMatrix DoubleMatrix::solve(const DoubleMatrix& B) const {
	assertSquare();  // We restruct solver to square matrices only
	if (getRows() != B.getRows()) {  // if A . B = C, rows of A must be equal to rows of C
		throw("error: solving matrix equation is impossible");
	}
    DoubleMatrix S(getColumns(), B.getColumns());
    S.arma_mat = arma::solve(arma_mat, B.arma_mat);  // Solves A * X = B
	return S;
}

/** get the inverse matrix S directly or by solving the equation M S = I */
DoubleMatrix DoubleMatrix::getInverse() const {
	assertSquare();
    // Direct inversion using .i() or inv()
    DoubleMatrix result(getRows(), getColumns());
    result.arma_mat = arma::inv(arma_mat);  // or A.i()
    // // Another way
    // int n = getRows();
    // arma::mat I = arma::eye<arma::mat>(n, n);  // Identity matrix
    // result.arma_mat = arma::solve(arma_mat, I);  // Solves A * X = I
	return result;
}

/**
 * Assert that the matrix has a certain length.
 * @throws SizeException
 */
void DoubleMatrix::checkLength(int l) const {
	if (getLength() != l) {
		ostringstream oss;
		oss << "cppblas error: Matrix does not have the necessary length (" << getLength() << " != " << l << ").";
		string msg = oss.str();
		throw(msg);
	}
}

/**
 * Asserts that the matrix has a certain number of rows.
 * @throws SizeException
 */
void DoubleMatrix::checkRows(int r) const {
	if (getRows() != r) {
		ostringstream oss;
		oss << "cppblas error: Matrix does not have the necessary number of rows (" << getRows() << " != " << r << ").";
		string msg = oss.str();
		throw(msg);
	}
}

/**
 * Asserts that the amtrix has a certain number of columns.
 * @throws SizeException
 */
void DoubleMatrix::checkColumns(int c) const {
	if (getColumns() != c) {
		ostringstream oss;
		oss << "cppblas error: Matrix does not have the necessary number of columns (" << getColumns() << " != " << c << ").";
		string msg = oss.str();
		throw(msg);
	}
}

/** Checks whether two matrices have the same size. */
bool DoubleMatrix::sameDimensions(const DoubleMatrix& a) const {
	return getRows() == a.getRows() && getColumns() == a.getColumns();
}

/** Throws SizeException unless two matrices have the same dimensions. */
void DoubleMatrix::assertSameDimensions(const DoubleMatrix& a) const {
	if (!sameDimensions(a)) {
		throw ("cppblas error: Matrices must have the same dimensions .");
	}
}

/** Checks whether two matrices can be multiplied (that is, number of columns of
 * this must equal number of rows of a. */
bool DoubleMatrix::multipliesWith(const DoubleMatrix& a) const {
	return getColumns() == a.getRows();
}

/** Throws SizeException unless matrices can be multiplied with one another. */
void DoubleMatrix::assertMultipliesWith(const DoubleMatrix& a) const {
	if (!multipliesWith(a)) {
		throw ("cppblas error: Number of columns of left matrix must be equal to number of rows of right matrix.");
	}
}

/** Checks whether two matrices have the same length. */
bool DoubleMatrix::sameLength(const DoubleMatrix& a) const {
	return getLength() == a.getLength();
}

/** Throws SizeException unless matrices have the same length. */
void DoubleMatrix::assertSameLength(const DoubleMatrix& a) const {
	if (!sameLength(a)) {
		ostringstream oss;
		oss << "cppblas error: Matrices must have same length (is: " << getLength() << " and " << a.getLength() << ")";
		string msg = oss.str();
		throw(msg);
	}
}

/** Checks whether the matrix is empty. */
bool DoubleMatrix::isEmpty() const {
	return (getColumns() == 0 || getRows() == 0);
}

/** Checks whether the matrix is square. */
bool DoubleMatrix::isSquare() const {
	return (getColumns() == getRows());
}

/** Throw SizeException unless matrix is square. */
void DoubleMatrix::assertSquare() const {
	if (!isSquare()) {
		throw ("cppblas error: Matrix must be square!");
	}
}

/** Throw SizeException unless matrix is a vector. */
void DoubleMatrix::assertVector() const {
	if (!isVector()) {
		throw ("cppblas error: Matrix must be a vector!");
	}
}

/** Checks whether the matrix is a vector. */
bool DoubleMatrix::isVector() const {
	return (getColumns() == 1 || getRows() == 1);
}

/** Checks whether the matrix is a row vector. */
bool DoubleMatrix::isRowVector() const {
	return getRows() == 1;
}

/** Checks whether the matrix is a column vector. */
bool DoubleMatrix::isColumnVector() const {
	return getColumns() == 1;
}

/** Test whether a matrix is scalar. */
bool DoubleMatrix::isScalar() const {
	return getLength() == 1;
}

/** Return the first element of the matrix. */
double DoubleMatrix::scalar() const {
	return get(0);
}

/** get a string representation of the matrix.
 * p_precision (if >0) specifies the precision */
std::string DoubleMatrix::toString(int p_precision) const {
	std::ostringstream oss;
	if (p_precision > 0)
		oss << std::fixed << std::setprecision(p_precision);
	for (int r = 0; r<getRows(); r++) {
		for (int c=0;c<getColumns();c++)
			oss << get(r,c) << " ";
		oss << "\n";
	}
	string msg = oss.str();
	return msg;
}

/** write the matrix to a file. p_precision (if >0) specifies the precision */
void DoubleMatrix::write_to_file(std::string p_file_name, int p_precision) const {
	std::ofstream pfile(p_file_name);
	if (pfile.is_open()) {
		pfile << toString(p_precision);
		pfile.close();
	}
	else {
		std::cerr << "Unable to open file to save the matrix";
		exit(0);
	}

}


/** Returns the elements of the matrix in a std::vector of doubles
 * in row major order. This is useful for functions that requires
 * the data in vector format, e.g. in C++ random generators, the EMD
 * functions, and also if we need to sort the elements (e.g. in matrix
 * inversion estimation with projection invp).
 **/
std::vector<double> DoubleMatrix::toVector() const {
    std::vector<double> vec = arma::conv_to<std::vector<double>>::from(arma_mat.t());
	return vec;
}


// non-member functions
std::ostream& operator<<(std::ostream &strm, const DoubleMatrix &a) {
	return strm << a.toString();
}



