/*
 *	
 *  DoubleMatrix.h
 *  Created on: Jun 11, 2019. Upgraded to use Armadillo on Mar 1, 2026
 *      
 */

#ifndef LIBMAT_DOUBLEMATRIX_H_
#define LIBMAT_DOUBLEMATRIX_H_

/**
 * Important Notes:
 *
 * 1. With the exception of constructors, any function that should return an object e.g. (DoubleMatrix::add)
 *    is implemented by creating an (empty) object in the caller, and passing it by reference (as classname& result) to the function.
 *
 * 2. Objects are created only by constructors. Static functions that returns an object (as in Java) are
 *    replaced by creating an empty object and modifying it by a normal function. For example to create
 *    a row vector Ones, we construct a matrix [m] of one row and then use m.fill(1.0).
 *
 * 3. Every function that modifies the content of the matrix (e.g. put, fill, etc) returns a reference to
 *    this object by returning [*this]. This allows chained evaluation statements as a.fill(b).putRow(r). etc.
 *    The in-place operations e.g. addi(other) we return result (*this) by reference instead of value. This
 *    does not produce out-of-scope problem because the returned object is the same one invoked from the caller.
 *    Then returning by reference in this case is safe and also more efficient.
 *    Note that if we want to write c = a.fill(b), the class of c (DoubleMatrix in this case) has to allow
 *    assignment (note that a.fill(b) is treated as an object of type DoubleMatrix). This is why this operation
 *    is allowed for primitive types.
 *
 * 4. Notes on references: a reference is just an alias (another name) for the object. the reference is its referent.
 *    DoubleMatrix a,x;           (allowed indeed)
 *    DoubleMatrix& b,c,d;        (not allowed: a reference must be initialized as an alias to some variable)
 *    DoubleMatrix& b = a;        (allowed: reference (alias) b is initialized to refer to the object a)
 *    x = b                       (allowed if the class of x (DoubleMatrix) allows assignments. If x is a primitive
 *                                 type, e.g. int or double, it is allowed. Recall that b an alias to the variable
 *                                 of the same type)
 *    DoubleMatrix& c = b;        (allowed and efficient: initialize the reference c to the object having reference b).
 *    DoubleMatrix& c = b; c = d; (allowed if DoubleMatrix allows assignment. Note this operation changes the underlying
 *                                 object (referenced by) c. The reference stays the same)
 *    int a,d; int& c = a; c = d; (allowed since c is reference to primitive data type.)
 *
 * 5. For array arguments/results we pass a pointer to the array together with the length.
 *
 * 6. Any function that returns array e.g. DoubleMatrix.columnArgmins is implemented as a void function
 *    that receive a pointer to an array, and then populates it.
 *
 * 7. Try to avoid: Any class that dynamically allocate data, e.g. DoubleMatrix must have a destructor
 *    to to delete these data
 *
 * 8. member functions that takes big object arguments, these arguments are passed by constant reference
 *    as [const &Class arg]. "const" is used to avoid changing <arg>.
 *
 * 9. functions that should not change the object internal data should be declared "const".
 *
 *
 */

#include <string>
#include <vector> // must be included because it is used as argument for some functions
#include <armadillo>


class DoubleRow;
class DoubleColumn;

class DoubleMatrix {

public:
	
    /** Number of rows. */
	// int rows;
	/** Number of columns */
	// int columns;
	/** The actual data stored by rows (that is, row 0, row 1...). */
	// double *data;  //rows are contiguous
	arma::mat arma_mat;

	/** constructors */

	DoubleMatrix();  // default constructor. Must be defined because we have other constructors.

	/**
	 * Creates a matrix of n rows and m columns.
	 */
	DoubleMatrix(int rows, int columns);

	/**
	 * Creates an n times m matrix filled from the given data source pointer.
	 * Data is read sequentially from the source. If the boolean parameter is
	 * true the obtained data are stored in the matrix row by row
	 * (used to read from C++ arrays). Otherwise, the obtained data are stored
	 * in the matrix filling column by column (used to integrate with libraries
	 * that stores data as column by column, e.g. Armadillo library)
	 *
	 * In the caller procedure use &a[0] if a is one-dimensional array,
	 * &a[0][0] if a is two-dimensional array.
	 */
	DoubleMatrix(int p_rows, int p_columns, const double* p_data_source, bool p_fill_row_by_row = true);

	/** copy constructor */
	DoubleMatrix(const DoubleMatrix& other);

	/** move constructor: used by the compiler to move rvalue DoubleMatrix (other) to this.
	 * Important: implementing this move constructor causes implicit deletion of the default
	 * copy assignment operator. Therefore we should implement it explicitly
	 * note: rvalue reference (&&) is a C++11 extension */
	DoubleMatrix(DoubleMatrix&& other);

	/** The assignment operator, used for assigning objects. Works for both copy and move assignments. */
	DoubleMatrix& operator=(DoubleMatrix other);

	/** fill the matrix with random values uniformly in 0..1. */
	DoubleMatrix& fillRand();

	/** fill the matrix with the given value */
	DoubleMatrix& fillWith(double v);

	/** make this an identity (eye) matrix. It must be square matrix*/
	DoubleMatrix& makeItEye();

	/**
	 * Make it a diagonal matrix of arbitrary shape and set the diagonal according
	 * to a passed vector d. length of d must be smaller than rows or columns.
	 */
	DoubleMatrix& makeItDiag(const DoubleMatrix& d);

	/** member methods */

	/**
	 * Concatenates this matrix to the given matrix b horizontally. Matrices must have identical
	 * numbers of rows.
	 */
	DoubleMatrix concatHorizontally(const DoubleMatrix& b) const;

	/**
	 * Concatenates this matrix to the given matrix b vertically. Matrices must have identical
	 * numbers of columns.
	 */
	DoubleMatrix concatVertically(const DoubleMatrix& b) const;


	/**************************************************************************
	 * restricting the matrix to specified rows	and columns
	 */

	/** Restrict the matrix in-place to the rows passed in the vector rowsIndexes */
	DoubleMatrix& selectiRows(const std::vector<int>& rowsIndexes);

	/** Restrict the matrix in-place to the columns passed in the vector rowsIndexes */
	DoubleMatrix& selectiColumns(const std::vector<int>& columnsIndexes);

	/** get row having index r of this matrix. */
	DoubleRow selectRow(int r) const;

	/** Restrict the matrix to the rows passed in the vector rowsIndexes*/
	DoubleMatrix selectRows(const std::vector<int>& rowsIndexes) const;

	/** Set the column c of this matrix. */
	DoubleColumn selectColumn(int c) const;

	/** Restrict the matrix to the columns passed in the vector rowsIndexes */
	DoubleMatrix selectColumns(const std::vector<int>& columnsIndexes) const;



	/**************************************************************************
	 * Basic operations (copying, resizing, element access)
	 */

	/** Get transposed copy of this matrix. */
	DoubleMatrix transpose() const;

	/** Resize the matrix. All elements will be set to zero. */
	DoubleMatrix& resize(int newRows, int newColumns);

	/** Reshape the matrix. Number of elements must not change. */
	DoubleMatrix& reshape(int newRows, int newColumns);

	/** Swap two columns of a matrix. */
	DoubleMatrix& swapColumns(int i, int j);

	/** Swap two rows of a matrix. */
	DoubleMatrix& swapRows(int i, int j);

	/** Set matrix element */
	DoubleMatrix& put(int rowIndex, int columnIndex, double value);

	/** Retrieve matrix element */
	double get(int rowIndex, int columnIndex) const;

	/** Get the linear index of an element */
	int index(int rowIndex, int columnIndex) const;

	/** Compute the row index of a linear index. */
	int rowIndex(int i) const;

	/** Compute the column index of a linear index. */
	int columnIndex(int i) const;

	/** Get a matrix element (linear indexing). */
	double get(int i) const;

	/** Set a matrix element (linear indexing). */
	DoubleMatrix& put(int i, double v);

	/** Get number of rows. */
	int getRows() const;

	/** Get number of columns. */
	int getColumns() const;

	/** Get total number of elements. */
	int getLength() const;

	/** return the diagonal of this matrix as a column vector */
	DoubleColumn getDiag() const;

	/** Add matrix a to this (in-place). */
	DoubleMatrix& addi(const DoubleMatrix& a);

	/** Add scalar v (element-wise) to a matrix (in-place). */
	DoubleMatrix& addi(double v);

	/** Subtract matrix a from this (in-place). */
	DoubleMatrix& subi(const DoubleMatrix& a);

	/** Subtract scalar v (element-wise) from this (in-place). */
	DoubleMatrix& subi(double v);

	/** Elementwise multiplication (in-place): multiply this by a and return this*/
	DoubleMatrix& muli(const DoubleMatrix& a);

	/** Elementwise multiplication with a scalar (in-place): multiply this by a and return this. */
	DoubleMatrix& muli(double v);

	/** Matrix-matrix multiplication (in-place): this = this.a */
	DoubleMatrix& mmuli(const DoubleMatrix& a);

	/** Matrix-multiply by a matrix: result = this * a   */
	DoubleMatrix mmul(const DoubleMatrix& a) const;

	/** Elementwise division (in-place). */
	DoubleMatrix& divi(const DoubleMatrix& a);

	/** Elementwise division with a scalar (in-place). */
	DoubleMatrix& divi(double a);

	/** Negate each element (in-place). */
	DoubleMatrix& negi();

	/** Negate each element, and set the result to matrix result */
	DoubleMatrix neg() const;

	/** Maps zero to 1.0 and all non-zero values to 0.0 (in-place). */
	DoubleMatrix& noti();

	/** Maps zero to 1.0 and all non-zero values to 0.0. Set the result matrix*/
	DoubleMatrix nott() const;

	/** Maps zero to 0.0 and all non-zero values to 1.0 (in-place). */
	DoubleMatrix& truthi();

	/** Maps zero to 0.0 and all non-zero values to 1.0. Set the result matrix*/
	DoubleMatrix truth() const;

	/** Returns the minimal element of the matrix. */
	double min() const;

	/**
	 * Returns the linear index of the minimal element. If there are
	 * more than one elements with this value, the first one is returned.
	 * if the matrix is empty or have all NaN numbers returns -1
	 */
	int argmin() const;

	/**
	 * Computes the minimum between this matrix and the matrix a. Returns the smaller of the
	 * corresponding elements in `this' matrix (in-place).
	 */
	DoubleMatrix& mini(const DoubleMatrix& a);

	/**
	 * Computes the minimum between 'this' matrix and matrix a. Set result to
	 * the smaller of the corresponding elements in a new matrix.
	 */
	DoubleMatrix min(const DoubleMatrix& a) const;

	/**
	 * Computes the minimum between every element of this matrix and the value v.
	 * Returns the smaller of the corresponding elements (in-place).
	 */
	DoubleMatrix& mini(double v);

	/**
	 * Computes the minimum between every element of this matrix and the value v.
	 * Returns the smaller of the corresponding elements in result.
	 */
	DoubleMatrix min(double v) const;

	/** Returns the maximal element of the matrix. */
	double max() const;

	/**
	 * Returns the linear index of the maximal element. If there are
	 * more than one elements with this value, the first one is returned.
	 * if the matrix is empty or have all NaN numbers returns -1
	 */
	int argmax() const;

	/**
	 * Computes the maximum between this matrix and the matrix a. Returns the smaller of the
	 * corresponding elements (in-place).
	 */
	DoubleMatrix& maxi(const DoubleMatrix& a);

	/**
	 * Computes the maximum between 'this' matrix and matrix a. Set the matrix result
	 * containing the maximum of the corresponding elements.
	 */
	DoubleMatrix max(const DoubleMatrix& a) const;

	/**
	 * Computes the maximum between every element of this matrix and the value v.
	 * Returns the smaller of the corresponding elements (in-place).
	 */
	DoubleMatrix& maxi(double v);

	/**
	 * Computes the maximum between every element of this matrix and the value v.
	 * Set the matrix results to the bigger of the corresponding elements.
	 */
	DoubleMatrix max(double v) const;

	/** Computes the sum of all elements of the matrix. */
	double sum() const;

	/** Computes the product of all elements of the matrix */
	double prod() const;

	/**
	 * Computes the mean value of all elements in the matrix,
	 * that is, <code>x.sum() / x.length</code>.
	 */
	double mean() const;

	/**
	 * Computes the cumulative sum, that is, the sum of all elements
	 * of the matrix up to a given index in linear addressing (in-place).
	 */
	DoubleMatrix& cumulativeSumi();

	/**
	 * Computes the cumulative sum, that is, the element at index i contains the sum of all elements
	 * up to i in linear addressing. Set the result matrix.
	 */
	DoubleMatrix cumulativeSum() const;

	/** Returns scalar product of this with matrix a: i.e. the sum of element-wise products */
	double dot(const DoubleMatrix& a) const;

	/**
	 * Computes the projection coefficient of other on this.
	 *
	 * The returned scalar times <tt>this</tt> is the orthogonal projection
	 * of <tt>other</tt> on <tt>this</tt>.
	 */
	double project(const DoubleMatrix& a) const;

	/**
	 * The Euclidean norm of the matrix as vector, also the Frobenius
	 * norm of the matrix.
	 */
	double norm2() const;

	/**
	 * The maximum norm of the matrix (maximal absolute value of the elements).
	 */
	double normmax() const;

	/**
	 * The 1-norm of the matrix (sum of absolute values of elements).
	 */
	double norm1() const;

	/**
	 * Returns the squared distance between this and a.
	 */
	double squaredDistance(const DoubleMatrix& a) const;

	/**
	 * Returns the (Euclidean) distance.
	 */
	double distance2(const DoubleMatrix& a) const;

	/**
	 * Returns the (1-norm) distance between this and a.
	 */
	double distance1(const DoubleMatrix& a) const;

	/** Return a row vector containing the sum of every column */
	DoubleRow columnSums() const;

	/** returns a row vector that contains the mean of every column. */
	DoubleRow columnMeans() const;

	/** returns a column vector that contains the sum of every row. */
	DoubleColumn rowSums() const;

	/** returns a column vector that contains the mean of every row. */
	DoubleColumn rowMeans() const;

	/** Set column c of this matrix to the column vector a */
	DoubleMatrix& putColumn(int c, const DoubleColumn& a);

	/** Set row r of this matrix to the row vector a */
	DoubleMatrix& putRow(int r, const DoubleRow& a);

	/** return a row vector that contains the column-wise minimums. */
	DoubleRow columnMins() const;

	/**
	 * This function receives a pointer to an array, and populates it with the index of
	 * minimal element per column. Use this with caution to create the array before passing it.
	 * */
	void columnArgmins(int* argmins) const;

	/** returns a row vector that contains the column-wise maximums. */
	DoubleRow columnMaxs() const;

	/**
	 * This function receives a pointer to an array, and populates it with the index of
	 * maximal element per column. Use this with caution. Array must be created properly
	 * before passing it.
	 * */
	void columnArgmaxs(int* argmaxs) const;

	/** returns a column vector that contains the row-wise minimums. */
	DoubleColumn rowMins() const;

	/**
	 * This function receives a pointer to an array, and populates it with the index of
	 * minimal element per row. Use this with caution to create the array before passing it.
	 * */
	void rowArgmins(int* argmins) const;

	/** returns a column vector that contains the row-wise maximums. */
	DoubleColumn rowMaxs() const;

	/**
	 * This function receives a pointer to an array, and populates it with the index of
	 * maximal element per row. Use this with caution to create the array before passing it.
	 * */
	void rowArgmaxs(int* argmaxs) const;


	/**************************************************************************
	 * Elementwise Functions
	 */

	/** Add a row vector to all rows of the matrix (in place). */
	DoubleMatrix& addiRowVector(const DoubleRow& a);

	/** Add a row vector to all rows of the matrix, and return the result */
	DoubleMatrix addRowVector(const DoubleRow& a) const;

	/** Add a column vector to all columns of the matrix (in-place). */
	DoubleMatrix& addiColumnVector(const DoubleColumn& a);

	/** Add a column vector to all columns of the matrix, and return the result */
	DoubleMatrix addColumnVector(const DoubleColumn& a) const;

	/** Subtract a row vector a from all rows of the matrix (in place). */
	DoubleMatrix& subiRowVector(const DoubleRow& a);

	/** Subtract a row vector from all rows of the matrix, and return the result */
	DoubleMatrix subRowVector(const DoubleRow& a) const;

	/** Subtract a vector from all columns of the matrix (in-place). */
	DoubleMatrix& subiColumnVector(const DoubleColumn& a);

	/** Subtract a column vector from all columns of the matrix, and return the result */
	DoubleMatrix subColumnVector(const DoubleColumn& a) const;

	/** Multiply row r by scalar scale (in-place). */
	DoubleMatrix& muliRow(int r, double scale);

	/** Multiply column c by a scalar scale (in-place). */
	DoubleMatrix& muliColumn(int c, double scale);

	/** Elementwise-multiply all columns with a column vector a (in-place). */
	DoubleMatrix& muliColumn(const DoubleColumn& a);

	/** Elementwise-multiply all columns with a vector, and return the result */
	DoubleMatrix mulColumn(const DoubleColumn& a) const;

	/** Elementwise-multiply all rows with a row vector (in-place). */
	DoubleMatrix& muliRow(const DoubleRow& a);

	/** Elementwise-multiply all rows with a row vector, and return the result */
	DoubleMatrix mulRow(const DoubleRow& a) const;

	/** Divide all rows by a row vector (in-place). */
	DoubleMatrix& diviRow(const DoubleRow& a);

	/** Divide all rows by a row vector, and return the result */
	DoubleMatrix divRow(const DoubleRow& a) const;

	/** Divide all columns by a column vector (in-place). */
	DoubleMatrix& diviColumn(const DoubleColumn& a);

	/** Divide all columns by a column vector, and return the result */
	DoubleMatrix divColumn(const DoubleColumn& a) const;

	/** Add matrix a and return the result. */
	DoubleMatrix add(const DoubleMatrix& a) const;

	/** Add a scalar v and return the result. */
	DoubleMatrix add(double v) const;

	/** Subtract matrix a and return the result. */
	DoubleMatrix sub(const DoubleMatrix& a) const;

	/** Subtract scalar v and return the result. */
	DoubleMatrix sub(double v) const;

	/** Elementwise divide by a matrix and return the result */
	DoubleMatrix div(const DoubleMatrix& a) const;

	/** Elementwise divide by a scalar and set to result */
	DoubleMatrix div(double v) const;

	/** Elementwise multiply by a matrix and return the result. */
	DoubleMatrix mul(const DoubleMatrix& a) const;

	/** Elementwise multiply by a scalar and return the result */
	DoubleMatrix mul(double v) const;

	/** Test for "less than" (in-place). */
	DoubleMatrix& lti(const DoubleMatrix& a);

	/** Test for "less than", and return the result */
	DoubleMatrix lt(const DoubleMatrix& a) const;

	/** Test for "less than" against a scalar (in-place). */
	DoubleMatrix& lti(double value);

	/** test for "less than" against a scalar, and return the result */
	DoubleMatrix lt(double value) const;

	/** Test for "greater than" (in-place). */
	DoubleMatrix& gti(const DoubleMatrix& a);

	/** Test for "greater than", and return the result */
	DoubleMatrix gt(const DoubleMatrix& a) const;

	/** Test for "greater than" against a scalar (in-place). */
	DoubleMatrix& gti(double value);

	/** test for "greater than" against a scalar, and return the result */
	DoubleMatrix gt(double value) const;

	/** Test for "less than or equal" (in-place). */
	DoubleMatrix& lei(const DoubleMatrix& a);

	/** Test for "less than or equal", and return the result */
	DoubleMatrix le(const DoubleMatrix& a) const;

	/** Test for "less than or equal" against a scalar (in-place). */
	DoubleMatrix& lei(double value);

	/** Test for "less than or equal" against a scalar, and return the result */
	DoubleMatrix le(double value) const;

	/** Test for "greater than or equal" (in-place). */
	DoubleMatrix& gei(const DoubleMatrix& a);

	/** Test for "greater than or equal", and return the result. */
	DoubleMatrix ge(const DoubleMatrix& a) const;

	/** Test for "greater than or equal" against a scalar (in-place). */
	DoubleMatrix& gei(double value);

	/** Test for "greater than or equal" against a scalar, and return the result. */
	DoubleMatrix ge(double value) const;

	/** Test for equality (in-place). */
	DoubleMatrix& eqi(const DoubleMatrix& a);

	/** Test for equality, and return the result */
	DoubleMatrix eq(const DoubleMatrix& a) const;

	/** Test for equality against a scalar (in-place). */
	DoubleMatrix& eqi(double value);

	/** test for equality against a scalar, and return the result */
	DoubleMatrix eq(double value) const;

	/** Test for inequality (in-place). */
	DoubleMatrix& nei(const DoubleMatrix& a);

	/** Test for inequality, and return the result */
	DoubleMatrix ne(const DoubleMatrix& a) const;

	/** Test for inequality against a scalar (in-place). */
	DoubleMatrix& nei(double value);

	/** test for inequality against a scalar, and return the result */
	DoubleMatrix ne(double value) const;

	/** Compute elementwise logical and (in-place). */
	DoubleMatrix& andi(const DoubleMatrix& a);

	/** Compute elementwise logical and, and return the result */
	DoubleMatrix _and(const DoubleMatrix& a) const;

	/** Compute elementwise logical and against a scalar (in-place). */
	DoubleMatrix& andi(double value);

	/** Compute elementwise logical and against a scalar, and return the result. */
	DoubleMatrix _and(double value) const;

	/** Compute elementwise logical or (in-place). */
	DoubleMatrix& ori(const DoubleMatrix& a);

	/** Compute elementwise logical or, and return the result */
	DoubleMatrix _or(const DoubleMatrix& a) const;

	/** Compute elementwise logical or against a scalar (in-place). */
	DoubleMatrix& ori(double value);

	/** Compute elementwise logical or against a scalar, and return the result */
	DoubleMatrix _or(double value) const;

	/** Compute elementwise logical xor (in-place). */
	DoubleMatrix& xori(const DoubleMatrix& a);

	/** Compute elementwise logical xor, and return the result */
	DoubleMatrix _xor(const DoubleMatrix& a) const;

	/** Compute elementwise logical xor against a scalar (in-place). */
	DoubleMatrix& xori(double value);

	/** Compute elementwise logical xor against a scalar, and return the result */
	DoubleMatrix _xor(double value) const;

	//*************** LAPACK functions (solving linear equations) *******************//


	// /** Solves the linear equation (*this)*X = B.
	//  * X can be one or multiple columns and B has the same dimensions as X*/
	DoubleMatrix solve(const DoubleMatrix& B) const;

	/** get the inverse matrix M which satisfies this.M = identity matrux */
	DoubleMatrix getInverse() const;

	/**
	 * Assert that the matrix has a certain length.
	 * @throws SizeException
	 */
	void checkLength(int l) const;

	/**
	 * Asserts that the matrix has a certain number of rows.
	 * @throws SizeException
	 */
	void checkRows(int r) const;

	/**
	 * Asserts that the amtrix has a certain number of columns.
	 * @throws SizeException
	 */
	void checkColumns(int c) const;

	/** Checks whether two matrices have the same dimensions. */
	bool sameDimensions(const DoubleMatrix& a) const;

	/** Throws SizeException unless two matrices have the same dimensions. */
	void assertSameDimensions(const DoubleMatrix& a) const;

	/** Checks whether two matrices can be multiplied (that is, number of columns of
	 * this must equal number of rows of a. */
	bool multipliesWith(const DoubleMatrix& a) const;

	/** Throws SizeException unless matrices can be multiplied with one another. */
	void assertMultipliesWith(const DoubleMatrix& a) const;

	/** Checks whether two matrices have the same length. */
	bool sameLength(const DoubleMatrix& a) const;

	/** Throws SizeException unless matrices have the same length. */
	void assertSameLength(const DoubleMatrix& a) const;

	/** Checks whether the matrix is empty. */
	bool isEmpty() const;

	/** Checks whether the matrix is square. */
	bool isSquare() const;

	/** Throw SizeException unless matrix is square. */
	void assertSquare() const;

	/** Throw SizeException unless matrix is a vector. */
	void assertVector() const;

	/** Checks whether the matrix is a vector. */
	bool isVector() const;

	/** Checks whether the matrix is a row vector. */
	bool isRowVector() const;

	/** Checks whether the matrix is a column vector. */
	bool isColumnVector() const;

	/** Test whether a matrix is scalar. */
	bool isScalar() const;

	/** Return the first element of the matrix. */
	double scalar() const;

	/** get a string representation of the matrix. p_precision (if >0) specifies the precision */
	std::string toString(int p_precision=0) const;

	/** write the matrix to a file. p_precision (if >0) specifies the precision */
	void write_to_file(std::string p_file_name, int p_precision=0) const;

	/**
	 * Returns the elements of the matrix in a std::vector of doubles
	 * in row major order. Used if we require to sort the elements.
	 */
	std::vector<double> toVector() const;


	// Important:
	// 1. this method is NOT a member of the class DoubleMatrix. But we declare it
	//    inside the class to give it access to the members of DoubleMatrix. We declare it
	//    as [friend] function to give it access to the private members. Otherwise it would
	//    have access only to the public members.
	// 2. Note also that [operator<<] is a non-member function, i.e. does not belong to any
	//    other class. Therefore we don't write its implementation as ClassName::operator<<.
	// 3. Since the argument a is passed as read-only [const], we can only call its members
	//    that are declared [const], i.e. which do not change the state of the object. Here
	//    the body of operator<< calls a.get(r,c) which is declared as [const]. Also note
	//    that [const] functions, e.g. get(r,c) can call only const functions, and const
	//    functions can be called from any function (const or not).
	// 4. It is a good practice to declare [const] when possible to minimize errors.
	friend std::ostream& operator<<(std::ostream& strm, const DoubleMatrix& a);

private:
	/** Put (copy) matrix a into this. The two matrices must have the same dimensions.
	 *  Note that it is in-place operation
	 */
	DoubleMatrix& put(const DoubleMatrix& a);

};  // class definition


#endif /* LIBMAT_DOUBLEMATRIX_H_ */
