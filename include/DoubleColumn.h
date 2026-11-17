/*
 * DoubleColumn.h
 *
 *  Created on: Jul 8, 2019 - Updated on Mar 1, 2026 to use Armadillo
 *      
 */

#ifndef MAT_DOUBLECOLUMN_H_
#define MAT_DOUBLECOLUMN_H_

#include "DoubleMatrix.h"


// definition of CoulmnVector
class DoubleColumn : public DoubleMatrix {
public:
	/** Default constructor. Must be defined because we have other consrructors
	 * this calls the default constructor of the base class DoubleMatrix */
	DoubleColumn(){};

	/** Creates an vector of size n filled with 0. */
	DoubleColumn(int n);

	/** a conversion copy constructor from a matrix to a column. the given matrix must consist of a single column */
	DoubleColumn(const DoubleMatrix& other);

	/** a conversion move constructor: used by the compiler to convert (by move) rvalue DoubleMatrix (other) to this
	 * other must consist of a single column
	 * note: rvalue reference (&&) is a C++11 extension */
	DoubleColumn(DoubleMatrix&& other);

	/**
	 * Creates a column vector of size n filled with the given data_array.
	 * In the caller procedure use &a[0] to pass the array a[].
	 */
	DoubleColumn(int n, const double *data_array);

	/**
	 * The following operation assigns a matrix (if it is a column) to a column vector,
	 * and throws an error if the operand is not a proper column vector.
	 * This is important because the Linear Algebraic operations are defined for the super class DoubleMatrix.
	 * If the result of an operation is a column vector, it is assigned to a vector via the following function.
	 **/
//	DoubleColumn& operator=(DoubleMatrix other);

	/* move constructor that allows conversion from DoubleMatrix to DoubleColumn; */
	//DoubleColumn(const DoubleMatrix&& m);

	/** Resize the vector. All elements will be set to zero. */
	void resize(int n);

};




#endif /* MAT_DOUBLECOLUMN_H_ */
