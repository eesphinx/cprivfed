/*
 * DoubleRow.h
 *
 *  Created on: Jun 27, 2019
 *      
 */

#ifndef MAT_DOUBLEROW_H_
#define MAT_DOUBLEROW_H_

#include "DoubleMatrix.h"   // we include (and not forward declare) this because it is a super class


// definition of RowVector
class DoubleRow : public DoubleMatrix {
public:

	// /** Default constructor: this calls the default constructor of the base class DoubleMatrix */
	DoubleRow();  // default constructor. Must be defined because we have other constructors.

	/** Creates an row vector of size n filled with 0 */
	DoubleRow(int n);

	/** a conversion copy constructor from a matrix to a row. the given matrix must consist of a single row */
	DoubleRow(const DoubleMatrix& other);

	/** a conversion move constructor: used by the compiler to convert (by move) rvalue DoubleMatrix (other) to this
	 * other must consist of a single row
	 * note: rvalue reference (&&) is a C++11 extension */
	DoubleRow(DoubleMatrix&& other);

	/**
	 * Creates a row vector of size n filled with the given data_array.
	 * In the caller procedure use &a[0] to pass the array a[].
	 */
	DoubleRow(int n, const double *data_array);

	/**
	 * The following operation assigns a matrix (if it is a row) to a row vector,
	 * and throws an error if the operand is not a proper row vector.
	 * This is important because the Linear Algebraic operations are defined for the super class DoubleMatrix.
	 * If the result of an operation is a row vector, it is assigned to a vector via the following function.
	 **/
	// DoubleRow& operator=(DoubleMatrix other);
	// note the above overrides because the argument and return var are instances of DoubleMatrix as in the superclass


	/* a special copy constructor that is invoked explicitly or implicilitly convert a DoubleMatrix to DoubleRow; */
	//DoubleRow(const DoubleMatrix& m);

	/* a special move constructor that is invoked explicitly or implicilitly convert a DoubleMatrix to DoubleRow; */
	//DoubleRow(const DoubleMatrix&& m);

	/** Resize the vector. All elements will be set to zero. */
	void resize(int n);  // note it only shadows (not override) the parent member resize because they have different arguments

};





#endif /* MAT_DOUBLEROW_H_ */
