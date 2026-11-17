/*
 * DoubleRow.cpp
 *
 *  Created on: Jun 27, 2019 - Updated on Mar 1, 2026 to use Armadillo
 *      Author: ehab
 */
#include <iostream>
#include "DoubleRow.h"

// Necessary because we there are other constructors. 
DoubleRow::DoubleRow(){
	// std::cout << "DoubleRow default constructor\n";
}
/**
 * Creates an vector of size n filled with 0.
 */
DoubleRow::DoubleRow(int n)
:DoubleMatrix(1,n){
// std::cout << "DoubleRow base constructor with n" <<"\n";
}

// The copy and move constructors are handled by the super class "DoubleMatrix"

// The following are conversion copy/move constructors (i.e., construct a RowVector from DoubleMatrix)
// These are not implemented by the super class, hence are implemented in the following.

/** a conversion copy constructor from a matrix to a row. the given matrix must consist of a single row */
DoubleRow::DoubleRow(const DoubleMatrix& other)
:DoubleMatrix(other.getRows(), other.getColumns()) {
	//todo: debug
	// std::cout << "conversion copy constructor from matrix to row\n";
	if (getRows() > 1) {
		std::cerr << "unable to convert a multi-rows matrix into DoubleRow.";
		exit(1);
	}
	arma_mat = other.arma_mat;
}

/** a conversion move constructor: used by the compiler to convert (by move) 
 * rvalue DoubleMatrix (other) to this */
DoubleRow::DoubleRow(DoubleMatrix&& other)
: DoubleMatrix(other.getRows(), other.getColumns()) {
	//todo: debug
    // std::cout << "conversion move constructor from matrix to row\n";
	if (getRows() > 1) {
		std::cerr << "unable to convert a multi-rows matrix into DoubleRow.";
		exit(0);
	}
	arma_mat = other.arma_mat;
	//std::cout << "matrix constructed: " << this->toString() << "\n";
}


/**
 * Creates a row vector of size n filled with the given data_array.
 * In the caller procedure use &a[0] to pass the array a[].
 */
DoubleRow::DoubleRow(int n, const double *data_array)
:DoubleMatrix(1,n,data_array) {
}


/** Resize the vector. All elements will be set to zero. */
void DoubleRow::resize(int n) {   // note it only shadows (not override) the parent member resize because they have different arguments
DoubleMatrix::resize(1, n);
}

/** Since there is no assignement operator defined for DoubleRow, the 
 * combiler implicitely creates default copy and move assignement operators 
 * of the form DoubleRow& DoubleRow::operator=(DoubleRow other), and they 
 * just invoke the copy/move assignements of the super class.  
 */

/** important */

/** We may think that we need a specific assignement operator for DoubleRow  
 * to ensure that only Rows are assigned to a DoubleRow object. However, we 
 * don't need that because the combiler creates a default assignement operation 
 * having the signature: 
 * 		DoubleRow& DoubleRow::operator=(DoubleRow other)
 * Note that the argement is (DoubleRow). When the combiler encounters 
 * "DoubleRow v = some_matrix", it first tries to convert some_matrix into 
 * DoubleRow, in order to allow assignement. Since there is a (copy/move) 
 * conversion constructors Matrix -> Row, one of them is invoked (depending on 
 * whether some_matrix is lvarue/rvalue). Inside these constructors, the check 
 * on the dimension is performed. In particular, it allows only converting 
 * one-row matrices, hence also allowing assigning these one-row matrices.  
*/

