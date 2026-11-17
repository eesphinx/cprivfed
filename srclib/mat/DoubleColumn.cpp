/*
 * DoubleColumn.cpp
 *
 *  Created on: Jul 8, 2019 - Updated on Mar 1, 2026 to use Armadillo
 *      Author: ehab
 */
#include <iostream>
#include "DoubleColumn.h"

/**
 * Creates an vector of size n filled with 0.
 */
DoubleColumn::DoubleColumn(int n)
:DoubleMatrix(n,1){
}

/** a conversion copy constructor from a matrix to a column. the given matrix must consist of a single column */
DoubleColumn::DoubleColumn(const DoubleMatrix& other)
:DoubleMatrix(other.getRows(), other.getColumns()) {
	//todo: debug
	//std::cout << "conversion copy constructor from matrix to column\n";
	if (getColumns() > 1) {
		std::cerr << "unable to convert a multi-columns matrix into DoubleColumn.";
		exit(0);
	}
	arma_mat = other.arma_mat;
}

/** a conversion move constructor: used by the compiler to convert (by move) rvalue DoubleMatrix (other) to a column */
DoubleColumn::DoubleColumn(DoubleMatrix&& other)
: DoubleMatrix(other.getRows(), other.getColumns()) {
	//todo: debug
	//std::cout << "conversion move constructor from matrix to a column\n";
	if (getColumns() > 1) {
		std::cerr << "unable to convert a multi-columns matrix into DoubleColumn.";
		exit(0);
	}
	arma_mat = other.arma_mat;
	//std::cout << "matrix constructed: " << this->toString() << "\n";
}

/**
 * Creates a column vector of size n filled with the given data_array.
 * In the caller procedure use &a[0] to pass the array a[].
 */
DoubleColumn::DoubleColumn(int n, const double *data_array)
:DoubleMatrix(n,1,data_array) {
}

/**
 * The following operation assigns a matrix (if it is a column) to a column vector,
 * and throws an error if the operand is not a proper column vector.
 * This is important because the Linear Algebraic operations are defined for the super class DoubleMatrix.
 * If the result of an operation is a column vector, it is assigned to a vector via the following function.
 **/
//DoubleColumn& DoubleColumn::operator=(DoubleMatrix other) {
//	if (!other.isColumnVector()) {
//		std::cerr << "error: unable to assign multiple columns to a column vector \n";
//		exit(-1);
//	}
//	std::cout << "column assignment operation \n";
//	rows = other.rows;
//	columns = other.columns;
//	std::swap(data, other.data);   // safe swap of the value of pointers data and other.data.
//	return *this;
//}


/** Resize the vector. All elements will be set to zero. */
void DoubleColumn::resize(int n) {
	DoubleMatrix::resize(n, 1);
}



