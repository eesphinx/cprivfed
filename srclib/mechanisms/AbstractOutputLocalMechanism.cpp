/*
 * AbstractOutputLocalMechanism.cpp
 *
 *  Created on: Sep 24, 2019
 *      
 */

#include "DoubleColumn.h"
#include "AbstractOutputLocalMechanism.h"


/** returns a mechanism column (out) as a column vector: executed only if its subclass does not
 * already have a mechanism matrix. If it has a matrix, it should override this method for performance.
 * p_out is specified in linear order (y*hsize + x) */
DoubleColumn AbstractOutputLocalMechanism::get_mech_column(const Observation* p_out_ptr) const {
	DoubleColumn retColVec(hsize*vsize);
	for (int in=0; in<hsize*vsize;in++)
		retColVec.put(in,getProb(in, p_out_ptr));
	return retColVec;
}




