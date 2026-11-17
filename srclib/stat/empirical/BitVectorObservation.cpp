/*
 * BitVectorObservation.cpp
 *
 *  Created on: Sep 23, 2019
 *      Author: ehab
 */

#include <iostream>
#include "BitVectorObservation.h"

// constructor
BitVectorObservation::BitVectorObservation(std::vector<bool> p_bitVector)
:bitVector(p_bitVector) {  // assumes a copy constructor exists for vector<bool>
}


/** get a string representation on the observation */
std::string BitVectorObservation::toString() const {
	std::string ret= "";
	for (size_t i=0; i<bitVector.size(); i++)
		if (bitVector[i]) ret=ret + "1 "; else ret = ret+ "0 ";
	return ret;
}

/** replicate the observation and return a pointer to the replicate */
// This uses the default copy constructor assuming it exists for vector<boolean>
BitVectorObservation* BitVectorObservation::clone() const {
	return new BitVectorObservation(*this);
}

/** returns true if the two observations are equal **/
bool BitVectorObservation::operator ==(const Observation& other) const {
	// in the following we use const because other is declared as const
	if (const BitVectorObservation* other_ptr = dynamic_cast<const BitVectorObservation*>(&other))
		return (this->bitVector == other_ptr->bitVector);
	else {
		std::cerr << "unable to compare BitVectorObservation objects";
		exit(0);
	}
}


