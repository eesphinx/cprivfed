/*
 * GridObservation.cpp
 *
 *  Created on: Sep 23, 2019
 *      Author: ehab
 */

#include <iostream>
#include "GridObservation.h"

/** constructor from a given location */
GridObservation::GridObservation(int p_location)
:location(p_location) {
}


/** get a string representation on the observation */
std::string GridObservation::toString() const {
	return std::to_string(location);
}

/** replicate the observation and return a pointer to the replicate */
GridObservation* GridObservation::clone() const {
	return new GridObservation(*this);
	// this uses the default copy constructor provided by c++ since
	// its members have also copy constructors (numeric types have)
}

/** returns true if the two observations are equal **/
bool GridObservation::operator ==(const Observation& other) const {
	// in the following we use const because other is declared as const
	if (const GridObservation* other_ptr = dynamic_cast<const GridObservation*>(&other))
		return (this->location == other_ptr->location);
	else {
		std::cerr << "unable to compare GridObservation objects";
		exit(0);
	}
}


