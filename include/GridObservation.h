/*
 * GridObservation.h
 *
 *  Created on: Sep 23, 2019
 *      
 */

/*
 * This is a location in a grid represented by its linear order.
 */

#ifndef SRC_STAT_EMPIRICAL_GRIDOBSERVATION_H_
#define SRC_STAT_EMPIRICAL_GRIDOBSERVATION_H_

#include "Observation.h"

class GridObservation : public Observation {
public:

	int location; // location in a grid specified in the linear order (row major).

	/** constructor */
	GridObservation(int p_location);

	/** get a string representation on the observation */
	std::string toString() const;

	/** replicate the observation and return a pointer to the replicate */
	GridObservation* clone() const;

	/** returns true if the two observations are equal **/
	bool operator ==(const Observation& other) const;

};



#endif /* SRC_STAT_EMPIRICAL_GRIDOBSERVATION_H_ */
