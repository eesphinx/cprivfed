/*
 * BitVectorObservation.h
 *
 *  Created on: Sep 23, 2019
 *      
 */

#ifndef SRC_STAT_EMPIRICAL_BITVECTOROBSERVATION_H_
#define SRC_STAT_EMPIRICAL_BITVECTOROBSERVATION_H_

/*
 * This is a bit vector that results from obfuscating a location using a OneTimeRappor mechanism.
 */
#include <vector>
#include "Observation.h"

class BitVectorObservation : public Observation {
public:
	std::vector<bool> bitVector;  // observed bitVector

	// constructor
	BitVectorObservation(std::vector<bool> p_bitVector);

	/** get a string representation for the observation */
	std::string toString() const;

	/** replicate the observation and return a pointer to the replicate.
	 * This is called from EmpDatum constructor and the clone is deleted in the EmpDatum destructor */
	BitVectorObservation* clone() const;

	/** returns true if the two observations are equal **/
	bool operator ==(const Observation& other) const;

};

#endif /* SRC_STAT_EMPIRICAL_BITVECTOROBSERVATION_H_ */
