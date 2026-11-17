/*
 * Observation.h
 *
 *  Created on: Sep 23, 2019
 *      
 */

#ifndef SRC_STAT_EMPIRICAL_OBSERVATION_H_
#define SRC_STAT_EMPIRICAL_OBSERVATION_H_

/*
 * This is an abstract Observation class that is instantiated to e.g. GridObservation and
 * BitVectorObservation.
 */

#include <string>

class Observation {
public:

	virtual ~Observation(){};

	/** get a string representation on the observation */
	virtual std::string toString() const=0;

	/**
	 *  The following function create a clone of this object and returns a pointer to the clone.
	 *  This is needed because pointers to Observation instances are stored as members in other
	 *  classes (EmpDatum). The constructors of those classes clone the passed pointer to
	 *  a new Observation instance, so that the new Observation object survives after the original
	 *  object is destroyed, and also to keep the original object alive when the EmpDatum is destructed
	 *  (note that the destructor of EmpDatum deletes the member pointers).
	 *  We declare it virtual to call the implementation of this function in the derived classes (e.g. GridObservation)
	 *  if exist. This is necessary if the derived class has additional members. Then the implementation of clone() in
	 *  that derived class copies these additional members. We make it pure virtual (using "const=0")
	 *  to force subclasses to implement it, then it is abstract class. Therefore its instances can be passed by reference.
	 *  not by value.
	 *
	 * replicate the observation and return a pointer to the replicate.
	 * This is invoked from EmpDatum constructor and the clone is deleted in EmpDatum destructor
	 * */
	virtual Observation* clone() const = 0;

	/** must be implemented by every subclass of Observation */
	virtual bool operator ==(const Observation& other) const =0;

	bool operator !=(const Observation& other) const; // need not to be virtual since it is implemented abstractly
};



#endif /* SRC_STAT_EMPIRICAL_OBSERVATION_H_ */
