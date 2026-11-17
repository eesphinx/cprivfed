/*
 * EmpDatum.h
 *
 *  Created on: Sep 23, 2019
 *      
 */

#ifndef SRC_STAT_EMPIRICAL_EMPDATUM_H_
#define SRC_STAT_EMPIRICAL_EMPDATUM_H_


/**
 * a record that includes an observation, and the frequency of seeing this observation
 */

//#include <string>

// we may need a copy and move constructors to allow relocating objects when they are items in a vector,
// in particular the EmpData vector.

class Observation;
class AbstractOutputLocalMechanism;
class EmpDatum {
public:
	Observation* obs_ptr;  // an observation
	int frequency;         // number of times this observation is seen

	/** construct observation-frequency datum from an abstract Observation (which is cloned) and frequency */
	EmpDatum(Observation* p_obs, int p_frequency);

	/** construct a datum from a given data line in the format: Observation | frequency
	 * Note that the mechanism p_mech_ptr is necessary to parse the observation */
	EmpDatum(std::string p_data_line, const AbstractOutputLocalMechanism* p_local_mech_ptr);

	/** The following copy and move constructors are implemented since this class has an abstract (pointer) member obs_ptr
	 * 'other.obs_ptr' is cloned in the copy constructor,  and is swapped in the move constructor. These constructors
	 * are necessary to make vectors of EmpDatum (e.g. in EmpData objects) which perform copies/moves of EmpDatum objects
	 * in which case the abstract pointer member obs_ptr has to be handled */

	/** copy constructor */
	EmpDatum(const EmpDatum& other);

	/** move constructor: used by the compiler to move rvalue EmpDatum (other) to this.
	 * Important: Implementing this move constructor causes implicit deletion of the
	 * default copy assignment (operator=). Therefore if this assignment is required
	 * e.g. for allowing assignments for  std::vector<EmpDatum>, we have to implement
	 * this operator explicitly.
	 * note: rvalue reference (&&) is a C++11 extension.
	 **/
	EmpDatum(EmpDatum&& other);

	/** The assignment operator, used for assigning objects. Works for both copy and move assignments. */
	EmpDatum& operator=(EmpDatum other);

	/** destructor: to delete the observation addressed by obs_ptr */
	~EmpDatum();

	/** return a string representation for an EmpDatum in the form [observation | frequency] */
	std::string toString() const;

};



#endif /* SRC_STAT_EMPIRICAL_EMPDATUM_H_ */
