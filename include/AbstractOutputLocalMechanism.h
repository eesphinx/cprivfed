/*
 * LocalMechanism.h
 *
 *  Created on: Sep 24, 2019
 *      
 */

#ifndef SRC_MECHANISMS_ABSTRACTOUTPUTLOCALMECHANISM_H_
#define SRC_MECHANISMS_ABSTRACTOUTPUTLOCALMECHANISM_H_

/**
 * This is an abstract class describing any local mechanism that works on a full grid input space.
 * However the output space is arbitrary, e.g. a grid space in the case of GridOutputLocalMechanism
 * or bitVector (in the case of Rappor). Note that the functions are taking arguments as abstract
 * Observation objects. Note that all mechanisms in this framework must work on a full grid space
 * defined by hsize, vsize. However, the space obfuscated by these mechanisms can be any subset
 * of this full grid space. In other words, we define all mechanisms using FullGridSpace, but
 * we apply them to subclasses of AbstractGridSpace.
 * */

#include <string>
#include "Mechanism.h"

class DoubleColumn;
class Observation;

class FullGridSpace;
class AbstractOutputLocalMechanism : public Mechanism {
public:
	/** basic mechanism attributes */
	//int hsize,vsize;      // The grid for which the mechanism is designed
	//double cell_side;     // The length of the side of one square cell in the grid.
	// The cell side is important only for the mechanisms that are defined by privacy parameters that depends
	// on the distance. In particular geo-indistinuishability mechanisms must be defined with eps and
	// cell_side. However cell_side is not relevant to some mechanisms e.g. kRR, Rappor.


	///* constructor */
	AbstractOutputLocalMechanism(int p_hsize, int p_vsize)
	   : Mechanism(p_hsize, p_vsize) {};   
	   // constructors are not inherited automatically, we have to forward it explicitly to allow 
	   // subclasses to invoke it (they cannot invoke their grandparent constructors directly)

	///* constructs a mechanism with a Full grid input space. However the mechanism output space is generic. */
	//AbstractOutputLocalMechanism(const FullGridSpace& p_space);

	/** make a virtual destructor in case we delete a Mechanism pointer (Mechanism *p),
	 * so the proper destructor of the derived object is invoked */
	//virtual ~AbstractOutputLocalMechanism(){};

	/** returns the probability for input location (specified in linear order)
	 * and an output observation. This must be implemented by every subclass */
	virtual double getProb(int p_in, const Observation* p_out_ptr) const = 0;  //We use a pointer instead of reference to validate casting to subclasses

	/** returns a mechanism column for a given observation: must be implemented by every sub-mechanism */
	virtual DoubleColumn get_mech_column(const Observation* p_out_ptr) const=0;

	/** obfuscate the given location p_loc (in the linear order) into an observation.
	 * We need to use a pointer because Observation is abstract.
	 * Whenever this function is called, the returned pointer should be deleted after use */
	virtual Observation* obfuscateLoc(int p_loc) const = 0;

	/** returns a pointer to newly created observation object constructed from the given string
	 * Note that the string representation of the observation depends on the mechanism, hence
	 * We need to use a pointer because Observation is abstract.
	 * this function must be implemented by every subclass */
	virtual Observation* readObservation(const std::string p_observation_str) const = 0;

	//void add_extra_info(std::string p_extra_info);

	/** return basic info string from the member mechanism_info */
	//virtual std::string get_info() const = 0;
	
	/** return a string consisting of the type and privacy level of the mechanism
	 * The meaning of the privacy level depends on the specific mechanism 
	 * This function must be implemented by all mechanisms
	 * p_precision is number of decimal places for the privacy level. 
	*/
	//virtual std::string getTypeAndPriv(const int p_precision) const = 0;


	/** returns the mechanism type as a string*/
	//std::string get_type_str() const;


	/** constructs the mechanism from an info string. It returns a pointer because the return class is abstract */
	//static AbstractOutputLocalMechanism* createMechanism(const std::string p_info);

		/** get a string representation for a mechanism type */
	//static std::string mechTypeToString(MechanismType p_mech_type);

//private:
//	std::string mechanism_extra_info;

};


#endif /* SRC_MECHANISMS_ABSTRACTOUTPUTLOCALMECHANISM_H_ */
