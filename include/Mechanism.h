/*
 * Mechanism.h
 *
 *  Created on: Aug 25, 2025
 *      
 */

#ifndef SRC_MECHANISMS_MECHANISM_H_
#define SRC_MECHANISMS_MECHANISM_H_

/**
 * This is an abstract class describing any mechanism (central or local) that works on a full grid input space.
 * It works on a fullGridSpace. However the output space is arbitrary, e.g. a grid space in the case of 
 * GridOutputLocalMechanism or bitVector (in the case of Rappor). Note that the functions are taking arguments 
 * as abstract Observation objects. Note that all mechanisms in this library must work on a full grid space
 * defined by hsize, vsize. However, the space obfuscated by these mechanisms can be any subset
 * of this full grid space. In other words, we define all mechanisms using FullGridSpace, but
 * we apply them to subclasses of AbstractGridSpace.
 * */

#include <string>

class DoubleColumn;
class Observation;

/** mechanism types */
enum MechanismType {IDENT, VLLAP, PLAP, VLGEOM, PGEOM, 
					TAXICAB, EXP, TC, KRR, PMAT, RMECH, 
					SMECH, BOTRAPPOR, MAXADVRISK,  
					CLAP, CGEOM, CGAUSS};  // These are the central Laplace, Geometric, and Gauss mechanisms

class FullGridSpace;
class Mechanism {
public:
	/** basic mechanism attributes */
	int hsize,vsize;      // The grid for which the mechanism is designed
	//double cell_side;     // The length of the side of one square cell in the grid.
	// The cell side is important only for the mechanisms that are defined by privacy parameters that depends
	// on the distance. In particular geo-indistinuishability mechanisms must be defined with eps and
	// cell_side. However cell_side is not relevant to some mechanisms e.g. kRR, Rappor.


	/* constructor */
	Mechanism(int p_hsize, int p_vsize);

	/* constructs a mechanism with a Full grid input space. However the mechanism output space is generic. */
	Mechanism(const FullGridSpace& p_space);

	/** make a virtual destructor in case we delete a Mechanism pointer (Mechanism *p),
	 * so the proper destructor of the derived object is invoked */
	virtual ~Mechanism(){};

	/** return basic info string from the member mechanism_info */
	virtual std::string get_info() const = 0;
	
	/** return a string consisting of the type and privacy level of the mechanism
	 * The meaning of the privacy level depends on the specific mechanism 
	 * This function must be implemented by all mechanisms
	 * p_precision is number of decimal places for the privacy level. 
	*/
	virtual std::string getTypeAndPriv(const int p_precision) const = 0;


	/** returns the mechanism type as a string*/
	std::string getTypeStr() const;


	/** constructs the mechanism from an info string. It returns a pointer because the return class is abstract */
	static Mechanism* createMechanism(const std::string p_info);


	/** get a string representation for a mechanism type */
	static std::string mechTypeToString(MechanismType p_mech_type);

};


#endif /* SRC_MECHANISMS_MECHANISM_H_ */
