/**
 * A OneTimeBasicRappor mechanism on a grid of hsize*vsize. It is a mechanism from k to 2^k where k = hsize*vsize.
 * The input location p (in linear order) is mapped to a k-bit vector, in which all bits are zero except the p-th 
 * bit set to 1. Then every bit is obfuscated independently by a 2-RR mechanism to yield a k-bits out_bitVector.
 * The parameter \eps describes the 2-RR mechanism that obfuscates every input bit. 
 * The mechanism is called one time because the obfuscation is done once (no longitudinal protection is done), called 
 * basic because it does not apply Bloom filters (every input corresponds to a unique one bit in the k-bit Vector B)
 * See [Erlingsson. et al 2014:RAPPOR]. 
 * The rows (inputs) of the mechanism matrix are ordered in linear order. 
 */

#include <iostream>
#include <sstream>
#include <cmath>
#include "DoubleRow.h"
#include "DoubleColumn.h"
#include "Distributions.h"
#include "FullGridSpace.h"
#include "BitVectorObservation.h"
#include "BasicOneTimeRappor.h"
#include "GSettings.h"

BasicOneTimeRappor::BasicOneTimeRappor(int p_hsize, int p_vsize, double p_ldp_eps)
:AbstractOutputLocalMechanism(p_hsize, p_vsize), ldp_eps(p_ldp_eps){
}

BasicOneTimeRappor::BasicOneTimeRappor(const FullGridSpace& p_space, double p_ldp_eps)
:AbstractOutputLocalMechanism(p_space.hsize, p_space.vsize), ldp_eps(p_ldp_eps){   // we set cell_side =1;
}

BasicOneTimeRappor::BasicOneTimeRappor(const std::string p_info)
:AbstractOutputLocalMechanism(0,0) {   // superclass must be constructed: hsize,vsize will be modified later.
	std::string mech_type_str = AbstractOutputLocalMechanism::mechTypeToString(BOTRAPPOR);
	std::istringstream iss(p_info);
	std::string p_mech_type_str; double p_ldp_eps;
	if (iss >> p_mech_type_str >> hsize >> vsize >> p_ldp_eps) {
		if (p_mech_type_str != mech_type_str) {
			std::cerr << "unable to construct " << mech_type_str << " mechanism from info string:\n "
					<< "info string must start with " << mech_type_str;
			exit(0);
		}
		// construction work
		ldp_eps = p_ldp_eps;
	}
	else {
		std::cerr << "unable to parse the mechanism info to construct "<< mech_type_str <<" mechanism";
		exit(0);
	}
}

/** returns the probability for input location (specified in the linear order) and an output
 * bitVector observation. If the given observation is not a bitVector observation, an error
 * is produced */
double BasicOneTimeRappor::getProb(int p_in, const Observation* p_out_ptr) const {
	if(const BitVectorObservation* bitVec_out_ptr = dynamic_cast<const BitVectorObservation*>(p_out_ptr))
		return getProb(p_in, bitVec_out_ptr->bitVector); // safe to call
	else {
		std::cerr << "unable to cast a generic observation to bit vector observation to get the probability";
		exit(0);
	}
}

/** returns a mechanism column for a given bitVector observation: if p_out_ptr doe not point to an object
 * of this type, an error is generated */
DoubleColumn BasicOneTimeRappor::get_mech_column(const Observation* p_out_ptr) const {
	if(const BitVectorObservation* grid_out_ptr = dynamic_cast<const BitVectorObservation*>(p_out_ptr))
    {
        std::vector<bool> out = grid_out_ptr->bitVector; // safe to call
        return get_mech_column(out);
    }
	else {
		std::cerr << "unable to cast a generic observation to bit vector observation to get a mechanism column";
		exit(0);
	}
}



/** obfuscate the given location p_loc (in the linear order) into an observation (BitVectorObservation)
 * The returned observation consists of a bit vector having length hsize*vsize */
BitVectorObservation* BasicOneTimeRappor::obfuscateLoc(int p_loc) const {
	int vecLen = hsize*vsize;
	std::vector<bool> bitVector(vecLen);
	DoubleRow distribution(2);
	// Let Xj be bit j in the vector representing original input
	// and Yj be bit j in the output vector after obfuscation
	distribution.put(0, exp(0.5 * ldp_eps)/(1+exp(0.5*ldp_eps)));  // probability of Yj=Xj
	distribution.put(1, 1/(1+ exp(0.5 * ldp_eps))); // probability of Yj!=Xj

	std::vector<int> bitRetention = Distributions::sampleDiscrete(distribution, vecLen); // 0 indicates no flip, 1 indicates flip
	for (int j=0; j<vecLen;j++)
		bitVector[j]= (bitRetention[j]==0) ? false : true; // handling the bit at j which is originally false.
	bitVector[p_loc]= (bitRetention[p_loc]==0) ? true : false; // handling the bit at p_loc which is originally true.
	return new BitVectorObservation(bitVector);  // create a new observation and return a pointer to it.
}


/** returns a pointer to newly created observation object constructed from the given string
 * Note that the string representation of the observation depends on the mechanism, hence
 * this function must be implemented by every subclass */
BitVectorObservation* BasicOneTimeRappor::readObservation(std::string p_observation_str) const {
	std::istringstream iss(p_observation_str);
	std::vector<bool>  bitVect;
	std::string token; // expected "0" or "1"
	for (int i=0; i< hsize*vsize; i++)  // read a vector of length hsize*vsize
		if (iss >> token) {
			if (token =="1") bitVect.push_back(true);
			else if (token == "0") bitVect.push_back(false);
			else {
				std::cerr << "unable to read bit vector observation from a string: invalid format";
				exit(0);
			}
		}
	return new BitVectorObservation(bitVect);
}

/**
 * get the probability of a k-size bitVector from an input location in_x, in_y. The input location X is mapped to a bit vector bitVector(X)
 * which has the X-th bit set to 1 and other bits set to 0. Then the probability of an output bitVector (p_bitVector) is
 * computed according to the bit values in this vector.
 */
double BasicOneTimeRappor::getProb(int p_in, const std::vector<bool> p_out_bitVector) const {
	double c1 = exp(0.5*ldp_eps)/(1 + exp(0.5*ldp_eps));  // probability that the bit remains the same
	double c2 = 1/(1 + exp(0.5*ldp_eps));  // probability that the bit is changed

	double prob;
	// compute the probability due the in-th bit. Note that p_in-th bit in the original bitVector is 1.
	if (p_out_bitVector[p_in])
		prob = c1;
	else prob = c2;
	// compute the probabilities due other bits. Note that these bits in the original bitVector are 0.
	for (int i= 0; i<p_in;i++)
		if (p_out_bitVector[i])
			prob = prob * c2;
		else prob = prob * c1;
	for (int i= p_in+1; i<hsize*vsize;i++)
		if (p_out_bitVector[i])
			prob = prob * c2;
		else prob = prob * c1;

	return prob;
}

/** returns a mechanism column that corresponds to the observable (vector) p_outBitVector */
DoubleColumn BasicOneTimeRappor::get_mech_column(const std::vector<bool> p_outBitVector) const {
	DoubleColumn retCol(hsize*vsize);
	for (int i=0;i<hsize*vsize; i++)
		retCol.put(i, getProb(i,p_outBitVector));
	return retCol;
}

/** return basic info string of the mechanism */
std::string BasicOneTimeRappor::get_info() const {
	return AbstractOutputLocalMechanism::mechTypeToString(BOTRAPPOR) + " " + std::to_string(hsize) + " " +
			std::to_string(vsize) + " "+ std::to_string(ldp_eps);
}

std::string BasicOneTimeRappor::getTypeAndPriv(const int p_precision) const {
	return AbstractOutputLocalMechanism::mechTypeToString(BOTRAPPOR) + 
	GSettings::to_string_with_precision(ldp_eps, p_precision);
}



