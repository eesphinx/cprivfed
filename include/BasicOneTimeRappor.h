/*
 * BasicOneTimeRappor.h
 *
 *  Created on: Jul 3, 2019
 *      
 */

#ifndef PRIV_PLAN_MECHANISMS_BASICONETIMERAPPOR_H_
#define PRIV_PLAN_MECHANISMS_BASICONETIMERAPPOR_H_

#include <vector>
#include "BitVectorObservation.h"
#include "AbstractOutputLocalMechanism.h"


class DoubleColumn;
class FullGridSpace;

/**
 * A OneTimeBasicRappor mechanism on a grid of hsize*vsize. it is a mechanism from k to 2^k where k = hsize*vsize.
 * The input location p (in linear order) is mapped to a k-bit vector, in which all bits are zero except the p-th
 * bit set to 1. Then every bit is obfuscated independently by a 2-RR mechanism to yield a k-bits out_bitVector.
 * The parameter \eps describes the 2-RR mechanism that obfuscates every input bit.
 * The mechanism is called one time because the obfuscation is done once (no longitudinal protection is done), called
 * basic because it does not apply Bloom filters (every input corresponds to a unique one bit in the k-bit Vector B)
 * See [Erlingsson. et al 2014:RAPPOR].
 * The rows (inputs) of the mechanism matrix are ordered in linear order.
 */
class BasicOneTimeRappor : public AbstractOutputLocalMechanism {
public:
    double ldp_eps;

    BasicOneTimeRappor(int p_hsize, int p_vsize, double p_ldp_eps);

    BasicOneTimeRappor(const FullGridSpace& p_space, double p_ldp_eps);

    BasicOneTimeRappor(const std::string p_info);

	/** returns the probability for input location (in the linear order) and an output observation. */
	double getProb(int p_in, const Observation* p_out_ptr) const override;

	/** returns a mechanism column for a given bitVector observation: if p_out_ptr doe not point to an object
	 * of this type, an error is generated */
	DoubleColumn get_mech_column(const Observation* p_out_ptr) const override;

	/** obfuscate the given location p_loc (in the linear order) into an observation (BitVectorObservation)
	 * The returned observation consists of a bit vector having length hsize*vsize */
    BitVectorObservation* obfuscateLoc(int p_loc) const override;

    /** returns a pointer to newly created observation object constructed from the given string
     * Note that the string representation of the observation depends on the mechanism, hence
     * this function must be implemented by every subclass */
    BitVectorObservation* readObservation(const std::string p_observation_str) const override;

	/** return basic info string of the mechanism */
	std::string get_info() const override;
	std::string getTypeAndPriv(const int p_precision) const override;

private:
	/**
     * get the probability of a k-size bitVector from an input location in (in the linear order).
     * The input location X is mapped to a bit vector bitVector(X) which has the X-th bit set
     * to 1 and other bits set to 0. Then the probability of an output bitVector (p_bitVector) is
     * computed according to the bit values in this vector. This function is called only from the
     * other generic version.
     */
    double getProb(int p_in, const std::vector<bool> p_bitVector) const;

	/** returns a mechanism column that corresponds to the observable (vector) p_outBitVector */
    DoubleColumn get_mech_column(const std::vector<bool> p_outBitVector) const;


};




#endif /* PRIV_PLAN_MECHANISMS_BASICONETIMERAPPOR_H_ */
