/*
 * PRRespMechanism.h
 *
 *  Created on: Sep 7, 2019
 *      
 */

#ifndef SRC_PRIV_PLAN_MECHANISMS_PRRESPMECHANISM_H_
#define SRC_PRIV_PLAN_MECHANISMS_PRRESPMECHANISM_H_


#include "GridOutputLocalMechanism.h"

class PrivMetric;  // forward declaration because we have a pointer to this class (a member and parameter)
class Loss;        // forward declaration because we use a reference object in a parameter.

class PRRespMechanism: public GridOutputLocalMechanism {
public:
	/**
	 * This mechanism is parameterized by local differential privacy eps	 */
	PRRespMechanism(int p_hsize, int p_vsize, double p_ldp_eps);

	/** The following constructs a mechanism that satisfies a given loss value for the given prior
	 *  and given loss object. This method is used to compare mechanisms with respect
	 *  to statistical utility (using em algorithm). Important: we assume that loss(i,i)=0. */
	PRRespMechanism(const FullGridSpace& p_space, const Prior& p_prior, Loss* p_loss_ptr, double p_exp_loss);

	PRRespMechanism(const std::string p_info);

	/** obfuscate the given location p_loc (in the linear order) into an observation
	 * This overrides GridOutputMechanism to improve the performance, taking advantage
	 * of the uniform distribution over hsize*vsize -1 elements */
	GridObservation* obfuscateLoc(int p_loc) const override;

	double get_eps() const;

	std::string get_info() const override;
	std::string getTypeAndPriv(const int p_precision) const override;

private:

	double prob_diag;   // probability of yielding i given the secret is i
	double prob_non_diag;   // probability of yielding j given the secret is i, j noteq i

	double getProb(int p_in, int p_out) const override;
};


#endif /* SRC_PRIV_PLAN_MECHANISMS_PRRESPMECHANISM_H_ */
