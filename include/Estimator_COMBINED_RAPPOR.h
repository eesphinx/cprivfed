/*
 * Estimator_COMBINED_RAPPOR.h
 *
 *  Created on: Jul 26, 2020
 *      
 */

#ifndef SRC_STAT_ESTIMATORS_ESTIMATOR_COMBINED_RAPPOR_H_
#define SRC_STAT_ESTIMATORS_ESTIMATOR_COMBINED_RAPPOR_H_

/**
 *
 * This estimator is based on various Rappor mechanisms with different privacy parameters.
 * The distribution is estimated using a combined eps computed from the given mechanisms.
 * Similar to INV, INV_KRR, AVG_MECH_INV, AVG_MECH_INV_KRR
 * the post processing is done using normalization or projection onto the simplex using the method from the paper:
 * "Projection onto the probability simplex: An efficient algorithm with a simple proof, and an application" by Weiran Wang et al.
 */

#include "Estimator_RAPPOR.h"

class Estimator_COMBINED_RAPPOR {
public:

	/**
	 * Construct an estimator for the distribution on a given space of secrets,
	 * and a set of empirical data objects resulting from obfuscation using various mechanisms.
	 * These mechanisms have to be RAPPOR, otherwise the estimator fails to start.
	 * The combined eps is computed from the RAPPOR mechanisms in the empirical data.
	 */
	Estimator_COMBINED_RAPPOR (const FullGridSpace p_xspace, const std::vector<EmpData>& p_empdata);


	/**
	 * Returns the estimate using an Estimator_RAPPOR object running on the combined eps.
	 * The -ve components in the resulting vector are corrected to a valid distribution
	 * using normalization 'N' or projection 'P'.
	 */
	DoubleRow getMaximizer (char p_post_processing) const;

private:

	Estimator_RAPPOR estimator_rappor;
};




#endif /* SRC_STAT_ESTIMATORS_ESTIMATOR_COMBINED_RAPPOR_H_ */
