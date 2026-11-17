/*
 * Estimator_AVG_MECH_INV_KRR.h
 *
 *  Created on: Jul 26, 2020
 *      
 */

#ifndef SRC_STAT_ESTIMATORS_ESTIMATOR_AVG_MECH_INV_KRR_H_
#define SRC_STAT_ESTIMATORS_ESTIMATOR_AVG_MECH_INV_KRR_H_

/**
 * This is an instance of AVG_MECH_INV specific to KRR mechanism. It uses direct equations
 * to obtain the estimated distribution using a combined eps computed from the given
 * mechanisms (which produced the noisy empirical data). Similar to INV, AVG_MECH_INV, INV_KRR
 * the post processing is done using normalization or projection onto the simplex using the method from the paper:
 * "Projection onto the probability simplex: An efficient algorithm with a simple proof, and an application" by Weiran Wang et al.
 */

#include "Estimator_INV_KRR.h"

class Estimator_AVG_MECH_INV_KRR {
public:
	/** default constructor (necessary to be used as a member in Estimator_AVG_INV) */
	Estimator_AVG_MECH_INV_KRR(){};

	/**
	 * Construct an estimator for the distribution on a given space of secrets,
	 * and a set of empirical data objects resulting from obfuscation using various mechanisms.
	 * These mechanisms have to be kRR, otherwise the estimator fails to start.
	 * It constructs from the EmpData object the average (full) empirical distribution which
	 * may have zero components (for the cells that are not observed). P.S. This distribution
	 * is obtained by invoking a function in the xspace (for every empdata object) to encounter
	 * all locations. The combined eps is computed from the kRR mechanisms in the empirical data.
	 */
	Estimator_AVG_MECH_INV_KRR(const FullGridSpace p_xspace, const std::vector<EmpData>& p_empdata);


	/**
	 * Returns the estimate using a Estimator_INV_KRR running on the combined eps and the
	 * full empirical data distribution emp_pr_z. The -ve components in the resulting vector
	 * are corrected to a valid distribution using normalization 'N' or projection 'P'.
	 */
	DoubleRow getMaximizer (char p_post_processing) const;

private:

	Estimator_INV_KRR inv_krr;
};




#endif /* SRC_STAT_ESTIMATORS_ESTIMATOR_AVG_MECH_INV_KRR_H_ */
