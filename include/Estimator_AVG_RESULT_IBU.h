/*
 * Estimator_EM.h
 *
 *  Created on: Sep 8, 2019
 *      
 */

#ifndef SRC_STAT_ESTIMATORS_ESTIMATOR_AVG_RESULT_IBU_H_
#define SRC_STAT_ESTIMATORS_ESTIMATOR_AVG_RESULT_IBU_H_

/**
 * This class runs IBU on every given empirical data and returns the weighted average of resulting distributions
 */

#include <vector>
#include "Estimator_IBU.h"


class Estimator_AVG_RESULT_IBU {
public:

	/**
	 * IBU is applied to each emp data object and the resulting distributions are averaged.
	 */
	Estimator_AVG_RESULT_IBU (const FullGridSpace p_xspace, const std::vector<EmpData>& p_empdata);

	/**
	 * Starts with the uniform distribution on the secrets space and iterates until the
	 * total variation between successive estimates is within p_min_diff or p_max_iterations is reached.
	 */
	DoubleRow getMaximizer(double p_min_diff, int p_max_iterations) const;

private:

	std::vector<Estimator_IBU> ibus;

	// total number of observations for every empirical data object:
	std::vector<int> empdata_total_observations;
};




#endif /* SRC_STAT_ESTIMATORS_ESTIMATOR_AVG_RESULT_IBU_H_ */
