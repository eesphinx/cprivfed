#ifndef SRC_STAT_ESTIMATORS_ESTIMATOR_AVG_MECH_IBU_H_
#define SRC_STAT_ESTIMATORS_ESTIMATOR_AVG_MECH_IBU_H_

/** constructs an average mechanism from the individual mechanisms in p_empdata and run IBU on the average mechanism */
#include "Estimator_IBU.h"

class Estimator_AVG_MECH_IBU {
public:

	Estimator_AVG_MECH_IBU (const FullGridSpace p_xspace, const std::vector<EmpData>& p_empdata);

	/**
	 * Starts with p_init_pr_x and iterates until the total variation between two successive
	 * iterations is within p_min_diff or p_max_iterations is reached.
	 */
	DoubleRow getMaximizer(const DoubleRow& p_init_pr_x, double p_min_diff, int p_max_iterations) const;



	/**
	 * Starts with the uniform distribution on the secrets space and iterates until the
	 * total variation between successive estimates is within p_min_diff or p_max_iterations is reached.
	 */
	DoubleRow getMaximizer(double p_min_diff, int p_max_iterations) const;


	/**
	 * Starts with the uniform distribution on the secrets and iterates for the given number : p_iterations
	 */
	DoubleRow getMaximizer(int p_iterations) const;


	/**
	 * set the original distribution if needed to track the distance to the
	 * distribution at every iteration.
	 */
	void setOriginalDistributionOnSecrets(const DoubleRow& p_orig_dist);

private:

	Estimator_IBU ibu;
};




#endif /* SRC_STAT_ESTIMATORS_ESTIMATOR_AVG_MECH_IBU_H_ */
