/*
 * Estimator_INVP.h
 *
 *  Created on: Sep 9, 2019
 *      
 */

#ifndef STAT_ESTIMATORS_ESTIMATOR_AVG_RESULT_INV_H_
#define STAT_ESTIMATORS_ESTIMATOR_AVG_RESULT_INV_H_

/**
 *
 * The estimation is performed by applying INV to every empirical data object with post-processing.
 * Then the resulting distributions are averaged.
 */

#include "Estimator_INV.h"

class Estimator_AVG_RESULT_INV {
public:

	/**
	 * Construct an estimator for the distribution on a given space of secrets,
	 * and a set of empirical data objects resulting from obfuscation using various mechanisms.
	 */
	Estimator_AVG_RESULT_INV (const FullGridSpace p_xspace, const std::vector<EmpData>& p_empdata);


     /*
	 * using normalization 'N' or projection 'P'.
	 */
	DoubleRow getMaximizer (char p_post_processing) const;

private:
	/** a matrix-inversion estimator for each empirical data object */
	std::vector<Estimator_INV> invs;
	
	std::vector<int> empdata_obs_count; // number of all observations in every empdata[i]
	std::vector<std::string> mech_infos; // to save the mechanisms that produced the noisy data

	/** total number of observations (that produced emp_pr_z) for every empirical data object:
	 * used to average the resulting distributions */
//	std::vector<int> empdata_total_observations;

};




#endif /* STAT_ESTIMATORS_ESTIMATOR_AVG_RESULT_INV_H_ */
