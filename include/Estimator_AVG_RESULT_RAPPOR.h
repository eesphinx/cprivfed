/*
 * Estimator_AVG_RESULT_RAPPOR.h
 *
 *  Created on: Aug 3, 2020
 *      
 */

#ifndef SRC_STAT_ESTIMATORS_ESTIMATOR_AVG_RESULT_RAPPOR_H_
#define SRC_STAT_ESTIMATORS_ESTIMATOR_AVG_RESULT_RAPPOR_H_


/**
 *
 * The estimation is performed by applying RAPPOR estimator to every empirical data object with post-processing.
 * Then the resulting distributions are averaged.
 */

#include "Estimator_RAPPOR.h"

class Estimator_AVG_RESULT_RAPPOR {
public:

	/**
	 * Construct an estimator for the distribution on a given space of secrets,
	 * and a set of empirical data objects resulting from obfuscation using various mechanisms.
	 */
	Estimator_AVG_RESULT_RAPPOR (const FullGridSpace p_xspace, const std::vector<EmpData>& p_empdata);


     /*
	 * using normalization 'N' or projection 'P'.
	 */
	DoubleRow getMaximizer (char p_post_processing) const;

private:
	/** a RAPPOR estimator for each empirical data object */
	std::vector<Estimator_RAPPOR> raps;
	std::vector<double> weights;  // number of observations in empdata[i] divided by the total number of observations

};





#endif /* SRC_STAT_ESTIMATORS_ESTIMATOR_AVG_RESULT_RAPPOR_H_ */
