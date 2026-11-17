/*
 * Estimator_EMP.h
 *
 *  Created on: Dec 28, 2025
 *      
 */

#ifndef STAT_ESTIMATORS_ESTIMATOR_EMP_H_
#define STAT_ESTIMATORS_ESTIMATOR_EMP_H_

#include "DoubleRow.h"

/**
 * 
 * The estimation is performed by aggregating the given hisograms
 * represented as vectors (which may contain negative components). 
 * Then the resulting vector is post-processed using normalization 
 * or projection onto the simplex using the method from the paper:
 * "Projection onto the probability simplex: An efficient algorithm 
 * with a simple proof, and an application" by Weiran Wang et al.
 */

class Estimator_EMP {
public:

	/**
	 * Construct an estimator for the distribution on a given space of secrets,
	 * and a set of histograms resulting from obfuscation using central mechanisms. 
	 */
	Estimator_EMP (const std::vector<DoubleRow>& p_histograms);

	/**
	 * Returns an estimate by aggregating the underlying histograms, 
     * and post-processing using 
     * 'N': truncating negative components to 0 and then normalize, 
     * 'P': normalize then project onto the simplex, Or 
     * 'O': for only normalization without correcting -ve components.
	 */
	DoubleRow getEstimate (char p_post_processing) const;

	/** apply weights using the given variances before post-processing */
	DoubleRow getEstimate (const std::vector<double> p_list_variances, char p_post_processing) const;
    
	/** returns the pseudo estimate defined by the aggregae divided by p_count */
	DoubleRow getPseudoEstimate(int p_count) const;

	/** 
	 * returns the pseudo estimate defined by the weighted aggregation (obtained from the given variances) 
	 * divided by the weighted count computed from <subcounts> 
	 * The vectors <p_list_variances>, <p_subcounts> and internal <histograms> must have the same size. 
	*/
	DoubleRow getPseudoEstimate(const std::vector<double> p_list_variances, std::vector<int> p_subcounts) const;

	/** analytic Mean SQD without weighting */
    static double getPseudoEstimateAnalyticMeanSQD(const DoubleRow p_original_distribution, const std::vector<double> p_list_variances, int p_numPoints);


	/** analytic Mean SQD with weighting */
    static double getPseudoEstimateAnalyticMeanSQD(const DoubleRow p_original_distribution, const std::vector<double> p_list_variances, std::vector<int> p_subcounts);

private:
    std::vector<DoubleRow> histograms;

// 	/**
// 	 * Returns the aggregation of the given vectors, 
//      * and post-processing using 
//      * 'N': truncating negative components to 0 and then normalize, 
//      * 'P': normalize then project onto the simplex, Or 
//      * 'O': for only normalization without correcting -ve components.
// 	 */
// 	DoubleRow getAggregation (const std::vector<DoubleRow>& vecs, char p_post_processing) const;
};




#endif /* STAT_ESTIMATORS_ESTIMATOR_EMP_H_ */
