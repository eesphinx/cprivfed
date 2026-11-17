/*
 * Estimator_GIBU.h
 *
 *  Created on: Sep 8, 2019
 *      
 */

#ifndef SRC_STAT_ESTIMATORS_ESTIMATOR_GIBU_H_
#define SRC_STAT_ESTIMATORS_ESTIMATOR_GIBU_H_

/**
 * This class evaluates the maximum likelihood distribution on a given space of secrets,
 * using multiple empirical data objects. Each object describes a collection of locations and the (nonzero)
 * frequency of every location. Each empirical data object contains the mechanism that produced the empirical data.
 * A GIBU object consists of multiple IBU objects, each one is constructed from one empdata object.
 */

#include "DoubleRow.h"
#include "Estimator_IBU.h"

class EmpData;
class Estimator_GIBU {
public:

	/**
	 * Construct a GIBU estimator for computing the max-likelihood distribution
	 * on a given space of secrets, with respect to empirical grid data objects
	 * (with nonzero frequencies) resulting by obfuscation by a given mechanism
	 * Every Empdata object contains info about the mechanism that produced the data.
	 */
	Estimator_GIBU (const FullGridSpace p_xspace, const std::vector<EmpData>& p_list_empdata);

	/**
	 * Construct a GIBU estimator for given <p_empdata> objects, where
	 * the mechanisms are selected probabilstically according to the conditional 
	 * probabilities matrix <p_omega>, which maps the elements of xspace (rows) to 
	 * the mechanisms in the <p_list_empdata> (columns).
	 */
	Estimator_GIBU (const FullGridSpace p_xspace, const std::vector<EmpData>& p_list_empdata, const DoubleMatrix& p_omega);

	/**
	 * computes the maximum-likelihood estimate on the elements of the secrets space (xspace), using
	 * the available empirical distribution on observations emp_pr_z, and the reduced G matrix
	 * (pr_z_giv_x). This computation starts with a uniform distribution on xspace. This method is used
	 * for debugging as to see the behavior of the EM algorithm on a large number of iterations.
	 * @param p_iterations: number of iterations to run
	 * @param p_cycle: a number of iterations after which, the distance and likelihood are recorded.
	 * @param p_distType: type of the distance to be plotted
	 * @param p_xspace_ptr: a pointer to a generic grid space specifying hsize,vsize,cell_side, required if p_distType=EMD
	 * @param p_logLikelihoodFileName : file to store the progress of the likelihood with the iterations.
	 * @param p_DistanceFileName : to store the progress of the progress of the distance to original distribution.
	 * @param p_DifferenceFileName: to store differences (TV) between successive distributions in iterations (slope)
	 * @return a distribution on the secrets space, as a row vector
	 */
	DoubleRow getMaximizer(int p_iterations, int p_cycle,
			StatDistanceType p_distType, AbstractGridSpace* p_xspace_ptr,
			std::string p_logLikelihoodFileName, std::string p_DistanceFileName, std::string p_DifferenceFileName) const;

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

	/** returns the log_likelihood of a given distribution with respect to the mechanism matrix and per one observation
	 * P.S. The log_likelihood with respect to all observations is obtained by multiplying the result by number of observations */
	double get_log_likelihood_per_obs(const DoubleRow& p_pr_x) const;

private:

	std::vector<Estimator_IBU> ibus;

	/** total number of observations (that produced emp_pr_z) for every empirical data object:
	 * used to average the resulting IBU distributions at each GIBU iterations*/
	std::vector<int> empdata_total_observations;


	/** a row vector representing the original distribution on X to be compared to obtained distributions */
	DoubleRow original_distribution_on_secrets;

	/** one iteration of GIBU: refines p_pr_x */
	DoubleRow update_estimate(const DoubleRow& p_pr_x) const;


};




#endif /* SRC_STAT_ESTIMATORS_ESTIMATOR_GIBU_H_ */
