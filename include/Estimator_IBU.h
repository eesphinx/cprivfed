/*
 * Estimator_IBU.h
 *
 *  Created on: Sep 8, 2019
 *      
 */

#ifndef SRC_STAT_ESTIMATORS_ESTIMATOR_IBU_H_
#define SRC_STAT_ESTIMATORS_ESTIMATOR_IBU_H_

/**
 * The class evaluates the maximum likelihood distribution on a given space of secrets (xspace),
 * using empirical data object which describes a collection of locations and the (nonzero) frequency of every location,
 * and the mechanism that produced this empirical data.
 */

#include "StatDistance.h"
#include "DoubleRow.h"
#include "FullGridSpace.h"

class EmpData;
class Estimator_IBU {
public:
	/** default constructor because it is included as a member in e.g. Estimator_AVG_MECH_IBU, Estimator_GIBU */
	Estimator_IBU(){};

	/** Constructs an IBU object from a matrix A (of probabilities) and empirical distribution q */
	Estimator_IBU(const DoubleMatrix& A, const DoubleRow& q);

	/**
	 * Construct an IBU estimator for computing the max-likelihood distribution
	 * on a given space of secrets, with respect to empirical grid data
	 * (with nonzero frequencies) resulting by obfuscation by a given mechanism
	 * The mechanism information is included in the empdata object
	 */
	Estimator_IBU(const FullGridSpace p_xspace, const EmpData& p_empdata);


	/** This implements one iteration of IBU. transforms p_pr_x in the light of the nonzero empirical distribution p_pr_z_giv_x
	 * and the probability matrix p_pr_z_giv_x (restricted to elements in p_emp_pr_z).
	 * Important: This function is invoked from getMaximizer() functions and also from the larger GIBU class */
	DoubleRow update_estimate(const DoubleRow& p_pr_x) const;

	/**
	 * computes the maximum-likelihood estimate on the elements of the secrets space (xspace), using
	 * the available empirical data distribution emp_pr_z, and the reduced G matrix (pr_z_giv_x). This
	 * computation starts with a uniform distribution on xspace. This method is used for debugging as
	 * to see the behavior of the IBU algorithm on a large number of iterations.
	 * @param p_iterations: number of iterations to run
	 * @param p_cycle: a number of iterations after which, the distance and likelihood are recorded.
	 * @param p_distType: type of the distance to be plotted
	 * @param p_xspace_ptr: a pointer to a generic grid space specifying hsize,vsize,cell_side, required if p_distType=EMD
	 * @param p_total_observations : used to calculate the log-likelihood
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
	DoubleRow getMaximizer(const DoubleRow& p_init_pr_x, double p_min_tv, int p_max_iterations) const;

	/**
 	* Starts with p_init_pr_x and iterates until the increase of log_likelihood between two successive
 	* iterations is within p_min_logLikelihood_dff or p_max_iterations is reached.
 	*/
	DoubleRow getMaximizer_by_logLikelihood_per_obs(const DoubleRow& p_init_pr_x, 
	double p_min_logLikelihood_per_obs_dff, int p_max_iterations) const;


	/**
	 * Starts with the uniform distribution on the secrets space and iterates until the
	 * total variation between successive estimates is within p_min_diff or p_max_iterations is reached.
	 */
	DoubleRow getMaximizer(double p_min_tv, int p_max_iterations) const;


	/**
	 * Starts with the uniform distribution on the secrets and iterates for the given number : p_iterations
	 */
	DoubleRow getMaximizer(int p_iterations) const;


	/**
	 * set a reference distribution if needed to track the distance to this
	 * distribution at every iteration. This may be for example the original 
	 * distribution or an MLE computed by another method 
	 */
	void setReferencelDistributionOnSecrets(const DoubleRow& p_ref_distribution);




	/** total number of observations (that produced emp_pr_z) for the empirical data object
	 * used for calculating the likelihood and also for composing IBU iterations in GIBU */
	//int empdata_total_observations;

	/** returns the log_likelihood of a given distribution with respect to the mechanism matrix and per one observation
	 * P.S. The log_likelihood with respect to all observations, the result has to be multiplied by number of observations */
	double get_log_likelihood_per_obs(const DoubleRow& p_pr_x) const;

	/** returns the number of secrets on which a distribution is estimated. Used to initialize distributions */
	int get_number_of_secrets() const;

	/** returns the Fisher information matrix I using the given distribution and number of observations */
	DoubleMatrix getAsympFisher(const DoubleRow& p_pr_x, int observationCount) const;

	/** returns I^-1, that is the inverse of the fisher information matrix of the model
	 * The variance of the MLE = 1/observation * I^-1 */
	DoubleMatrix getAsympFisherInv(const DoubleRow& p_pr_x, int observationCount) const;

	/** returns an approximate variance, also called the expected mean squared error (MSE), also
	 * called the expected norm-2 error of the MLE for the given mechanism and data based on the
	 * asymptotic normality of the MLE. The argument p_pr_x may be the true distribution or an MLE
	 * obtained from IBU  */
	double getAsympRootVariance(const DoubleRow& p_pr_x, int observationCount) const;

	/** returns an approximate variance, also called the expected mean squared error (MSE), also
	 * called the expected norm-2 error of the MLE for the given mechanism and data based on the
	 * asymptotic normality of the MLE. The argument p_pr_x may be the true distribution or an MLE
	 * obtained from IBU  */
	double getApproxVariance1(const DoubleRow& p_pr_x, int observationCount) const;



private:

	/** Reduced G matrix for the empirical data object:rows represent elements of X,
	 * and columns correspond to elements of Z that have nonzero frequencies in the empdata,
	 * where the set Z depends on the empdata object */
	DoubleMatrix pr_z_giv_x;

	/** empirical distribution on Z for the empirical data object.
	 * This vector contains only frequencies of the seen observations, hence it has no zero components */
	DoubleRow emp_pr_z;

	/** A reference distribution on X to be compared to obtained distributions 
	 * (only for debug). This may be for example the original distribution, 
	 * or an MLE computed by another method */
	DoubleRow reference_distribution_on_secrets;

};




#endif /* SRC_STAT_ESTIMATORS_ESTIMATOR_IBU_H_ */
