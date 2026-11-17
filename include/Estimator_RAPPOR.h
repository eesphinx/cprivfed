/*
 * Estimator_RAPPOR.h
 *
 *  Created on: Sep 8, 2019
 *      
 */

#ifndef SRC_STAT_ESTIMATORS_ESTIMATOR_RAPPOR_H_
#define SRC_STAT_ESTIMATORS_ESTIMATOR_RAPPOR_H_

/**
 * This class computes an estimator for the distribution on a given space of secrets,
 * assuming a Rappor mechanism has be used for obfuscation. This estimator is the implementation of the
 * method in [Kairouz, et al, 2016:discrete distribution estimation under local privacy]. It supports
 * two post-processing methods to correct the negative components: normalization and projection
 **/

#include "DoubleRow.h"
#include "FullGridSpace.h"

class EmpData;

class Estimator_RAPPOR {
public:

	Estimator_RAPPOR(); // default constructor, used in Estimator_COMBINED_RAPPOR

	/**
	 * Construct an estimator for the distribution on a given space of secrets, assuming a Rappor
	 * mechanism has been used for obfuscation. This estimator implements the method in
	 * [Kairouz, et al, 2016:discrete distribution estimation under local privacy].
	 * @param p_xspace: the space of secrets
	 * @param p_empdata: the observed empirical data.
	 */
	Estimator_RAPPOR(const FullGridSpace p_xspace, const EmpData& p_empdata);


	/** construct a rappor estimator from ldp_eps and emp_bit_relative_frequencies */
	Estimator_RAPPOR(double p_ldp_eps, DoubleRow p_emp_bit_relative_frequencies);

	/**
	 * Returns an estimator on the elements of the secrets space (xspace), using the method in
	 * [Kairouz, et al, 2016:discrete distribution estimation under local privacy].
	 * Two methods for correcting negative components are supported: normalization
	 * (N) and projection (P).
	 */
	DoubleRow getEstimator (char p_post_processing) const;

	/** returns the number of secrets on which a distribution is estimated. Used to initialize distributions */
	int get_number_of_secrets() const;


private:

	//private GridSecretsSpace xspace;   // space of secrets X, which we require to estimate its distribution.

	/** Element j is the number of times that bit j is equal to 1 in the observed bitVectors
	 * divided by the total number (n) of bit vectors in the empirical data [See Kairouz 2016] */
	DoubleRow emp_bit_relative_frequencies;
	double ldp_eps;

};




#endif /* SRC_STAT_ESTIMATORS_ESTIMATOR_RAPPOR_H_ */
