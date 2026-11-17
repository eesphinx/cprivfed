/*
 * Estimator_KRR_INVN.h
 *
 *  Created on: Sep 10, 2019
 *      
 */

#ifndef SRC_STAT_ESTIMATORS_ESTIMATOR_INV_KRR_H_
#define SRC_STAT_ESTIMATORS_ESTIMATOR_INV_KRR_H_

/**
 * The purpose of this class is to evaluate an estimated distribution for KRR mechanism on a given
 * FullGridSecretSpace xspace. This is an instance of INV running on KRR mechanism.
 * In the estimation procedure the inverse matrix is evaluated using the closed form Equation (6)
 * in [Kairouz, et al 2016] and normalization afterwards.
 * This class uses empirical data or the 'full' distribution from it. Every datum describes an (observed)
 * location in the grid along with its (nonzero) frequency. If a distribution is used it is 'full' in
 * the sense that it is defined on the full grid (i.e. includes zero components for nonobserved locations).
 * The procedure uses also the kRR mechanism that produced the noisy data/distribution.
 */

#include "DoubleRow.h"   // this include the definition of "DoubleMatrix"
#include "FullGridSpace.h"

class PRRespMechanism;
class GridEmpData;


class Estimator_INV_KRR {
public:

	/** default constructor (necessary to be used as a member in Estimator_AVG_INV_KRR) */
	Estimator_INV_KRR();

	/**
	 * Construct an estimator on a given fullGrid space from a given 'full' noisy distribution
	 * on the secrets (this distribution is computed based on the empirical data). Note that
	 * the given distribution may contain 0 components for unseen cells. The estimation relies
	 * also on the kRR mechanism.
	 */
	Estimator_INV_KRR (const FullGridSpace p_xspace, PRRespMechanism p_mech, const DoubleRow& p_full_emp_distribution);



	/**
	 * Construct an estimator on a given fullGrid space from noisy data. The mechanism that produced
	 * these noisy data must be kRR, otherwise that estimator fails to start.
	 */
	Estimator_INV_KRR (const FullGridSpace p_xspace, const EmpData& p_empdata);


	/**
	 * Using matrix inversion method it computes an estimator on the elements of the secrets space (xspace), using
	 * the full empirical data distribution emp_pr_z, and the inverse matrix of kRR (See Eq 6 in Kairouz, 2016].
	 * If the resulting vector has -ve components, it is corrected to a valid distribution using normalization.
	 * @return a distribution on the secrets space, as a row vector
	 */
	DoubleRow getMaximizer (char p_post_processing) const;

private:

    DoubleRow emp_pr_z;            // a row vector, empirical distribution on Z
	double ldp_eps; // privacy parameter


};



#endif /* SRC_STAT_ESTIMATORS_ESTIMATOR_INV_KRR_H_ */
