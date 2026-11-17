/*
 * Estimator_AVG_MECH_INV.h
 *
 *  Created on: Sep 9, 2019
 *      
 */

#ifndef STAT_ESTIMATORS_ESTIMATOR_AVG_MECH_INV_H_
#define STAT_ESTIMATORS_ESTIMATOR_AVG_MECH_INV_H_

/**
 * This estimation method is based on inverting the mechanism matrix. Therefore it assumes that the given
 * mechanisms have the same alphabet of observables. Since all mechanisms in this software
 * (except Rappor) produce outputs in the full grid, we assume therefore that the input space is also a
 * full grid. The estimated distribution is computed using the inverse of the average mechanism and the
 * 'full' empirical distribution on secrets (i.e. including the zero components that correspond
 * to elements that are not seen).
 *
 * The estimation is performed by multiplying the full empirical distribution (on the space) by the inverse
 * matrix of the average mechanism. Then the resulting vector is post-processed using normalization or projection
 * onto the simplex using the method from the paper:
 * "Projection onto the probability simplex: An efficient algorithm with a simple proof, and an application" by Weiran Wang et al.
 */

#include "Estimator_INV.h"
#include "Estimator_AVG_MECH_INV_KRR.h"

class Estimator_AVG_MECH_INV {
public:

	/**
	 * Construct an estimator for the distribution on a given space of secrets,
	 * and a set of empirical data objects resulting from obfuscation using various mechanisms. We construct
	 * from the EmpData object the average (full) empirical distribution which may have zero components
	 * (for the cells that are not observed). P.S. This distribution is obtained by invoking a function in the
	 * xspace (for every empdata object) to encounter all locations.
	 * An average mechanism is obtained from the individual mechanisms of the given empirical data.
	 */
	Estimator_AVG_MECH_INV (const FullGridSpace p_xspace, const std::vector<EmpData>& p_empdata);


	/**
	 * Returns an MLE on the elements of the Full grid secrets space (xspace) using matrix inversion method.
	 * using the full empirical data distribution emp_pr_z (computed from xspace), and the mechanism matrix
	 * (pr_z_giv_x). The -ve components in the resulting vector are corrected to a valid distribution
	 * using normalization 'N' or projection 'P', Or 'O' for no-post processing. 
	 */
	DoubleRow getMaximizer (char p_post_processing) const;

private:

	Estimator_INV inv;
	Estimator_AVG_MECH_INV_KRR avg_mech_inv_krr;  // an estimator optimized for krr
	bool all_mechanisms_are_krr; //set to true if all mechanisms are kRR
};




#endif /* STAT_ESTIMATORS_ESTIMATOR_AVG_MECH_INV_H_ */
