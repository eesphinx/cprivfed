/*
 * Estimator_INVP.h
 *
 *  Created on: Sep 9, 2019
 *      
 */

#ifndef STAT_ESTIMATORS_ESTIMATOR_INV_H_
#define STAT_ESTIMATORS_ESTIMATOR_INV_H_

/**
 * This estimation method is based on inverting the mechanism matrix. Therefore it assumes that the given
 * mechanism maps the elements of the space to the same space. Since all mechanisms in this software
 * (except Rappor) produce outputs in the full grid, we assume therefore that the input space is also a
 * full grid. The estimated distribution is computed using the inverse of the given mechanism and the
 * 'full' empirical distribution on secrets (i.e. including the zero components that correspond
 * to elements that are not seen).
 *
 * The estimation is performed by multiplying the full empirical distribution (on the space) by the inverse
 * matrix of the mechanism. Then the resulting vector is post-processed using normalization or projection
 * onto the simplex using the method from the paper:
 * "Projection onto the probability simplex: An efficient algorithm with a simple proof, and an application" by Weiran Wang et al.
 */

#include "DoubleRow.h"   // this include the definition of "DoubleMatrix"
#include "FullGridSpace.h"


class GridOutputLocalMechanism; // we just forward declare this class because it is used as a pointer
class GridEmpData;

class Estimator_INV {
public:

	/** default constructor (necessary to be used as a member in Estimator_AVG_INV) */
	Estimator_INV(){};
	/**
	 * Construct an estimator for the distribution on a given space of secrets from a given 'full' noisy
	 * distribution on the secrets (this distribution is assumed to be the empirical distribution based
	 * on observations). Note that the given distribution may contain 0 components which correspond to
	 * unseen cells.
	 * P.S. We use a pointer to the mechanism because it is abstract (contains pure virtual functions).
	 */
	Estimator_INV (const FullGridSpace p_xspace, const GridOutputLocalMechanism* p_local_mech_ptr, const DoubleRow& p_full_emp_distribution);

	/**
	 * Construct an estimator for computing the max-likelihood distribution on a given space of secrets,
	 * an obfuscation mechanism, and empirical data resulting by obfuscation. here we construct from
	 * GridEmpData a (full) empirical distribution which may have zero components (for the cells that are
	 * not observed). P.S. This distribution is obtained by invoking a function in the xspace to encounter
	 * all locations. The obfuscation mechanism is obtained from the empirical data.
	 */
	Estimator_INV (const FullGridSpace p_xspace, const EmpData& p_empdata);

	/** Constructs an INV object from a mechanism matrix A and a full empirical distribution q: used by AVG_MECH_INV estimator */
	Estimator_INV(const DoubleMatrix& A, const DoubleRow& q);

	/**
	 * Returns an MLE on the elements of the Full grid secrets space (xspace) using matrix inversion method.
	 * using the full empirical data distribution emp_pr_z (computed from xspace), and the mechanism matrix
	 * (pr_z_giv_x). The -ve components in the resulting vector are corrected to a valid distribution
	 * using normalization 'N', projection 'P', or no-processinh 'O'.
	 */
	DoubleRow getMaximizer (char p_post_processing='O') const;

	/** returns the number of secrets on which a distribution is estimated. Used to initialize distributions */
	int get_number_of_secrets() const;

private:
	DoubleMatrix pr_z_giv_x;   // Reduced G matrix:rows represent elements of X, columns correspond to elements of Z
	DoubleRow emp_pr_z;     // a row vector, the full empirical distribution on Z (including zero components for unseen locations)

};




#endif /* STAT_ESTIMATORS_ESTIMATOR_INV_H_ */
