
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

#include "iostream"
#include "DoubleRow.h"
#include "GridOutputLocalMechanism.h"
#include "EmpData.h"
#include "Distributions.h"
#include "Estimator_INV.h"


/**
 * Construct an estimator for the distribution on a given space of secrets from a given 'full' noisy
 * distribution on the secrets (this distribution is assumed to be the empirical distribution based
 * on observations). Note that the given distribution may contain 0 components which correspond to
 * unseen cells.
 * P.S. We use a pointer to the mechanism because it is abstract (contains pure virtual functions).
 */
Estimator_INV::Estimator_INV (const FullGridSpace p_xspace, const GridOutputLocalMechanism* p_local_mech_ptr, const DoubleRow& p_full_emp_distribution) {
	
	// validate that the secret space is compatible with the mechanism and get the mechanism matrix
	if (p_local_mech_ptr->hsize != p_xspace.hsize || p_local_mech_ptr->vsize != p_xspace.vsize) {
		std::cerr << "unable to construct INV estimator: the grid of given mechanism ("<< 
		          p_local_mech_ptr->getTypeStr() <<
				  ") that produced the empirical data is incompatible with the given space";
		exit(1);
	}

	// validate also that the mechanism is compatible with the given full noisy distribution
	if ((p_full_emp_distribution.getLength() != p_local_mech_ptr->hsize*p_local_mech_ptr->vsize)) {
		std::cerr << "The emprical distribution is incompatible with the mechanism working grid \n";
		exit(0);
	}
	pr_z_giv_x = p_local_mech_ptr -> get_mech_matrix();
	emp_pr_z = p_full_emp_distribution;
}


/**
 * Construct an estimator for computing the max-likelihood distribution on a given space of secrets,
 * an obfuscation mechanism, and empirical data resulting by obfuscation. here we construct from empdata
 * a (full) empirical distribution which may have zero components (for the cells that are not observed).
 * P.S. This distribution is obtained by invoking a function in the xspace to encounter all locations.
 * The obfuscation mechanism is obtained from the empirical data.
 */
Estimator_INV::Estimator_INV (const FullGridSpace p_xspace, const EmpData& p_empdata) {

	Mechanism* mech_ptr = Mechanism::createMechanism(p_empdata.getObfMechInfo());
	GridOutputLocalMechanism* local_mech_ptr = dynamic_cast<GridOutputLocalMechanism*>(mech_ptr);
	if (!local_mech_ptr) {  // the given mechanism is not local
		std::cerr << "unable to construct INV estimator: the mechanism of data object is not local grid output mechansim";
		delete(mech_ptr); mech_ptr=nullptr;
		exit(1);
	}
	if (local_mech_ptr->hsize != p_xspace.hsize || local_mech_ptr->vsize != p_xspace.vsize) {
		std::cerr << "unable to construct INV estimator: the grid of mechanism ("<< local_mech_ptr->getTypeStr()
				<<") that produced the empirical data is incompatible with the given space";
		delete(mech_ptr); mech_ptr=nullptr;
		exit(1);
	}
	// Now everything should be fine
	pr_z_giv_x = local_mech_ptr-> get_mech_matrix();
	emp_pr_z = p_xspace.get_distribution_from_emp_data(p_empdata);
	delete(mech_ptr); mech_ptr=nullptr;
	//P.S. the full empirical distribution: contains an probability for every element in xspace including those that do not appear in p_emp_data.
}

/** Constructs an INV object from a mechanism matrix A and a full empirical distribution q */
Estimator_INV::Estimator_INV(const DoubleMatrix& A, const DoubleRow& q) {
	if (A.getColumns() != q.getLength()) {
		std::cerr << "error: given A and q must have the same columns";
		exit(0);
	}
	pr_z_giv_x = A; emp_pr_z = q;
}

/**
 * Returns an MLE on the elements of the Full grid secrets space (xspace) using matrix inversion method.
 * using the full empirical data distribution emp_pr_z (computed from xspace), and the mechanism matrix
 * (pr_z_giv_x). The -ve components in the resulting vector are corrected to a valid distribution
 * using normalization 'N', projection 'P', or no-processing 'O'.
 */
DoubleRow Estimator_INV::getMaximizer (char p_post_processing) const {
	//std::cout << "empirical = " << emp_pr_z << "\n";
	DoubleMatrix inv = pr_z_giv_x.getInverse();
	DoubleRow pr_x = emp_pr_z.mmul(inv);
	//pr_x.write_to_file("emp_dist_times_inv.txt");    // just to inspect the product before post-processing
	//std::cout << "estimated = " << pr_x.toString() << "\n";
	//std::cout << "prob sum = " << pr_x.sum() << "\n";

	if (p_post_processing == 'O') return pr_x;
	if (p_post_processing == 'N')
		Distributions::normalize(pr_x);
	else if (p_post_processing == 'P')
		Distributions::project_onto_simplex(pr_x);
	else {
		std::cerr << "invalid post-processing parameter ("<< p_post_processing << ") for the INV estimator";
		exit(0);
	}

	return pr_x;
}

/** returns the number of secrets on which a distribution is estimated. Used to initialize distributions */
int Estimator_INV::get_number_of_secrets() const {
	return pr_z_giv_x.getRows();
}

