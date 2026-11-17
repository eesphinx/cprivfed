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

#include <iostream>
#include <cmath>
#include "Distributions.h"
#include "PRRespMechanism.h"
#include "EmpData.h"
#include "Estimator_INV_KRR.h"

/** default constructor (necessary to be used as a member in Estimator_AVG_INV_KRR) */
Estimator_INV_KRR::Estimator_INV_KRR()
:ldp_eps(0){}

/**
 * Construct an estimator on a given fullGrid space from a given 'full' noisy distribution
 * on the secrets (this distribution is computed based on the empirical data). Note that
 * the given distribution may contain 0 components for unseen cells. The estimation relies
 * also on the kRR mechanism.
 */
Estimator_INV_KRR::Estimator_INV_KRR (const FullGridSpace p_xspace, PRRespMechanism p_mech, const DoubleRow& p_full_emp_distribution)
:emp_pr_z(0), ldp_eps(p_mech.get_eps()) {

	// validate that the secret space is compatible with the mechanism and get the mechanism matrix
	if ((p_xspace.hsize != p_mech.hsize)||(p_xspace.vsize != p_mech.vsize)) {
		std::cerr << "The grid of the given full grid space is incompatible with the kRR mechanism grid \n";
		exit(0);
	}

	// validate also that the mechanism is compatible with the given full noisy distribution
	if ((p_full_emp_distribution.getLength() != p_mech.hsize*p_mech.vsize)) {
		std::cerr << "The emprical distribution is incompatible with the kRR mechanism working grid \n";
		exit(0);
	}

	// the empirical distribution: every element of emp_pr_z including zero entries.
	emp_pr_z = p_full_emp_distribution;
}



/**
 * Construct an estimator on a given fullGrid space from noisy data. The mechanism that produced
 * these noisy data must be kRR, otherwise that estimator fails to start.
 */
Estimator_INV_KRR::Estimator_INV_KRR (const FullGridSpace p_xspace, const EmpData& p_empdata)
:emp_pr_z(0), ldp_eps(0) {

	Mechanism* mech_ptr = Mechanism::createMechanism(p_empdata.getObfMechInfo());
	PRRespMechanism* krr_mech_ptr = dynamic_cast<PRRespMechanism*>(mech_ptr);
	if (!krr_mech_ptr) {  // the given mechanism is not local
		std::cerr << "unable to construct INV_KRR estimator: data must be produced by a KRR mechanism";
		delete(mech_ptr); mech_ptr=nullptr;
		exit(1);
	}
	if (krr_mech_ptr->hsize != p_xspace.hsize || krr_mech_ptr->vsize != p_xspace.vsize) {
		std::cerr << "unable to construct INV_KRR estimator: the grid of mechanism ("<< krr_mech_ptr->getTypeStr()
				<<") that produced the empirical data is incompatible with the given space";
		delete(mech_ptr); mech_ptr=nullptr;
		exit(1);
	}
	// Now everything should be fine
	ldp_eps = krr_mech_ptr->get_eps();
	emp_pr_z = p_xspace.get_distribution_from_emp_data(p_empdata);
	delete(mech_ptr); mech_ptr=nullptr;
	//P.S. the full empirical distribution emp_pr_z contains an probability for every element in xspace including those that do not appear in p_emp_data.
}

/**
 * Using matrix inversion method it computes an estimator on the elements of the secrets space (xspace), using
 * the full empirical data distribution emp_pr_z, and the inverse matrix of kRR (See Eq 6 in Kairouz, 2016].
 * If the resulting vector has -ve components, it may returned as it is with post-processing 
 * option 'O', or it is corrected to a valid distribution using normalization 'N' or projection 'P'.
 * @return a distribution on the secrets space, as a row vector
 */
DoubleRow Estimator_INV_KRR::getMaximizer (char p_post_processing) const {
	int k = emp_pr_z.getLength();
	double c1 = (exp(ldp_eps)+k-1)/(exp(ldp_eps)-1);
	double c2 = 1.0/(exp(ldp_eps)-1);
	DoubleRow pr_x = (emp_pr_z.mul(c1)).sub(c2);

	if (p_post_processing == 'O') return pr_x;
	if (p_post_processing == 'N')
	// normalization
		Distributions::normalize(pr_x);
	//System.out.println("euclidean distance to the -ve distribution = "+ pr_x.distance2(pr_x1));
	else if (p_post_processing == 'P')
		Distributions::project_onto_simplex(pr_x);
	/*
	System.out.println("After projection:\n"+pr_x);
	System.out.println("Total probability = "+pr_x. sum());
	System.out.println("Size =" + pr_x.length);
	for (int i=0;i<pr_x.length;i++)
		if (pr_x.get(i)<0)
			System.out.println("-ve element at "+i +":"+pr_x.get(i));
	 */
	//System.out.println("euclidean distance to the -ve distribution = "+ pr_x.distance2(pr_x1));
	else {
		std::cerr << "invalid post-processing parameter ("<< p_post_processing << ") for the KRR-INV estimator";
		exit(0);
	}
	return pr_x;
}
