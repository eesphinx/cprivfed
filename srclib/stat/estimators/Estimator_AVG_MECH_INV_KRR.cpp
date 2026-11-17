/*
 * Estimator_AVG_MECH_INV_KRR.cpp
 *
 *  Created on: Jul 26, 2020
 *      Author: ehab
 */

#include <iostream>
#include "DoubleRow.h"
#include "DoubleColumn.h"
#include "PRRespMechanism.h"
#include "FullGridSpace.h"
#include "EmpData.h"
#include "Distributions.h"
#include "Prior.h"
#include "Estimator_AVG_MECH_INV_KRR.h"




/**
 * This is an instance of AVG_MECH_INV specific to KRR mechanism. It uses direct equations
 * to obtain the estimated distribution using a combined eps computed from the given
 * mechanisms (which produced the noisy empirical data). Similar to INV, AVG_MECH_INV, INV_KRR
 * the post processing is done using normalization or projection onto the simplex using the method from the paper:
 * "Projection onto the probability simplex: An efficient algorithm with a simple proof, and an application" by Weiran Wang et al.
 */

	/**
	 * Construct an estimator for the distribution on a given space of secrets,
	 * and a set of empirical data objects resulting from obfuscation using various mechanisms.
	 * These mechanisms have to be kRR, otherwise the estimator fails to start.
	 * We construct
	 * from the EmpData object the average (full) empirical distribution which may have zero components
	 * (for the cells that are not observed). P.S. This distribution is obtained by invoking a function in the
	 * xspace (for every empdata object) to encounter all locations.
	 * The combined eps is computed from the kRR mechanisms in the empirical data.
	 */
Estimator_AVG_MECH_INV_KRR::Estimator_AVG_MECH_INV_KRR (const FullGridSpace p_xspace, const std::vector<EmpData>& p_list_empdata) {

	// combine empdata objects that have the same mechanism 
	std::vector<EmpData> list_empdata = EmpData::compress(p_list_empdata);

	DoubleRow emp_pr_z(p_xspace.getSize());
	emp_pr_z.fillWith(0);
	double param = 0;

	// obtain the total observations from the given data objects
	int total_observations = 0;
	for (size_t j=0; j< list_empdata.size(); j++)
		total_observations = total_observations + list_empdata[j].get_totalObservations();

	for (size_t j=0; j< list_empdata.size(); j++) {
		if (list_empdata[j].getSize() == 0) {
			std::cerr << "unable to construct AVG_MECH_IBU: data object " << j << " is empty";
			exit(1);
		}
		Mechanism* mech_ptr = Mechanism::createMechanism(list_empdata[j].getObfMechInfo());
		PRRespMechanism* krr_mech_ptr = dynamic_cast<PRRespMechanism*>(mech_ptr);
		if (!krr_mech_ptr) {  // the given mechanism is not a krr
			std::cerr << "unable to construct AVGE_INV_KRR: the mechanism of data object " << j << " is not a kRR";
			delete(mech_ptr); mech_ptr=nullptr;
			exit(1);
		}
		// validate that the secret space is compatible with the mechanism and get the mechanism matrix
		if ((p_xspace.hsize != krr_mech_ptr->hsize)||(p_xspace.vsize != krr_mech_ptr->vsize)) {
			std::cerr << "The grid of the given full grid space is incompatible with the kRR mechanism grid \n";
			delete(mech_ptr); mech_ptr = nullptr;
			exit(1);
		}
		// Now everything should be fine
		// accumulate the variable param which will be used to obtain the combined_eps.
		double ratio = (double) list_empdata[j].get_totalObservations()/total_observations;
		param = param + ratio/(p_xspace.hsize*p_xspace.vsize - 1 + exp( krr_mech_ptr->get_eps() ) );
		emp_pr_z.addi( p_xspace.get_distribution_from_emp_data(list_empdata[j]).mul(ratio)  );
		delete(mech_ptr); mech_ptr=nullptr;
		//P.S. the full empirical distribution emp_pr_z contains an probability for every element in xspace including those that do not appear in p_emp_data.
	}
	double comb_ldp_eps = log(1.0/param - p_xspace.hsize*p_xspace.vsize +1);
	// construct a kRR mechanism using comb_ldp_eps and then use it for estimation
	PRRespMechanism krr_mech(p_xspace.hsize, p_xspace.vsize, comb_ldp_eps);
	inv_krr = Estimator_INV_KRR(p_xspace, krr_mech, emp_pr_z);
}


/**
 * Returns an MLE on the elements of the Full grid secrets space (xspace) using matrix inversion method.
 * using the full empirical data distribution emp_pr_z (computed from xspace), and the mechanism matrix
 * (pr_z_giv_x). The -ve components in the resulting vector are corrected to a valid distribution
 * using normalization 'N' or projection 'P'.
 */
DoubleRow Estimator_AVG_MECH_INV_KRR::getMaximizer (char p_post_processing) const {
	return inv_krr.getMaximizer(p_post_processing);
}





