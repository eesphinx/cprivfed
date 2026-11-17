/*
 * Estimator_COMBINED_RAPPOR.cpp
 *
 *  Created on: Jul 26, 2020
 *      Author: ehab
 */

#include <iostream>
#include <cmath>
#include "DoubleRow.h"
#include "DoubleColumn.h"
#include "BasicOneTimeRappor.h"
#include "EmpData.h"
#include "Distributions.h"
#include "Prior.h"
#include "Estimator_COMBINED_RAPPOR.h"




/**
 * This estimator is based on various Rappor mechanisms with different privacy parameters.
 * The distribution is estimated using a combined eps computed from the given mechanisms.
 * Similar to INV, INV_KRR, AVG_MECH_INV, AVG_MECH_INV_KRR
 * the post processing is done using normalization or projection onto the simplex using the method from the paper:
 * "Projection onto the probability simplex: An efficient algorithm with a simple proof, and an application" by Weiran Wang et al.
 */

	/**
	 * Construct an estimator for the distribution on a given space of secrets,
	 * and a set of empirical data objects resulting from obfuscation using various mechanisms.
	 * These mechanisms have to be RAPPOR, otherwise the estimator fails to start.
	 * The combined eps is computed from the RAPPOR mechanisms in the empirical data.
	 */
Estimator_COMBINED_RAPPOR::Estimator_COMBINED_RAPPOR (const FullGridSpace p_xspace, const std::vector<EmpData>& p_list_empdata) {

	// combine empdata objects that have the same mechanism 
	std::vector<EmpData> list_empdata = EmpData::compress(p_list_empdata);

	DoubleRow emp_bit_relative_frequencies(p_xspace.getSize());
	emp_bit_relative_frequencies.fillWith(0);
	double param = 0; // used to calculate the combined eps

	// obtain the total observations from the given data objects
	int total_observations = 0;
	for (size_t j=0; j< list_empdata.size(); j++)
		total_observations = total_observations + list_empdata[j].get_totalObservations();

	for (size_t j=0; j< list_empdata.size(); j++) {
		Mechanism* mech_ptr = Mechanism::createMechanism(list_empdata[j].getObfMechInfo());
		BasicOneTimeRappor* rap_mech_ptr = dynamic_cast<BasicOneTimeRappor*>(mech_ptr);
		if (!rap_mech_ptr) {  // the given mechanism is not local
			std::cerr << "unable to construct COMBINED_RAPPOR estimator: the mechanism of data object " 
			<< j << " is not RAPPOR";
			delete(mech_ptr); mech_ptr=nullptr;
			exit(1);
		}
		if (rap_mech_ptr->hsize != p_xspace.hsize || rap_mech_ptr->vsize != p_xspace.vsize) {
			std::cerr << "unable to construct an COMBINED_RAPPOR estimator: "
				<< "The grid of mechanism ("<< rap_mech_ptr->getTypeStr()<<")"
				<< "that produced the empirical data object " << j << "is not compatible with the given space";
			delete(mech_ptr); mech_ptr=nullptr;
			exit(1);
		}
		// Now everything should be fine

		//std::cout << "rappor " << j << " with eps = " << rappor_mech_ptr ->ldp_eps << "\n";
		double ratio = (double) list_empdata[j].get_totalObservations()/total_observations;
		// accumulate the variable param which will be used to obtain the combined_eps.
		param = param + ratio/(1 + exp( 0.5*rap_mech_ptr->ldp_eps  ) );
		emp_bit_relative_frequencies.addi(list_empdata[j].get_relative_frequency_of_bits().mul(ratio));
		delete(mech_ptr); mech_ptr=nullptr;
	}
	double comb_ldp_eps = 2.0 * log(1.0/param - 1);
	estimator_rappor = Estimator_RAPPOR(comb_ldp_eps, emp_bit_relative_frequencies);
	//std::cout << "combined rappor with eps = " << comb_ldp_eps << "\n";
}


/**
 * Returns an MLE on the elements of the Full grid secrets space (xspace) using matrix inversion method.
 * using the full empirical data distribution emp_pr_z (computed from xspace), and the mechanism matrix
 * (pr_z_giv_x). The -ve components in the resulting vector are corrected to a valid distribution
 * using normalization 'N' or projection 'P'.
 */
DoubleRow Estimator_COMBINED_RAPPOR::getMaximizer (char p_post_processing) const {
	return estimator_rappor.getEstimator(p_post_processing);
}




