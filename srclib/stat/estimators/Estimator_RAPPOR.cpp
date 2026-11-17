/*
 * Estimator_RAPPOR.cpp
 *
 *  Created on: Sep 8, 2019
 *      Author: ehab
 */


/**
 * This class computes an estimator for the distribution on a given space of secrets,
 * assuming a Rappor mechanism has be used for obfuscation. This estimator is the implementation of the
 * method in [Kairouz, et al, 2016:discrete distribution estimation under local privacy]. It supports
 * two post-processing methods to correct the negative components: normalization and projection
 **/

#include <iostream>
#include <cmath>
#include "EmpData.h"
#include "AbstractOutputLocalMechanism.h"
#include "BasicOneTimeRappor.h"
#include "Distributions.h"
#include "Estimator_RAPPOR.h"

Estimator_RAPPOR::Estimator_RAPPOR()
:ldp_eps(0) {}

/**
 * Construct an estimator for the distribution on a given space of secrets, assuming a Rappor
 * mechanism has been used for obfuscation. This estimator implements the method in
 * [Kairouz, et al, 2016:discrete distribution estimation under local privacy].
 * @param p_xspace_ptr : pointer to the space of secrets (use ptr because it is abstract class)
 * @param p_empdata: the observed empirical data.
 */
Estimator_RAPPOR::Estimator_RAPPOR(const FullGridSpace p_xspace,
		const EmpData& p_empdata)
:emp_bit_relative_frequencies(0) {

	Mechanism* mech_ptr = Mechanism::createMechanism(p_empdata.getObfMechInfo());
	BasicOneTimeRappor* rap_mech_ptr = dynamic_cast<BasicOneTimeRappor*>(mech_ptr);
	if (!rap_mech_ptr) {  // the given mechanism is not local
		std::cerr << "unable to construct RAPPOR estimator: the mechanism of data object is not RAPPOR";
		delete(mech_ptr); mech_ptr=nullptr;
		exit(1);
	}
	if (rap_mech_ptr->hsize != p_xspace.hsize || rap_mech_ptr->vsize != p_xspace.vsize) {
		std::cerr << "unable to construct a RAPPOR estimator: "
			<< "The grid of mechanism ("<< rap_mech_ptr->getTypeStr()<<")"
			<< "that produced the empirical data object is not compatible with the given space";
		delete(mech_ptr); mech_ptr=nullptr;
		exit(1);
	}
	// Now everything should be fine

	ldp_eps = rap_mech_ptr->ldp_eps;
	delete(mech_ptr); mech_ptr=nullptr;

	// the element j is just Tj/n where Tj is the number of bitVectors that have bit j equal to 1.
	// get for every bit Yj the empirical probability that Yj=1.
	emp_bit_relative_frequencies = p_empdata.get_relative_frequency_of_bits();
}

/** construct a rappor estimator from the emp_bit_relative_frequencies */
Estimator_RAPPOR::Estimator_RAPPOR(double p_ldp_eps, DoubleRow p_emp_bit_relative_frequencies)
:emp_bit_relative_frequencies(p_emp_bit_relative_frequencies), ldp_eps(p_ldp_eps) {
}


/**
 * Returns an estimator on the elements of the secrets space (xspace), using the method in
 * [Kairouz, et al, 2016:discrete distribution estimation under local privacy], which uses
 * the inverse matrix of the k-RR mechanism for every bit (i.e. the kRR is applied on the
 * binary space). Two methods for correcting negative components are supported: normalization
 * (N) and projection (P).
 */
DoubleRow Estimator_RAPPOR::getEstimator (char p_post_processing) const {
	double c1 = exp(0.5*ldp_eps)+1;
	double c2 = exp(0.5*ldp_eps)-1;
	DoubleRow pr_x = emp_bit_relative_frequencies.mul(c1/c2).sub(1/c2);

	if (p_post_processing == 'N') {
		// normalization
		Distributions::normalize(pr_x);
		//System.out.println("euclidean distance to the -ve distribution = "+ pr_x.distance2(pr_x1));
	}

	else if (p_post_processing == 'P') {
		// projection
		Distributions::project_onto_simplex(pr_x);
	}
	else {
		std::cerr << "invalid post-processing parameter ("<< p_post_processing << ") for RAPPOR estimator";
		exit(0);
	}

	return pr_x;
}


/** returns the number of secrets on which a distribution is estimated. Used to initialize distributions */
int Estimator_RAPPOR::get_number_of_secrets() const {
	return emp_bit_relative_frequencies.getLength();
}

