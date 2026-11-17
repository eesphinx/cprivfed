
#include <iostream>
#include <sstream>
#include <string>
#include <cmath>
#include "Prior.h"
#include "Loss.h"
#include "Distributions.h"
#include "DoubleRow.h"
#include "FullGridSpace.h"
#include "PRRespMechanism.h"
#include "GSettings.h"

/*
 * Constructs a KRR mechanism that satisfies eps-local differential privacy.
 */
PRRespMechanism::PRRespMechanism(int p_hsize, int p_vsize, double p_ldp_eps)
:GridOutputLocalMechanism(p_hsize, p_vsize) {
	double t= exp(p_ldp_eps);
	prob_diag = t/(hsize*vsize-1+t);
	prob_non_diag = 1/(hsize*vsize-1+t);
}

/** The following constructs a mechanism that satisfies a given loss value for the given prior
 *  and given loss object. This method is used to compare mechanisms with respect
 *  to statistical utility (using em algorithm). Important: we assume that loss(i,i)=0. */
PRRespMechanism::PRRespMechanism(const FullGridSpace& p_space, const Prior& p_prior, Loss* p_loss_ptr, double p_exp_loss)
:GridOutputLocalMechanism(p_space.hsize, p_space.vsize) {
	int k = hsize*vsize;
	// prob_diag = 1 - (k-1) * prob_non_diag.
	// exp_loss = sum_{i,j} pi_i A_ij L_ij = c1 * prob_diag + c2 * prob_non_diag
	// Then prob_diag, prob_non_diag are obtained by solving the above two equations.
	// c1 = sum_i pi_i L_ii
	// c2 = sum_i pi_i sum_(j n= i) L_ij

	double c1 = 0;
	for (int i=0; i<k; i++)
		c1 = c1 + p_prior.getProb(i) *
		p_loss_ptr->get_loss(p_space.getx(i), p_space.gety(i),p_space.getx(i), p_space.gety(i));

	double c2 = 0;
	for (int i=0; i<k; i++) {
		double sum=0;
		for (int j=0; j<k; j++)
			if (j!=i)
				sum = sum + p_loss_ptr->get_loss( p_space.getx(i), p_space.gety(i), p_space.getx(j), p_space.gety(j) );
		c2 = c2 + p_prior.getProb(i) * sum;
	}

	prob_non_diag = (c1 - p_exp_loss)/( (k-1)*c1 - c2 );
	prob_diag = 1 - (k-1) * prob_non_diag;

	// special case: for uniform prior and with L(i,i)=0, we have exp_loss = 1/k * non_diag_prob * sumL(i,j). Then
	//prob_non_diag = (k*p_exp_loss)/p_loss.get_loss_matrix(hsize, vsize, p_space.cell_side).sum();
	//prob_diag = 1 - (k-1)*prob_non_diag;

	if (prob_diag <= 0) {
		std::cerr << "Unable to construct RR mechanism: the given loss " <<p_exp_loss << " is too large" ;
		exit(0);
	}
	std::cout << "KRR with loss = " << p_exp_loss << ", eps = " << get_eps() << "\n";
}

PRRespMechanism::PRRespMechanism(const std::string p_info)
: GridOutputLocalMechanism(0,0){
	std::string mech_type_str = AbstractOutputLocalMechanism::mechTypeToString(KRR);
	std::istringstream iss(p_info);
	std::string p_mech_type_str; double p_ldp_eps;
	if (iss >> p_mech_type_str >> hsize >> vsize >> p_ldp_eps) {
		if (p_mech_type_str != mech_type_str) {
			std::cerr << "unable to construct " << mech_type_str << " mechanism from info string:\n "
					<< "info string must start with " << mech_type_str;
			exit(0);
		}
		// construction work
		double t= exp(p_ldp_eps);
		prob_diag = t/(hsize*vsize-1+t);
		prob_non_diag = 1.0/(hsize*vsize-1+t);
	}
	else {
		std::cerr << "unable to parse the mechanism info to construct "<< mech_type_str <<" mechanism";
	}
}

/** obfuscate the given location p_loc (in the linear order) into an observation
 * This overrides GridOutputMechanism to improve the performance, taking advantage
 * of the uniform distribution over hsize*vsize -1 elements */
GridObservation* PRRespMechanism::obfuscateLoc(int p_loc) const {
	//todo: debug

	//std::cout << "obfuscating loc = "<<p_loc << "\n";
	DoubleRow bin_distribution(2);
	bin_distribution.put(0, prob_diag);
	bin_distribution.put(1, (hsize*vsize -1)*prob_non_diag);
	//std::cout << "binary distribution = " << bin_distribution.toString() << "\n";
	std::vector<int> bin = Distributions::sampleDiscrete(bin_distribution, 1);
	//std::cout << "returned bin = " << bin[0] << "\n";
	if (bin[0]==0) {
		//std::cout << "returning " << p_loc << "\n";
		return new GridObservation(p_loc);
	}
	else {
		int s = Distributions::sampleUniformInt(0, hsize*vsize-2);
		if (s>=p_loc) s++;
		//std::cout << "returning " << s << "\n";
		//exit(0);
		return new GridObservation(s);
	}
}

double PRRespMechanism::getProb(int p_in, int p_out) const {
	if (p_in == p_out ) return prob_diag;
	else return prob_non_diag;
}

double PRRespMechanism::get_eps() const {
	return log(prob_diag/prob_non_diag);
}

/** get the basic info string of the mechanism */
std::string PRRespMechanism::get_info() const {
	double eps = log(prob_diag/prob_non_diag);
	return AbstractOutputLocalMechanism::mechTypeToString(KRR)+ " "+ std::to_string(hsize)+" "+
			std::to_string(vsize) +" "+std::to_string(eps);
}

/** get the type and privacy of the mechanism */
std::string PRRespMechanism::getTypeAndPriv(const int p_precision) const {
	double eps = log(prob_diag/prob_non_diag);
	return AbstractOutputLocalMechanism::mechTypeToString(KRR)+
	GSettings::to_string_with_precision(eps, p_precision);

}
