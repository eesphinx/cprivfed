/**
 * Truncated Taxicab mechanism
 * A mechanism is defined by two independent Truncated linear geometric mechanisms
 * randomizing on x and y dimensions respectively: P(out|in) = P(out_x|in_x) * P(out_y|in_y).
 * This mechanism satisfies d-privacy on a grid, where d is the manhattan distance.
 * This mechanism uses parameter plambda and satisfies (sqrt(2)*plambda)-geo-indistinguishability.
 * This mechanism is also the planar tight-constraints mechanism with respect to Manhattan privacy.
 * This mechanism is also the  truncated planar geometric mechanism with respect to Manhattan privacy.
 *  */

#include <iostream>
#include <sstream>
#include <cmath>
#include "DoubleMatrix.h"
#include "Loss.h"
#include "FullGridSpace.h"
#include "TTaxicabMechanism.h"
#include "GSettings.h"

	/* constructs the mechanism on a full grid space p_space, with
	 * parameter p_lambda and p_cell_side. Note that this class has a cell_side attribute
	 * since the mechanism probabilities are determined by the distance between inputs and outputs.
	 * This is different from e.g. TCRectMechanism which is generic with respect to privacy metric
	 * and therefore is constructed by privMetric parameter.*/
TTaxicabMechanism::TTaxicabMechanism(int p_hsize, int p_vsize, double p_cell_side, double p_lambda)
:GridOutputLocalMechanism(p_hsize, p_vsize), cell_side(p_cell_side), lambda(p_lambda) {
}

TTaxicabMechanism::TTaxicabMechanism(const FullGridSpace& p_space, double p_lambda)
: GridOutputLocalMechanism(p_space.hsize, p_space.vsize), cell_side(p_space.cell_side), lambda(p_lambda) {
}

/**
 * Constructs a mechanism that satisfies a given loss value for the given prior
 * with respect to the given loss object. This method is used to investigate mechanisms
 * that satisfy a fixed user utility.
 */
TTaxicabMechanism::TTaxicabMechanism(const FullGridSpace& p_space, const Prior& p_prior, const Loss* p_loss_ptr, double p_exp_loss)
: GridOutputLocalMechanism(p_space.hsize, p_space.vsize), cell_side(p_space.cell_side) {
	double lambda0 = 0.001; double lambda1 = 20;    // initial bounds for eps
	double u;
	do {
		lambda = (lambda0+lambda1)/2;
		// compute utility
		u = get_exp_loss_for_prior(p_prior, p_loss_ptr, p_space.cell_side);
		//std::cout << "Geometric: " << "eps = " << eps << ", u = " << u << std::endl;
		if (p_exp_loss<u) lambda0=lambda; else lambda1= lambda;
	}
	while ( abs(p_exp_loss - u)>0.001 );
	std::cout << "Taxicab with loss = " << p_exp_loss << ", lambda = " << lambda << "\n";
}

TTaxicabMechanism::TTaxicabMechanism(const std::string p_info)
:GridOutputLocalMechanism(0,0) {
	std::string mech_type_str = AbstractOutputLocalMechanism::mechTypeToString(TAXICAB);
	std::istringstream iss(p_info);
	std::string p_mech_type_str; int p_hsize, p_vsize; double p_cell_side; double p_lambda;
	if (iss >> p_mech_type_str >> p_hsize >> p_vsize >> p_cell_side >> p_lambda) {
		if (p_mech_type_str != mech_type_str) {
			std::cerr << "unable to construct " << mech_type_str << " mechanism from info string:\n "
					<< "info string must start with " << mech_type_str;
			exit(0);
		}
		// construction work
		hsize = p_hsize; vsize = p_vsize; cell_side = p_cell_side; lambda = p_lambda;
	}
	else {
		std::cerr << "unable to parse the mechanism info to construct "<< mech_type_str <<" mechanism";
		exit(0);
	}
}


/** returns the probability of an input and output cells (specified in the linear order) */
double TTaxicabMechanism::getProb(int p_in, int p_out) const {
	int p_in_x = p_in%hsize; int p_in_y = p_in/hsize;
	int p_out_x = p_out%hsize; int p_out_y = p_out/hsize;

	double alpha = exp(-lambda*cell_side);
	double prob_x,prob_y;

	if (p_out_x==0 || p_out_x==hsize-1)
		prob_x= (double)1/(1+alpha);
	else prob_x = (double) (1-alpha)/(1+alpha);
	prob_x = prob_x * pow(alpha, abs(p_out_x - p_in_x));

	if (p_out_y==0 || p_out_y==vsize-1)
		prob_y= (double)1/(1+alpha);
	else prob_y = (double) (1-alpha)/(1+alpha);
	prob_y = prob_y * pow(alpha, abs(p_out_y - p_in_y));

	return prob_x*prob_y;

}


std::string TTaxicabMechanism::get_info() const {
	return AbstractOutputLocalMechanism::mechTypeToString(TAXICAB) + " " + std::to_string(hsize) + " "+
			std::to_string(vsize) + " "+ std::to_string(cell_side) + " "+ std::to_string(lambda);
}

std::string TTaxicabMechanism::getTypeAndPriv(const int p_precision) const {
	return AbstractOutputLocalMechanism::mechTypeToString(TAXICAB) + 
	GSettings::to_string_with_precision(lambda, p_precision);
}
