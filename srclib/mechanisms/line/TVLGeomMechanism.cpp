#include <iostream>
#include <string>
#include <sstream>    // for reading from string as a file
#include <cmath>
#include "DoubleMatrix.h"
#include "Loss.h"
#include "TVLGeomMechanism.h"
#include "FullGridSpace.h"
#include "GSettings.h"
// A linear truncated geometric mechanism on a vertical line of cells: hsize = 1, vsize= size

/* constructs a linear truncated geometric mechanism on a vertically linear space p_space (hsize=1),
 * with (lambda*d)-privacy with distinguishability lambda (for every unit distance). Note that this
 * has a cell_side attribute since mechanism probabilities are determined by the distance between
 * inputs and outputs.
 */
TVLGeomMechanism::TVLGeomMechanism(int p_size , double p_cell_side, double p_lambda)
:GridOutputLocalMechanism(1, p_size), cell_side(p_cell_side), lambda(p_lambda) {
}

TVLGeomMechanism::TVLGeomMechanism(const FullGridSpace& p_space, double p_lambda)
:GridOutputLocalMechanism(p_space.hsize, p_space.vsize), cell_side(p_space.cell_side), lambda(p_lambda) {
	if (p_space.hsize !=1) {
		std::cerr << "unable to construct a linear truncated geometric mechanism: hsize must be 1";
		exit(0);

	}
}

/** Constructs a mechanism that satisfies a given loss value for the given prior (with respect to given loss)*/
TVLGeomMechanism::TVLGeomMechanism(int p_size, double p_cell_side, const Prior& p_prior, const Loss* p_loss_ptr, double p_exp_loss)
: GridOutputLocalMechanism(1, p_size), cell_side(p_cell_side) {
	double lambda0 = 0.01; double lambda1 = 20;    // initial bounds for lambda
	double u;
	do {
		lambda = (lambda0+lambda1)/2;
		// compute utility
		u = get_exp_loss_for_prior(p_prior, p_loss_ptr, cell_side);
		if (p_exp_loss<u) lambda0=lambda; else lambda1= lambda;
	}
	while ( abs(p_exp_loss - u)>0.001 );
	std::cout << "VLGEOM with loss = " << p_exp_loss << ", lambda = " << lambda << "\n";
}



TVLGeomMechanism::TVLGeomMechanism(const std::string p_info)
:GridOutputLocalMechanism(0,0) {  //base class has to be initialized here since there is no default constructor. members are set later
	std::string mech_type_str = AbstractOutputLocalMechanism::mechTypeToString(VLGEOM);
	std::istringstream iss(p_info);
	std::string p_mech_type_str; int p_size; double p_cell_side; double p_lambda;
	if (iss >> p_mech_type_str >> p_size >> p_cell_side >> p_lambda) {
		if (p_mech_type_str != mech_type_str) {
			std::cerr << "unable to construct " << mech_type_str << " mechanism from info string:\n "
					<< "info string must start with " << mech_type_str;
			exit(0);
		}
		// construction work
		hsize = 1; vsize = p_size; cell_side = p_cell_side; lambda = p_lambda;
	}
	else {
		std::cerr << "unable to parse the mechanism info to construct "<< mech_type_str <<" mechanism";
		exit(0);
	}
}

/**
 *  returns the entry of the truncated geometric mechanism for input p_in and output p_out.
 */
double TVLGeomMechanism::getProb(int p_in, int p_out) const {
	//  In this implementation I use the fact that this mechanism is also a TC mechanism
	//  Therefore it evaluates the diagonal element c[out, out] and then it evaluates the
	//  required probability as p[in,out]=c[out, out]*alpha^d(in, out).

	// for debug
	//std::cout << "getting LGeom probability from "<< p_in << " to " << p_out;
	double alpha = exp(-lambda*cell_side);
	double prob;

	if (p_out==0 || p_out==vsize-1)
		prob= (double)1/(1+alpha);
	else prob = (double) (1-alpha)/(1+alpha);
	prob = prob * pow(alpha, abs(p_out - p_in));

	return prob;
}

double TVLGeomMechanism::get_lambda() {
	return lambda;
}

std::string TVLGeomMechanism::get_info() const {
	return Mechanism::mechTypeToString(VLGEOM) + " " + std::to_string(vsize) + " " +
			std::to_string(cell_side) + " "+ std::to_string(lambda);
}

std::string TVLGeomMechanism::getTypeAndPriv(const int p_precision) const {
	return Mechanism::mechTypeToString(VLGEOM) + 
	       GSettings::to_string_with_precision(lambda, p_precision);
}
