#include <iostream>
#include <string>
#include <sstream>    // for reading from string as a file
#include <cmath>
#include "DoubleMatrix.h"
#include "Loss.h"
#include "TVLLapMechanism.h"
#include "FullGridSpace.h"
#include "GSettings.h"
// A linear truncated Laplace mechanism on a vertical line of cells: hsize = 1, vsize= size

/*
 * constructs Truncated linear Laplace mechanism which works on a linear space (hsize=1), with
 * (\lambda*d)-privacy (lambda is the distinguishability per distance unit). Note that this
 * must have a cell_side attribute since mechanism probabilities are determined by the distance
 * between inputs and outputs.
 */
TVLLapMechanism::TVLLapMechanism(int p_size , double p_cell_side, double p_lambda)
:GridOutputLocalMechanism(1, p_size), cell_side(p_cell_side), lambda(p_lambda) {
}

TVLLapMechanism::TVLLapMechanism(const FullGridSpace& p_space, double p_lambda)
:GridOutputLocalMechanism(p_space.hsize, p_space.vsize), cell_side(p_space.cell_side), lambda(p_lambda) {
	if (p_space.hsize !=1) {
		std::cerr << "unable to construct a linear truncated Laplace mechanism: hsize must be 1";
		exit(0);

	}
}

/** Constructs a mechanism that satisfies a given loss value for the given prior (with respect to given loss)*/
TVLLapMechanism::TVLLapMechanism(int p_size, double p_cell_side, const Prior& p_prior, const Loss* p_loss_ptr, double p_exp_loss)
: GridOutputLocalMechanism(1, p_size), cell_side(p_cell_side) {
	double lambda0 = 0.01; double lambda1 = 20;    // initial bounds for eps
	double u;
	do {
		lambda = (lambda0+lambda1)/2;
		// compute utility
		u = get_exp_loss_for_prior(p_prior, p_loss_ptr, cell_side);
		if (p_exp_loss<u) lambda0=lambda; else lambda1= lambda;
	}
	while ( abs(p_exp_loss - u)>0.001 );
	std::cout << "VLLAP with loss = " << p_exp_loss << ", lambda = " << lambda << "\n";
}



TVLLapMechanism::TVLLapMechanism(const std::string p_info)
:GridOutputLocalMechanism(0,0) {  //base class has to be initialized here since there is no default constructor. members are set later
	std::string mech_type_str = AbstractOutputLocalMechanism::mechTypeToString(VLLAP);
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
 *  returns the entry of the truncated laplace mechanism for input p_in and output p_out.
 */
double TVLLapMechanism::getProb(int p_in, int p_out) const {
	// Laplace pdf f(x,z) = lambda/2 e^(-eps*|z-x|)
	// conditional probability of every line segment is the integral of the pdf over the segment.
	double p;
	double alpha = exp(-lambda*cell_side);
	int in=p_in; int out=p_out;
	if (out<in) {
		in=vsize-1-in;    // using the fact that Laplace is symmetric
		out=vsize-1-out;  // now we are sure that in<out
	}
	if (out==in)
		p=1-sqrt(alpha);
	else
		p=0.5*((1-alpha)/sqrt(alpha))*pow(alpha,out-in);
	if (out==vsize-1)
		p= p + 0.5*sqrt(alpha)*pow(alpha,out-in);  // truncation
	return p;
}

double TVLLapMechanism::get_lambda() {
	return lambda;
}

std::string TVLLapMechanism::get_info() const {
	return Mechanism::mechTypeToString(VLLAP) + " " + std::to_string(vsize) + " " +
			std::to_string(cell_side) + " "+ std::to_string(lambda);
}
std::string TVLLapMechanism::getTypeAndPriv(const int p_precision) const {
	return Mechanism::mechTypeToString(VLLAP) + 
	GSettings::to_string_with_precision(lambda, p_precision);
}
