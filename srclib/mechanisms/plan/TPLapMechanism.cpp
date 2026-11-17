
/**
 * Truncated planar Laplace mechanism
 * The rows and columns of the mechanism matrix represent respectively the inputs and outputs
 * of the mechanism, where rows/columns are in the S order
 * i.e. the index of (x,y) is y*hsize+x.
 *
 * The mechanism is constructed by Monte Carlo simulation, i.e. by generating points according
 * to the pdf for every input and compute the probability of every cell as the proportion of
 * the inclosed points relative to the total generated points.
 */

#include <iostream>
#include <sstream>
#include <cmath>
#include "Distributions.h"
#include "DoubleRow.h"
#include "DoubleColumn.h"
#include "Loss.h"
#include "FullGridSpace.h"
#include "TPLapMechanism.h"
#include "GSettings.h"


/**
 * Constructs a TLap mechanism on a full grid space p_space (hsize, vsize, cell_side) and
 * geo-indistinguishability parameter eps. Note that this mechanism must have cell_side
 * attributes because the mechanism probabilities are determined by the distance between
 * input and output cells.
 */
TPLapMechanism::TPLapMechanism(int p_hsize, int p_vsize, double p_cell_side, double p_eps)
:GridOutputLocalMechanism(p_hsize, p_vsize), cell_side(p_cell_side), eps(p_eps),
 mechanism_matrix(p_hsize*p_vsize,p_hsize*p_vsize),
 num_sampled_points_per_input(30000) {
	inc_prob = (double)1/num_sampled_points_per_input;
	fill_matrix_lap();
}


/**
 * Constructs a TLap mechanism on a full grid space p_space and geo-indistinguishability parameter eps.
 * This mechanism also has cell_side attributes because the mechanism probabilities are determined
 * by the distance between input and output cells.
 */
TPLapMechanism::TPLapMechanism(const FullGridSpace& p_space, double p_eps)
:GridOutputLocalMechanism(p_space.hsize, p_space.vsize), cell_side(p_space.cell_side), eps(p_eps),
 mechanism_matrix(p_space.hsize*p_space.vsize,p_space.hsize*p_space.vsize),
 num_sampled_points_per_input(100000) {

	inc_prob = (double)1/num_sampled_points_per_input;
	fill_matrix_lap();
}

/**
 * Constructs a mechanism that satisfies a given loss value for the uniform prior
 * (with respect to the given loss object). This method is used to compare mechanisms with respect
 * to statistical utility (using EM algorithm). Important: we assume that loss(i,i)=0.
 */
TPLapMechanism::TPLapMechanism(const FullGridSpace& p_space, const Loss& p_loss, double p_exp_loss)
: GridOutputLocalMechanism(p_space.hsize, p_space.vsize), cell_side(p_space.cell_side),
  mechanism_matrix(p_space.hsize*p_space.vsize,p_space.hsize*p_space.vsize),
  num_sampled_points_per_input(100000) {

	double eps0=0.1; double eps1= 2/p_exp_loss;   // for p_exp_loss = 0.6, use 1.45, 1.0->0.85
	double u;
	do {
		eps = (eps0+eps1)/2;
		fill_matrix_lap();
		// compute utility
		DoubleMatrix loss_mat = p_loss.get_loss_matrix(hsize,vsize,cell_side); // get the loss matrix
		u = mechanism_matrix.muli(loss_mat).sum()/(hsize*vsize);
		//u = this.get_exp_loss_for_prior(new Prior(hsize,vsize,cell_side),p_loss  );
		if (p_exp_loss<u) eps0=eps; else eps1= eps;
	}
	while (abs(p_exp_loss - u)>0.001);
}

TPLapMechanism::TPLapMechanism(const std::string p_info)
:GridOutputLocalMechanism(0,0), mechanism_matrix(0,0), num_sampled_points_per_input(100000) {
	std::string mech_type_str = AbstractOutputLocalMechanism::mechTypeToString(PLAP);
	std::istringstream iss(p_info);
	std::string p_mech_type_str;
	if (iss >> p_mech_type_str >> hsize >> vsize >> cell_side >> eps) {
		if (p_mech_type_str != mech_type_str) {
			std::cerr << "unable to construct " << mech_type_str << " mechanism from info string:\n "
					<< "info string must start with " << mech_type_str;
			exit(0);
		}
		// construction work
		mechanism_matrix = DoubleMatrix( hsize*vsize, hsize*vsize );
		//num_sampled_points_per_input = 100000;
		inc_prob = (double)1/num_sampled_points_per_input;
		fill_matrix_lap();
	}
	else {
		std::cerr << "unable to parse the mechanism info to construct "<< mech_type_str <<" mechanism";
		exit(0);
	}
}


/**
 * returns the mechanism matrix. Since this mechanism has already the mechanism matrix,
 * it overrides the member function in the super class to improve performance
 * */
DoubleMatrix TPLapMechanism::get_mech_matrix() const {
	return mechanism_matrix;
}

/**
 * returns a mechanism row distribution (probabilities) at p_in. Since this mechanism
 * has already the mechanism matrix, it overrides the member function in the super class.
 * This method is used by EmpData class in the EM framework to generate obfuscated data.
 * p_in is specified in linear order (y*hsize + x)
 */
DoubleRow TPLapMechanism::get_mech_row(int p_in) const {
	return mechanism_matrix.selectRow(p_in);
}

/**
 * returns a mechanism column (with index p_out) as a column vector: Since this mechanism has already
 * the mechanism matrix, it overrides the member function in the super class. p_out is
 * specified in linear order (y*hsize + x).
 */
DoubleColumn TPLapMechanism::get_mech_column(int p_out) const {
	return mechanism_matrix.selectColumn(p_out);
}

double TPLapMechanism::getProb(int p_in, int p_out) const {
	return mechanism_matrix.get(p_in, p_out);
}



/** fill the mechanism matrix elements with Laplace values */
void TPLapMechanism::fill_matrix_lap() {
	mechanism_matrix.fillWith(0);   // reset the matrix
	for (int in_x=0; in_x<hsize;in_x++)
		for (int in_y=0; in_y<vsize;in_y++)
			update_mech_matrix_for_input(in_x,in_y);
}

/**
 * compute the mechanism entries having the input cell (p_in_x, p_in_y). This is performed by
 * simulation: a large set of points are generated from (p_in_x, p_in_y), and then the conditional
 * probability of every output cell in the gris is computed
 * */
void TPLapMechanism::update_mech_matrix_for_input(int p_in_x,int p_in_y) {
	double r;
	double theta;
	double dx,dy;
	double out_x_d,out_y_d; // horizontal, vertical displacements from the corner of the grid.
	int out_x,out_y;

	// lambert function (-1) evaluation variables
	double p;
	double z;
	double w;
	double w1;
	double l1,l2;
	for (int i=0;i<num_sampled_points_per_input;i++) {  // sample max_points

		// get a random noise vector
		theta = 2*M_PI*  get_random();  // get the angle
		p = get_random();
		z = (p-1)/M_E;
		l1= log(-z);
		l2= log(-log(-z));
		w = l1- l2+ l2/l1 + l2*(-2+l2)/(2*l1*l1) +
				l2*(6-9*l2+2*l2*l2)/(6*l1*l1*l1) +
				l2*(-12+36*l2-22*l2*l2+3*l2*l2*l2)/(12*l1*l1*l1*l1)+
				l2*(60-300*l2+350*l2*l2-125*l2*l2*l2+12*l2*l2*l2*l2)/(60*l1*l1*l1*l1*l1);
		// enhancing approximation of lambert function using Halley's iterative method,:
		for (int j=0;j<1;j++) {
			double num = w*exp(w)-z;
			double dom = exp(w)*(w+1) - (w+2)*(w*exp(w)-z)/(2*w+2);
			w1 = w - num/dom;
			w=w1;
		}
		r = (double)(-1)*(w+1)/eps;  // get the magnitude in km. Should be converted to orders of cell size
		// by dividing by cell_side

		dx = r*cos(theta)/cell_side;
		dy = r*sin(theta)/cell_side;
		out_x_d = 0.5+ p_in_x + dx;
		out_y_d = 0.5+ p_in_y + dy;
		out_x = (int)floor(out_x_d);
		out_y = (int)floor(out_y_d);

		// truncation
		if (out_x<0) out_x=0;
		if (out_x>hsize-1) out_x=hsize-1;
		if (out_y<0) out_y=0;
		if (out_y>vsize-1) out_y=vsize-1;

		// update the appropriate mechanism matrix element
		double t = mechanism_matrix.get(p_in_y*hsize+p_in_x,out_y*hsize+out_x);
		mechanism_matrix.put(p_in_y*hsize+p_in_x,out_y*hsize+out_x, t+inc_prob) ;
	}
}


double TPLapMechanism::get_random() {
	//todo: compare between the following two
	return (double) rand()/RAND_MAX;
	//return Distributions::sampleUniformDouble(0,1);
}

std::string TPLapMechanism::get_info() const {
	return AbstractOutputLocalMechanism::mechTypeToString(PLAP) + " " + std::to_string(hsize) + " "+
			std::to_string(vsize) + " "+ std::to_string(cell_side) + " "+ std::to_string(eps);
}

std::string TPLapMechanism::getTypeAndPriv(const int p_precision) const {
	return AbstractOutputLocalMechanism::mechTypeToString(PLAP) + 
	GSettings::to_string_with_precision(eps, p_precision);
}
