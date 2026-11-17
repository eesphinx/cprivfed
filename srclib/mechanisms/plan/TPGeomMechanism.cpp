/**
 * Truncated planar geometric mechanism on a rectangular grid.
 * The mechanism matrix has the form
 * (double) mechanism[i_x][i_y][j_x][j_y] where the first
 * coordinates represent the row (mechanism input) and the other two
 * coordinates represent the column (mechanism output).

 * Here we optimize the construction and sampling memory by using the symmetries.
 * Basically for every interior output (out_x, out_y), we do not need to store the
 * conditional probability of producing (out_x,out_y) from a given input (in_x, in_y)
 * because this probability is obtained easily as G.e^-eps*d(in,out). However, due to
 * the truncation, we need to store for every border output (out_x,0), and every input
 * (in_x, in_y) the probability of producing (out_x,0) given input (in_x, in_y). We
 * also need to do the same for vertical border outputs (0, out_y).
 */

#include <iostream>
#include <sstream>
#include <cmath>
#include "DoubleMatrix.h"
#include "Loss.h"
#include "FullGridSpace.h"
#include "TPGeomMechanism.h"
#include "GSettings.h"

	/* constructs a planar truncated geometric mechanism on a full grid space p_space, with
	 * geo-indistinguishability eps (for 1km). Note that this class has a cell_side attribute
	 * since mechanism probabilities are determined by the distance between inputs and outputs.
	 * This is different from e.g. TCRectMechanism which is generic with respect to privacy metric
	 * and therefore is constructed by privMetric parameter.*/
TPGeomMechanism::TPGeomMechanism(int p_hsize, int p_vsize, double p_cell_side, double p_eps)
:GridOutputLocalMechanism(p_hsize, p_vsize), cell_side(p_cell_side), eps(p_eps), vsmall(1e-13) {
	compute_g();
	fill_border_probs();
}

TPGeomMechanism::TPGeomMechanism(const FullGridSpace& p_space, double p_eps)
: GridOutputLocalMechanism(p_space.hsize, p_space.vsize), cell_side(p_space.cell_side), eps(p_eps), vsmall(1e-13) {
	//infinity = 200;
	compute_g();
	fill_border_probs();
}

/**
 * Constructs a mechanism that satisfies a given loss value for the given prior
 * with respect to the given loss object. This method is used to investigate mechanisms
 * that satisfy a fixed user utility.
 */
TPGeomMechanism::TPGeomMechanism(const FullGridSpace& p_space, const Prior& p_prior, const Loss* p_loss_ptr, double p_exp_loss)
: GridOutputLocalMechanism(p_space.hsize, p_space.vsize), cell_side(p_space.cell_side), vsmall(1e-10) {
	double eps0 = 0.01; double eps1 = 20;    // initial bounds for eps
	double u;
	do {
		eps = (eps0+eps1)/2;
		compute_g();
		fill_border_probs();
		// compute utility
		u = get_exp_loss_for_prior(p_prior, p_loss_ptr, p_space.cell_side);
		//std::cout << "Geometric: " << "eps = " << eps << ", u = " << u << std::endl;
		if (p_exp_loss<u) eps0=eps; else eps1= eps;
	}
	while ( abs(p_exp_loss - u)>0.001 );
	std::cout << "PGEOM with loss = " << p_exp_loss << ", eps = " << eps << "\n";
}

TPGeomMechanism::TPGeomMechanism(const std::string p_info)
:GridOutputLocalMechanism(0,0), vsmall(10e-10) {
	std::string mech_type_str = AbstractOutputLocalMechanism::mechTypeToString(PGEOM);
	std::istringstream iss(p_info);
	std::string p_mech_type_str; int p_hsize, p_vsize; double p_cell_side; double p_eps;
	if (iss >> p_mech_type_str >> p_hsize >> p_vsize >> p_cell_side >> p_eps) {
		if (p_mech_type_str != mech_type_str) {
			std::cerr << "unable to construct " << mech_type_str << " mechanism from info string:\n "
					<< "info string must start with " << mech_type_str;
			exit(0);
		}
		// construction work
		hsize = p_hsize; vsize = p_vsize; cell_side = p_cell_side; eps = p_eps;
		compute_g();
		fill_border_probs();
	}
	else {
		std::cerr << "unable to parse the mechanism info to construct "<< mech_type_str <<" mechanism";
		exit(0);
	}
}

/** destructor: frees the dynamically allocated border arrays */
TPGeomMechanism::~TPGeomMechanism() {
	if (hborder_prob) {   // this check is important to avoid deallocating data pointed to by null pointers
		for (int in_x=0; in_x < hsize; in_x++) {
			for (int in_y=0; in_y< vsize; in_y++) {
				delete [] hborder_prob[in_x][in_y];
				hborder_prob[in_x][in_y] = nullptr;
			}
			delete [] hborder_prob[in_x];
			hborder_prob[in_x] = nullptr;
		}
		delete [] hborder_prob;
		hborder_prob = nullptr;
	}
	if (vborder_prob) {
		for (int in_x=0; in_x < hsize; in_x++) {
			for (int in_y=0; in_y< vsize; in_y++) {
				delete [] vborder_prob[in_x][in_y];
				vborder_prob[in_x][in_y] = nullptr;
			}
			delete [] vborder_prob[in_x];
			vborder_prob[in_x] = nullptr;
		}
		delete [] vborder_prob;
		vborder_prob = nullptr;
	}
}


/** returns the probability of an input and output cells (specified in the linear order) */
double TPGeomMechanism::getProb(int p_in, int p_out) const {
	int p_in_x = p_in%hsize; int p_in_y = p_in/hsize;
	int p_out_x = p_out%hsize; int p_out_y = p_out/hsize;
	double ret=0;
	// check if the output is on the horizontal border, i.e. with p_out_y=0 or p_out_y=vsize-1
	if (p_out_y==0)
		ret=hborder_prob[p_in_x][p_in_y][p_out_x];
	else if (p_out_y==vsize-1)
		ret=hborder_prob[p_in_x][vsize-1-p_in_y][p_out_x];
	// check if the output is on the vertical border, i.e. with p_out_x=0 or p_out_x=hsize-1
	else if (p_out_x==0)
		ret=vborder_prob[p_in_x][p_in_y][p_out_y];
	else if (p_out_x==hsize-1)
		ret=vborder_prob[hsize-1-p_in_x][p_in_y][p_out_y];
	// otherwise the output is in the interior
	else
		ret = get_tgeom_prob(p_in_x, p_in_y, p_out_x, p_out_y);
	return ret;
}

/** computes the normalization constant */
void TPGeomMechanism::compute_g() {
	// Compute the normalization constant G = sum_{(i,j) in Z^2}  e^{- d sqrt(i^2+j^2)}
	// sum terms:
	// i=0, j=0  --> G0 = 1;    // 1 time
	// i>0, j=0  --> G1 = sum_{i>0} e^{-eps*d*i} = e^-eps*d / (1 - e^-eps*d) // 4 times
	// i=j>0     --> G2 = sum_{i>0} e^{-eps*d*sqrt(2)*i} = e^-eps*d*sqrt(2)/ ( 1- e^-eps*d*sqrt(2)) // 4 times
	// i>0,i>j>0 --> G3 = sum_{i>0,i>j>0} e^{-eps*d*sqrt(i^2+j^2)}  // 8 times
	g = 1;
	g = g + 4 *(double) exp(- eps*cell_side)/(1-exp(- eps*cell_side));
	g = g + 4 *(double) exp(- eps*cell_side * sqrt(2.0))/(1- exp(- eps*cell_side* sqrt(2.0)));

	bool g_min_precision_reached = false;
	double prob;

	int i=2; int j=1;
	while (!g_min_precision_reached) {
		prob = 8* exp(- eps*cell_side * sqrt(i*i+j*j));
		if (prob < vsmall)   // assuming precision to 8 decimal places, i.e. 10^-8
			g_min_precision_reached = true;
		g = g + prob;
		j++;
		if (j==i) {
			i++;j=1;
		}
	}
	g = (double) 1/g;
}

/** fill the mechanism matrix elements */
void TPGeomMechanism::fill_border_probs() {
	// allocate the arrays
	hborder_prob = new double**[hsize];
	for (int in_x = 0; in_x < hsize; in_x++) {
		hborder_prob[in_x] = new double*[vsize];
		for (int in_y=0; in_y < vsize; in_y++)
			hborder_prob[in_x][in_y] = new double[hsize];
	}
	vborder_prob = new double**[hsize];
	for (int in_x = 0; in_x < hsize; in_x++) {
		vborder_prob[in_x] = new double*[vsize];
		for (int in_y=0; in_y < vsize; in_y++)
			vborder_prob[in_x][in_y] = new double[vsize];
	}

	// compute the entries of the hborder_prob
	for (int in_x=0; in_x<hsize;in_x++)
		for (int in_y=0; in_y<vsize;in_y++)
			for (int out_x=0;out_x<hsize;out_x++)
				hborder_prob[in_x][in_y][out_x] = get_tgeom_prob(in_x,in_y,out_x,0);

	// compute the entries of the vborder_prob
	for (int in_x=0; in_x<hsize;in_x++)
		for (int in_y=0; in_y<vsize;in_y++)
			for (int out_y=0;out_y<vsize;out_y++)
				vborder_prob[in_x][in_y][out_y] = get_tgeom_prob(in_x,in_y,0,out_y);
	//hborder_prob=new double[hsize][vsize][hsize];
	//vborder_prob=new double[hsize][vsize][vsize];
}

/** returns the mechanism probability for input in_x, in_y and output out_x, out_y
 * Note that in this particular function we pass x,y coordinates instead of the linear index
 * because it handle borders */
double TPGeomMechanism::get_tgeom_prob(int p_in_x,int p_in_y,int p_out_x,int p_out_y) const {

	int imax = std::numeric_limits<int>::max();  // max integer value

	int out_x_start= p_out_x; // starting x for summation.
	int out_x_end =  p_out_x; // end x for summation.
	int out_y_start= p_out_y; // starting y for summation.
	int out_y_end =  p_out_y; // ending y for summation.

	int xstep=1, ystep=1;  // specifies whether the summations goes forward or backward.

	// adjust the boundaries of the sum if out is on the borders
	if (p_out_x==0) {
		out_x_end = out_x_end-imax;
		xstep=-1;
	};
	if (p_out_x==hsize-1) {
		out_x_end = imax;
		xstep =1;
	}
	if (p_out_y==0) {
		out_y_end = out_y_end-imax;
		ystep = -1;
	}
	if (p_out_y==vsize-1) {
		out_y_end = imax;
		ystep = 1;
	}

	// now make the summation to compute the required probability
	double dist; // distance between input and output positions
	double sum=0;

	int out_x_curr, out_y_curr;
	double xprob;
	double yprob;
	bool x_min_precision_reached = false;
	bool y_min_precision_reached = false;

	out_x_curr=out_x_start;
	while ( ((xstep>0 && out_x_curr <= out_x_end) || (xstep< 0 && out_x_curr >= out_x_end))
			&& !x_min_precision_reached) {
		xprob=0;
		out_y_curr=out_y_start;
		y_min_precision_reached = false;
		while (((ystep>0 && out_y_curr <= out_y_end) || (ystep< 0 && out_y_curr >= out_y_end))
				&& !y_min_precision_reached) {
			dist= sqrt(pow(out_x_curr- p_in_x, 2)+ pow(out_y_curr-p_in_y, 2));
			yprob = g*exp(-eps*cell_side*dist);
			xprob = xprob + yprob;
			out_y_curr=out_y_curr + ystep;
			if (yprob<vsmall)
				y_min_precision_reached = true;
		}
		out_x_curr = out_x_curr+xstep;
		sum = sum + xprob;
		if (xprob<vsmall)
			x_min_precision_reached = true;
	}

	return sum;
}

std::string TPGeomMechanism::get_info() const {
	return Mechanism::mechTypeToString(PGEOM) + " " + std::to_string(hsize) + " "+
			std::to_string(vsize) + " "+ std::to_string(cell_side) + " "+ std::to_string(eps);
}

std::string TPGeomMechanism::getTypeAndPriv(const int p_precision) const {
	return Mechanism::mechTypeToString(PGEOM) + 
	GSettings::to_string_with_precision(eps, p_precision);
}
