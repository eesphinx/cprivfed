#include <iostream>
#include <cmath>
#include "Loss.h"
#include "DoubleMatrix.h"

/**
 * This class specifies the necessary function (get_loss) that computes
 * the loss based on the coordinates of given cells. Every loss function
 * is assumed to be increasing. This function behaves as follows.
 * loss = get_loss_by_dist(distance) For all 0<=distance<= saturation_dist,
 * loss = max_loss for all saturation_dist < dist.
 */
Loss::Loss(double p_saturation_dist, double p_max_loss)
:saturation_dist(p_saturation_dist), max_loss(p_max_loss){
	// the max_loss must be at least equal to get_loss(saturation_dist)
	//if (p_max_loss < get_loss_by_dist(p_saturation_dist)) {
	//	std::cerr << "unable to create loss function: specified max loss is too low";
	//	exit(0);
}

double Loss::get_saturation_dist() const {
	return saturation_dist;
}
double Loss::get_max_loss() const {
	return max_loss;
}

/**
 * for a grid with dimensions hsize,vsize, cell_side returns a matrix representing the loss for every pair
 * of input and output cells. We assume that loss_matrix is a matrix of dimensions
 * (p_hsize*p_vsize) * (p_hsize*p_vsize) and locations in row and column are ordered in S (linear) order
 **/

DoubleMatrix Loss::get_loss_matrix(int p_hsize, int p_vsize, double p_cell_side) const {
	DoubleMatrix loss_matrix(p_hsize*p_vsize, p_hsize*p_vsize);
	double loss;
	int in_x,in_y,out_x,out_y;
	for (int in=0; in <p_hsize*p_vsize; in++)
		for (int out=0; out<p_hsize*p_vsize; out++) {
			in_x = in % p_hsize;
			in_y = in / p_hsize;
			out_x = out % p_hsize;
			out_y = out / p_hsize;
			loss = get_loss(in_x*p_cell_side, in_y*p_cell_side, out_x*p_cell_side, out_y*p_cell_side);
			loss_matrix.put(in, out, loss);
		}
	return loss_matrix;
}
