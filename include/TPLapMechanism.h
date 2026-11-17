/*
 * TPLapMechansim.h
 *
 *  Created on: Jul 9, 2019
 *      
 */

#ifndef PRIV_PLAN_MECHANISMS_TPLAPMECHANISM_H_
#define PRIV_PLAN_MECHANISMS_TPLAPMECHANISM_H_

/**
 * Truncated planar Laplacian mechanism
 * The rows and columns of the mechanism matrix represent respectively the inputs and outputs
 * of the mechanism, where rows/columns are in the S order
 * i.e. the index of (x,y) is y*hsize+x.
 *
 * The mechanism is constructed by Monte Carlo simulation, i.e. by generating points according
 * to the pdf for every input and compute the probability of every cell as the proportion of
 * the inclosed points relative to the total generated points.
 */

#include "GridOutputLocalMechanism.h"
#include "DoubleMatrix.h"

class TPLapMechanism : public GridOutputLocalMechanism {

public:
	double cell_side;
	double eps; // used to construct Laplace matrix

	/**
	 * Constructs a TLap mechanism on a full grid space p_space (hsize, vsize, cell_side) and
	 * geo-indistinguishability parameter eps. Note that this mechanism must have cell_side
	 * attributes because the mechanism probabilities are determined by the distance between
	 * input and output cells.
	 */
	TPLapMechanism(int p_hsize, int p_vsize, double p_cell_side, double p_eps);

	/**
	 * Constructs a TLap mechanism on a full grid space p_space and geo-indistinguishability parameter eps.
	 * This mechanism also has cell_side attributes because the mechanism probabilities are determined
	 * by the distance between input and output cells.
	 */
	TPLapMechanism(const FullGridSpace& p_space, double p_eps);

	/**
	 * Constructs a mechanism that satisfies a given loss value for the uniform prior
	 * (with respect to the given loss object). This method is used to compare mechanisms with respect
	 * to statistical utility (using EM algorithm). Important: we assume that loss(i,i)=0.
	 */
	TPLapMechanism(const FullGridSpace& p_space, const Loss& p_loss, double p_exp_loss);

	TPLapMechanism(const std::string p_info);
	/**
	 * returns the mechanism matrix. Since this mechanism has already the mechanism matrix,
	 * it overrides the member function in the super class to improve performance
	 * */
	DoubleMatrix get_mech_matrix() const override;

	/**
	 * returns a mechanism row distribution (probabilities) at p_in. Since this mechanism
	 * has already the mechanism matrix, it overrides the member function in the super class.
	 * This method is used by EmpData class in the EM framework to generate obfuscated data.
	 * p_in is specified in linear order (y*hsize + x)
	 */
	DoubleRow get_mech_row(int p_in) const override;

	std::string get_info() const override;
	std::string getTypeAndPriv(const int p_precision) const override;

private:

	double getProb(int p_in, int p_out) const override;

	/**
	 * returns a mechanism column (with index p_out) as a column vector: Since this mechanism has already
	 * the mechanism matrix, it overrides the member function in the super class. p_out is
	 * specified in linear order (y*hsize + x).
	 */
	DoubleColumn get_mech_column(int p_out) const override;

	DoubleMatrix mechanism_matrix;  // stores the entries of the mechanism since every entry takes long time to be computed

	const int num_sampled_points_per_input; //  = 5000;  // number of points generated for each input cell.
	double inc_prob; //= (double)1/max_points;

	/** fill the mechanism matrix elements with Laplace values */
	void fill_matrix_lap();

	/**
	 * compute the mechanism entries having the input cell (p_in_x, p_in_y). This is performed by
	 * simulation: a large set of points are generated from (p_in_x, p_in_y), and then the conditional
	 * probability of every output cell in the gris is computed
	 * */
	void update_mech_matrix_for_input(int p_in_x,int p_in_y);

	double get_random();


};




#endif /* PRIV_PLAN_MECHANISMS_TPLAPMECHANISM_H_ */
