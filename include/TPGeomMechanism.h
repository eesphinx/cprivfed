/*
 * TPGeomMechanism.h
 *
 *  Created on: Jul 9, 2019
 *      
 */

#ifndef PRIV_PLAN_MECHANISMS_TPGEOMMECHANISM_H_
#define PRIV_PLAN_MECHANISMS_TPGEOMMECHANISM_H_

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

#include "GridOutputLocalMechanism.h"

class Loss;

class TPGeomMechanism: public GridOutputLocalMechanism {
public:
	double cell_side;
	double eps;     // d-privacy epsilon

	/* constructs a planar truncated geometric mechanism on a full grid space p_space, with
	 * geo-indistinguishability eps (for 1km). Note that this class has a cell_side attribute
	 * since mechanism probabilities are determined by the distance between inputs and outputs.
	 * This is different from e.g. TCRectMechanism which is generic with respect to privacy metric
	 * and therefore is constructed by privMetric parameter.*/
	TPGeomMechanism(int p_hsize, int p_vsize, double p_cell_side, double p_eps);

	TPGeomMechanism(const FullGridSpace& p_space, double p_eps);

	TPGeomMechanism(const std::string p_info);
	/**
	 * Constructs a mechanism that satisfies a given loss value for the given prior
	 * (with respect to the given loss object). This method is used to investigate mechanisms
	 * that satisfy a fixed user utility.
	 */
	TPGeomMechanism(const FullGridSpace& p_space, const Prior& p_prior, const Loss* p_loss_ptr, double p_exp_loss);

	/** destructor: to free the dynamically allocated border arrays */
	~TPGeomMechanism();


	/** get the basic info string of the mechanism */
	std::string get_info() const override;
	std::string getTypeAndPriv(const int p_precision) const override;
	

private:

	/**
	 * probabilities for outputs on the horizontal border in the form hborder_prob[in_x][in_y][out_x]
	 * hborder_prob[in_x][in_y][out_x] is the probability of producing the output (out_x,0) from the
	 * input (in_x, in_y)
	 **/
	double*** hborder_prob;

	/**
	 * probabilities for outputs on the vertical border in the form hborder_prob[in_x][in_y][out_y]
	 * vborder_prob[in_x][in_y][out_y] is the probability of producing the output (0,out_y) from the
	 * input (in_x, in_y)
	 **/
	double*** vborder_prob;


	/** The normalization constant, also is probability that out = in */
	double g;
	const double vsmall;

	/** returns the probability of an input and output points (specified in the linear order) */
	double getProb(int p_in, int p_out) const override;

	/** computes the normalization constant */
	void compute_g();

	/** fill the mechanism matrix elements */
	void fill_border_probs();

	/** returns the mechanism probability for input in_x, in_y and output out_x, out_y */
	double get_tgeom_prob(int p_in_x,int p_in_y,int p_out_x,int p_out_y) const;


};




#endif /* PRIV_PLAN_MECHANISMS_TPGEOMMECHANISM_H_ */
