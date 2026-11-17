/*
 * TTaxicabMechanism.h
 *
 */

#ifndef PRIV_PLAN_MECHANISMS_TTAXICABMECHANISM_H_
#define PRIV_PLAN_MECHANISMS_TTAXICABMECHANISM_H_

/**
 * Truncated Approximate Planar Geometric mechanism
 * A mechanism is defined by two independent Truncated linear geometric mechanisms
 * randomizing on x and y dimensions respectively: P(out|in) = P(out_x|in_x) * P(out_y|in_y).
 * This mechanism satisfies d-privacy on a grid, where d is the manhattan distance.
 * This mechanism uses parameter plambda and satisfies (sqrt(2)*plambda)-geo-indistinguishability.
 * This mechanism is also the planar tight-constraints mechanism with respect to Manhattan privacy.
 * This mechanism is also the  truncated planar geometric mechanism with respect to Manhattan privacy.
 *  */

#include "GridOutputLocalMechanism.h"

class Loss;

class TTaxicabMechanism: public GridOutputLocalMechanism {
public:
	double cell_side;
	double lambda;

	/* constructs the mechanism on a full grid space p_space, with
	 * parameter p_lambda and p_cell_side. Note that this class has a cell_side attribute
	 * since the mechanism probabilities are determined by the distance between inputs and outputs.
	 * This is different from e.g. TCRectMechanism which is generic with respect to privacy metric
	 * and therefore is constructed by privMetric parameter.*/
	TTaxicabMechanism(int p_hsize, int p_vsize, double p_cell_side, double p_lambda);

	TTaxicabMechanism(const FullGridSpace& p_space, double p_lambda);

	TTaxicabMechanism(const std::string p_info);
	/**
	 * Constructs a mechanism that satisfies a given loss value for the given prior
	 * (with respect to the given loss object). This method is used to investigate mechanisms
	 * that satisfy a fixed user utility.
	 */
	TTaxicabMechanism(const FullGridSpace& p_space, const Prior& p_prior, const Loss* p_loss_ptr, double p_exp_loss);

	/** get the basic info string of the mechanism */
	std::string get_info() const override;
	std::string getTypeAndPriv(const int p_precision) const override;

private:

	/** returns the probability of an input and output points (specified in the linear order) */
	double getProb(int p_in, int p_out) const override;

};




#endif /* PRIV_PLAN_MECHANISMS_TTAXICABMECHANISM_H_ */
