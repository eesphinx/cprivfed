/*
 * TLLapMechanism.h
 *
 *  Created on: May 25, 2023
 *      
 */

/*
 * constructs Truncated linear Laplace mechanism which works on a linear space (hsize=1), with
 * (\lambda*d)-privacy (lambda is the distinguishability per distance unit). Note that this
 * must have a cell_side attribute since mechanism probabilities are determined by the distance
 * between inputs and outputs.
 */


#ifndef MECHANISMS_LINE_TVLLAPMECHANISM_H_
#define MECHANISMS_LINE_TVLLAPMECHANISM_H_

#include "GridOutputLocalMechanism.h"     // we have to include because it is the base class

class TVLLapMechanism : public GridOutputLocalMechanism {

public:

	TVLLapMechanism(int p_vsize , double p_cell_side, double p_lambda);

	TVLLapMechanism(const FullGridSpace& p_space, double p_lambda);

	/** Constructs a mechanism that satisfies a given loss value for the given prior (with respect to given loss)*/
	TVLLapMechanism(int p_size, double p_cell_side, const Prior& p_prior, const Loss* p_loss_ptr, double p_exp_loss);

	TVLLapMechanism(const std::string p_info);

	double get_lambda();

	std::string get_info() const override;
	std::string getTypeAndPriv(const int p_precision) const override;

private:
	double cell_side;
	double lambda;

	/**
	 *  returns the entry of the truncated laplace mechanism for input p_in and output p_out.
	 */
	double getProb(int p_in, int p_out) const override;
};




#endif /* MECHANISMS_LINE_TVLLAPMECHANISM_H_ */
