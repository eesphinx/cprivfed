/*
 * TLGeomMechanism.h
 *
 *  Created on: Jun 21, 2019
 *      
 */

/*
 * A linear truncated geometric mechanism on a vertically linear space p_space (hsize=1),
 * with (lambda*d)-privacy with distinguishability lambda (for every unit distance). Note that this
 * has a cell_side attribute since mechanism probabilities are determined by the distance between
 * inputs and outputs.
 */


#ifndef MECHANISMS_LINE_TVLGEOMMECHANISM_H_
#define MECHANISMS_LINE_TVLGEOMMECHANISM_H_

#include "GridOutputLocalMechanism.h"     // we have to include because it is the base class

class TVLGeomMechanism : public GridOutputLocalMechanism {

public:

	TVLGeomMechanism(int p_vsize , double p_cell_side, double p_lambda);

	TVLGeomMechanism(const FullGridSpace& p_space, double p_lambda);

	/** Constructs a mechanism that satisfies a given loss value for the given prior (with respect to given loss)*/
	TVLGeomMechanism(int p_size, double p_cell_side, const Prior& p_prior, const Loss* p_loss_ptr, double p_exp_loss);

	TVLGeomMechanism(const std::string p_info);

	double get_lambda();

	std::string get_info() const override;
	std::string getTypeAndPriv(const int p_precision) const override;

private:
	double cell_side;
	double lambda;

	/**
	 *  returns the entry of the truncated geometric mechanism for input p_in and output p_out.
	 */
	double getProb(int p_in, int p_out) const override;
};




#endif /* MECHANISMS_LINE_TVLGEOMMECHANISM_H_ */
