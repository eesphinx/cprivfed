/*
 * BinLoss.h
 *
 *  Created on: Jun 23, 2019
 *      
 */

#ifndef PRIV_PLAN_LOSS_BINLOSS_H_
#define PRIV_PLAN_LOSS_BINLOSS_H_

# include "Loss.h"
// loss is 0 within distance saturation_dist, and is max_loss for larger distances
// note: in the special case if saturation_dist=0, the loss is 0 only when dist=0,
// and max_loss otherwise.
class BinLoss: public Loss {

public:

	BinLoss(double p_saturation_dist, double p_max_loss);

	double get_loss(double i_x, double i_y, double j_x, double j_y) const;

	double get_loss_by_dist(double p_dist) const;

	lossType get_name() const;
};




#endif /* PRIV_PLAN_LOSS_BINLOSS_H_ */
