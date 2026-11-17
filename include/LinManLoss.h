/*
 * LinManLoss.h
 *
 *  Created on: Jun 23, 2019
 *      
 */

#ifndef PRIV_PLAN_LOSS_LINMANLOSS_H_
#define PRIV_PLAN_LOSS_LINMANLOSS_H_

# include "Loss.h"

class LinManLoss: public Loss {
public:

	LinManLoss(double p_saturation_dist, double p_max_loss);

	double get_loss(double i_x, double i_y, double j_x, double j_y) const override;

	double get_loss_by_dist(double p_man_dist) const override;


	lossType getType() const override;


};





#endif /* PRIV_PLAN_LOSS_LINMANLOSS_H_ */
