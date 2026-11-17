/*
 * LinGeoLoss.h
 *
 *  Created on: Jun 23, 2019
 *      
 */

#ifndef PRIV_PLAN_LOSS_LINGEOLOSS_H_
#define PRIV_PLAN_LOSS_LINGEOLOSS_H_

// package priv.plan.loss;

/** This class corresponds to the Euclidean distance between location.
 * @author ehab
 */

# include "Loss.h"


class LinGeoLoss: public Loss {
public:
	/**
	 * Constructs an Euclidean distance loss function, with a saturation distance at which the loss saturates to a fixed maximum loss.
	 * @param p_saturation_dist : Euclidean distance at which the loss saturates.
	 * @param p_max_loss : Maximum loss value
	 */
	LinGeoLoss(double p_saturation_dist, double p_max_loss);

	/**
	 * returns the loss given the real point is (i_x,i_y) and the noisy point is (j_x,j_y)
	 */
	double get_loss(double i_x, double i_y, double j_x, double j_y) const override;

	/**
	 * returns the loss between two points given the distance between them.
	 */
	double get_loss_by_dist(double p_geo_dist) const override;
	/**
	 * returns the identifier of this loss function
	 */
	lossType getType() const override;
};





#endif /* PRIV_PLAN_LOSS_LINGEOLOSS_H_ */
