//package priv.plan.loss;

/** This class corresponds to the Euclidean distance between location.
 * @author ehab
 */

#include <cmath>
#include "LinGeoLoss.h"

/**
 * Constructs an Euclidean distance loss function, with a saturation distance at which the loss saturates to a fixed maximum loss.
 * @param p_saturation_dist : Euclidean distance at which the loss saturates.
 * @param p_max_loss : Maximum loss value
 */
LinGeoLoss::LinGeoLoss(double p_saturation_dist, double p_max_loss)
:Loss(p_saturation_dist,p_max_loss){
}

	// set saturation distance, and max loss. In general they are independent
	// of each other.

/**
 * returns the loss given the real point is (i_x,i_y) and the noisy point is (j_x,j_y)
 */
double LinGeoLoss::get_loss(double i_x, double i_y, double j_x, double j_y) const {
	return get_loss_by_dist(sqrt(pow(i_x-j_x, 2) + pow(i_y-j_y, 2)));
}

/**
 * returns the loss between two points given the distance between them.
 */
double LinGeoLoss::get_loss_by_dist(double p_geo_dist) const {
	if (p_geo_dist<=saturation_dist)
		return p_geo_dist;
	else
		return max_loss;
}
/**
 * returns the identifier of this loss function
 */
lossType LinGeoLoss::getType() const {
	return LOSS_LIN_GEO;
}
