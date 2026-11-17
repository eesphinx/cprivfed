//package priv.plan.loss;


// loss is 0 within distance saturation_dist, and is max_loss for larger distances
// note: in the special case if saturation_dist=0, the loss is 0 only when dist=0, 
// and max_loss otherwise. 

# include <cmath>
# include "BinLoss.h"

BinLoss::BinLoss(double p_saturation_dist, double p_max_loss)
:Loss(p_saturation_dist, p_max_loss) {
}

double BinLoss::get_loss(double i_x, double i_y, double j_x, double j_y) const {
	return get_loss_by_dist(sqrt(pow(i_x-j_x, 2)+ pow(i_y-j_y, 2)));
}

double BinLoss::get_loss_by_dist(double p_dist) const {
	if (p_dist<=saturation_dist)
		return 0;
	else
		return max_loss;
}

lossType BinLoss::get_name() const {
	return LOSS_BIN;
}
