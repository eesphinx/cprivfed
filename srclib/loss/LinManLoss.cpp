//package priv.plan.loss;

#include <cmath>
#include "LinManLoss.h"

LinManLoss::LinManLoss(double p_saturation_dist, double p_max_loss)
:Loss(p_saturation_dist,p_max_loss) {
	}

double LinManLoss::get_loss(double i_x, double i_y, double j_x, double j_y) const {
	          return get_loss_by_dist(abs(i_x-j_x)+  abs(i_y-j_y));
		}

double LinManLoss::get_loss_by_dist(double p_man_dist) const {
		if (p_man_dist<=saturation_dist)
			return p_man_dist;
		else
			return max_loss;
	    }


lossType LinManLoss::getType() const {
		return LOSS_LIN_MAN;
	}
