/*
 * Min_Loss_Given_dp_Mechanism.h
 *
 *  Created on: Feb 6, 2020
 *      
 */

#ifndef SRC_MECHANISMS_PLAN_MIN_LOSS_GIVEN_DP_MECHANISM_H_
#define SRC_MECHANISMS_PLAN_MIN_LOSS_GIVEN_DP_MECHANISM_H_

#include "OptimizedMechanism.h"

class Min_Loss_Given_dp_Mechanism: public OptimizedMechanism {
public:

	/**
	 * constructs a dp-mechanisms that is optimal for a given privacy metric, a loss metric and a prior.
	 * The space is given to determine the distance between cells to compute the privacy and loss distances.
	 */
	Min_Loss_Given_dp_Mechanism(const Prior& p_prior, const FullGridSpace& p_xspace,
			const PrivMetric* p_privMetric, const Loss* p_lossMetric,
			const std::string p_mechFile);
};


#endif /* SRC_MECHANISMS_PLAN_MIN_LOSS_GIVEN_DP_MECHANISM_H_ */
