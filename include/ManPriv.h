/*
 * ManPriv.h
 *
 *  Created on: Jun 23, 2019
 *      
 */

#ifndef PRIV_PLAN_PRIVMETRICS_MANPRIV_H_
#define PRIV_PLAN_PRIVMETRICS_MANPRIV_H_

#include "PrivMetric.h"
// This class specifies the necessary function (get_privacy_between) that computes
// the privacy factor based on the coordinates of given cells.
class ManPriv : public PrivMetric {

public:

	ManPriv(double p_eps);

	double get_priv_distance(double i_x, double i_y, double j_x, double j_y) const;

	std::string get_name() const;

	/** returns a pointer to a new copy of the object */
	ManPriv* clone() const;  // return type is co-variant to the return type (PrivMetric) of the super class.
};



#endif /* PRIV_PLAN_PRIVMETRICS_MANPRIV_H_ */
