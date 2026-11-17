/*
 * MaxPriv.h
 *
 *  Created on: Jun 23, 2019
 *      
 */

#ifndef PRIV_PLAN_PRIVMETRICS_MAXPRIV_H_
#define PRIV_PLAN_PRIVMETRICS_MAXPRIV_H_

#include "PrivMetric.h"
// the distance (metric) between two points is the max coordinate difference
class MaxPriv: public PrivMetric {
public:
	MaxPriv(double p_eps);

	double get_priv_distance(double i_x, double i_y, double j_x, double j_y) const;

	std::string get_name() const;

	/** returns a pointer to a new copy of the object */
	MaxPriv* clone() const;  // return type is co-variant to the return type (PrivMetric) of the super class.

};




#endif /* PRIV_PLAN_PRIVMETRICS_MAXPRIV_H_ */
