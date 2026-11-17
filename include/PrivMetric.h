/*
 * PrivMetric.h
 *
 *  Created on: Jun 23, 2019
 *      
 */

#ifndef PRIV_PLAN_PRIVMETRICS_PRIVMETRIC_H_
#define PRIV_PLAN_PRIVMETRICS_PRIVMETRIC_H_

enum PrivacyType {PRIV_GEO, PRIV_MAX, PRIV_MAN};

#include <string>

// This class specifies the necessary function (get_privacy_between) that computes
// the privacy factor based on the coordinates of given cells.
class PrivMetric {
public:
	double eps;  // the privacy per 1 km

	PrivMetric(double p_eps);
	// make a virtual destructor in case we delete a object using a pointer, so the proper derived destructor is called
	virtual ~PrivMetric(){};

	// gets the privacy factor between points i=(i_x, i_y) and j=(j_x,j_y)
    virtual double get_priv_distance(double i_x, double i_y, double j_x, double j_y) const =0;

	virtual std::string get_name() const =0;

	/** we need this function because the Exponential and TC mechanisms store pointers to instances of this class
	 * clone is created in the construction of these mechanisms to avoid the problems when the original objects
	 * are destructed, for instance when the calling function terminates. */
	virtual PrivMetric* clone() const = 0;

	/** get a string representation for a mechanism type */
	static std::string privTypeToString(PrivacyType p_priv_type);

};


#endif /* PRIV_PLAN_PRIVMETRICS_PRIVMETRIC_H_ */
