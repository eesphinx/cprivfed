// the distance (metric) between two points is the max coordinate difference
#include <cmath>
#include <algorithm>
#include "MaxPriv.h"

MaxPriv::MaxPriv(double p_eps)
: PrivMetric(p_eps){
}

double MaxPriv::get_priv_distance(double i_x, double i_y, double j_x, double j_y) const {
	return eps * std::max(abs(i_x-j_x), abs(i_y-j_y));
}

std::string MaxPriv::get_name() const {
	return PrivMetric::privTypeToString(PRIV_MAX);
}

/** returns a pointer to a new copy of the object */
MaxPriv* MaxPriv::clone() const {
	return new MaxPriv(*this);
}
