// This class specifies the necessary function (get_privacy_between) that computes 
// the privacy factor based on the coordinates of given cells.
#include "ManPriv.h"
#include "Distance.h"

ManPriv::ManPriv(double p_eps)
:PrivMetric(p_eps){
}

double ManPriv::get_priv_distance(double i_x, double i_y, double j_x, double j_y) const {
	return eps*Distance::get_manhat_dist(i_x, i_y, j_x, j_y);
}

std::string ManPriv::get_name() const {
	return PrivMetric::privTypeToString(PRIV_MAN);
}

/** returns a pointer to a new copy of the object */
ManPriv* ManPriv::clone() const {
	return new ManPriv(*this);
}
