// This class specifies the necessary function (get_privacy_between) that computes 
// the privacy factor based on the coordinates of given cells.
#include "Distance.h"
#include "GeoPriv.h"

GeoPriv::GeoPriv(double p_eps)
: PrivMetric(p_eps) {
}

// gets the privacy factor between points i=(i_x, i_y) and j=(j_x,j_y)
double GeoPriv::get_priv_distance(double i_x, double i_y, double j_x, double j_y) const {
	return eps*Distance::get_euclid_dist(i_x, i_y, j_x, j_y);   // accessing a static member of Distance
}


std::string GeoPriv::get_name() const{
	return PrivMetric::privTypeToString(PRIV_GEO);
}

/** returns a pointer to a new copy of the object */
GeoPriv* GeoPriv::clone() const {
	return new GeoPriv(*this);
}

