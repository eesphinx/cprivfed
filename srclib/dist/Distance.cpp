// package priv.plan.dist;



// This class contains the necessary functions to compute the distances
// between two locations

#include <cmath>
#include "Distance.h"

double Distance::get_euclid_dist(double x1, double y1, double x2, double y2) {
	return sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2));
}

double Distance::get_manhat_dist(double i_x, double i_y, double j_x, double j_y){
	return abs(i_x-j_x)+ abs(i_y-j_y);
}
// returns the radial representation of an angle in degrees
double Distance::get_rad(double x) {
	return (x * M_PI/180);
}

// returns the distance in km between to points p1, p2 on earth using their
// geographic coordinates (lat/lng)
double Distance::get_distance_on_earth(double p1_lat, double p1_lng, double p2_lat, double p2_lng) {
	double R = 6378.137;              // Earth mean radius in km
	double dLat = get_rad(p2_lat - p1_lat);
	double dLong = get_rad(p2_lng - p1_lng);
	double a = sin(dLat / 2) * sin(dLat / 2) +
			cos(get_rad(p1_lat)) * cos(get_rad(p2_lat)) *
			sin(dLong / 2) * sin(dLong / 2);
	double c = 2 * atan2(sqrt(a), sqrt(1 - a));
	double d = R * c;
	return d; // returns the distance in km
}

// adds a south-to-north distance p_dy (in km) and return the new latitude
double Distance::adjust_latitude (double p_lat, double p_lng, double p_dy) {
	double R = 6378.137;              // Earth mean radius in km
	return p_lat + (p_dy/R)*(180/M_PI);
}

// adds a west-to-east distance p_dx (in km) and return the new longitude
double Distance::adjust_longitude (double p_lat, double p_lng, double p_dx) {
	double R = 6378.137;              // Earth mean radius in km
	return p_lng + (p_dx/R)*(180/M_PI)/cos(p_lat* M_PI/180);
}

