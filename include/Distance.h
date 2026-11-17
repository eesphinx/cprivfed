/*
 * Distance.h
 *
 *  Created on: Jun 23, 2019
 *      
 */

#ifndef PRIV_PLAN_DIST_DISTANCE_H_
#define PRIV_PLAN_DIST_DISTANCE_H_

// package priv.plan.dist;

class Distance {
public:

	// This class contains the necessary functions to compute the distances
	// between two locations

	static double get_euclid_dist(double x1, double y1, double x2, double y2);

	static double get_manhat_dist(double i_x, double i_y, double j_x, double j_y);


	// returns the distance in km between to points p1, p2 on earth using their
	// geographic coordinates (lat/lng)
	static double get_distance_on_earth(double p1_lat, double p1_lng, double p2_lat, double p2_lng);

	// adds a south-to-north distance p_dy (in km) and return the new latitude
	static double adjust_latitude (double p_lat, double p_lng, double p_dy);

	// adds a west-to-east distance p_dx (in km) and return the new longitude
	static double adjust_longitude (double p_lat, double p_lng, double p_dx);

private:
	// returns the radial representation of an angle in degrees
	static double get_rad(double x);


};


#endif /* PRIV_PLAN_DIST_DISTANCE_H_ */
