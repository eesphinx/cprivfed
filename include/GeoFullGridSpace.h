/*
 * GeoFullGridSpace.h
 *
 *  Created on: Sep 22, 2019
 *      
 */

#ifndef SPACES_GEOFULLGRIDSPACE_H_
#define SPACES_GEOFULLGRIDSPACE_H_

#include "FullGridSpace.h"

class GeoFullGridSpace : public FullGridSpace {

public:

	double south;
	double west;


	/** The constructor from a full grid space and the geographical coordinates */
	GeoFullGridSpace(FullGridSpace& p_space, double p_south, double p_west);

	/** The constructor specifies the grid dimensions, cell side, and the geographical coordinates */
	GeoFullGridSpace(int p_hsize, int p_vsize, double p_cell_side, double p_south, double p_west);

	/**
	 * the following two functions returns xy discrete coordinates for a point specified as (lat, lng) relative to
	 * a the longitude (west) and latitude (south). recall x and y increase in the east and north directions respectively
	 */
	int get_cell_x_of_lng(double p_lng) const;

	/** get the y coordinate of the cell enclosing a point with latitude p_lat in a full grid relative to (south, west) */
	int get_cell_y_of_lat(double p_lat) const;


};


#endif /* SRC_SPACES_GEOFULLGRIDSPACE_H_ */
