/**
 * This class describes a full grid overlaid on a rectangular geographical region. This region
 * is specified by a south boundary (p_south), i.e. minimum latitude and west boundary (p_west),
 * i.e. the minimum longitude. This class provides functions to obtain the losses with respect
 * to GDsets (i.e. check-ins) for a given GridOutputMechanism. The grid coordinates of the input points
 * in every utility function are determined by the region borders, and the region cell_side
 * (through get_cell functions). The grid coordinates of the output points are determined by the
 * grid info of the given mechanism (which may have a different granularity). Note that priors
 * for utility are defined on the input grid.
 * convention: x coordinate increases as we go towards the east direction and y coordinate
 * increases as we go towards the north.
 * */


#include <iostream>
#include "Distance.h"
#include "GeoFullGridSpace.h"



/** The constructor specifies the grid dimensions, cell side, and the geographical coorrdinates */
GeoFullGridSpace::GeoFullGridSpace(FullGridSpace& p_space, double p_south, double p_west)
: FullGridSpace(p_space), south(p_south), west(p_west)  {
}

/** The constructor specifies the grid dimensions, cell side, and the geographical coordinates */
GeoFullGridSpace::GeoFullGridSpace(int p_hsize, int p_vsize, double p_cell_side, double p_south, double p_west)
:FullGridSpace(p_hsize, p_vsize, p_cell_side), south(p_south), west(p_west) {
}

/**
 * the following two functions returns xy discrete coordinates for a point specified as (lat, lng) relative to
 * a the longitude (west) and latitude (south). recall x and y increase in the east and north directions respectively
 */
int GeoFullGridSpace::get_cell_x_of_lng(double p_lng) const {
	// get the cell of the given point
	double xdist = Distance::get_distance_on_earth(south, west, south, p_lng);
	int in_x =  (int) (xdist/cell_side);
	if (in_x >= hsize) {
		std::cerr << "warning: input point is slightly out of grid: xdist = " << xdist << "\n";
		in_x--;
	}
	return in_x;
}

/** get the y coordinate of the cell enclosing a point with latitude p_lat in a full grid relative to (south, west) */
int GeoFullGridSpace::get_cell_y_of_lat(double p_lat) const {
	// get the cell of the given point on the grid specified by p_cell_side
	double ydist = Distance::get_distance_on_earth(south, west, p_lat, west);
	int in_y =  (int) (ydist/cell_side);
	if (in_y >= vsize) {
		std::cerr << "warning: input point is slightly out of grid: ydist = " << ydist << "\n";
		in_y--;
	}
	return in_y;
}

