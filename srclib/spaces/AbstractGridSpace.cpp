#include "AbstractGridSpace.h"
/**
 * Describes an abstract space X of secrets, as a general subset of a grid of locations
 * The grid is specified by hsize,vsize, cell_side. The subset may be the entire grid as
 * in (FullGridSecretsSpace) or a discrete set of cells, or
 * a convex hull of locations. Every case of these is implemented by a
 * subclass of SecretsSpace.
 **/

// default constructor
AbstractGridSpace::AbstractGridSpace(int p_hsize, int p_vsize, double p_cell_side)
:hsize(p_hsize), vsize(p_vsize), cell_side(p_cell_side) {
}


/** return the continuous (double) x coordinate for a given location in the linear order (linear_index)
 * This is used for computing the loss between two points, hence computing the EMD distance */
double AbstractGridSpace::getx(int linear_index) const {
	return (linear_index % hsize)*cell_side;
}

/** return the continuous y coordinate for a given location in the linear order (linear_index) */
double AbstractGridSpace::gety(int linear_index) const {
	return (linear_index/hsize)*cell_side;
}
