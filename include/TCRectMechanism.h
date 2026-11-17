/*
 * TCRectMechanism.h
 *
 *  Created on: Feb 7, 2020
 *      
 */

#ifndef SRC_MECHANISMS_PLAN_TCRECTMECHANISM_H_
#define SRC_MECHANISMS_PLAN_TCRECTMECHANISM_H_

/* Tight-constraints mechanism on a rectangular grid (using automorphism symmetries) */

#include "GridOutputLocalMechanism.h"
#include "DoubleColumn.h"

/* Internal class to describe the continuous coordinates of a point relative to the grid center, unit distance is cell_side */
class Point {
public:
	double x; double y;
	Point(double p_x, double p_y) {
		x=p_x; y=p_y;
	}
};


class PrivMetric;

class TCRectMechanism: public GridOutputLocalMechanism {
public:

	/** this constructor creates the TC mechanism specified by the hsize, vsize of the grid,
	 * the cell_side (in km), and a pointer to privacy parameter */
	TCRectMechanism(int p_hsize, int p_vsize, double p_cell_side, const PrivMetric* p_priv);

	TCRectMechanism(const std::string p_info);

	/** to delete the pointer member priv */
	~TCRectMechanism();

	std::string get_info() const override;
	std::string getTypeAndPriv(const int p_precision) const override;

private:

	/* constructs all data structures required to evaluate the probability for in and out cells */
	void initialize();

	/* get the element of sym_phi for two classes p_i, p_j */
	double get_sym_phi_entry(int p_i, int p_j) const;

	/* the coordinates of in and out points are in the original grid (relative to the corner) */
	double getProb(int p_in, int p_out) const override;

	/*  get the index of the class that contains a given point p described by its continuous coordinates relative to the center of the grid, e.g. (0.5,0.5) */
	int get_sym_class(Point p) const;

	/* get the point that represents the class having the given index */
	Point get_sym_class_rep(int p_class_index) const;

	/* get the number of classes in the grid */
	int get_sym_classes_count() const;

	/* a pointer to a a privacy metric (created or cloned in the constructor) used to evaluate the mechanism probabilities */
	PrivMetric* priv;

	/* important for this mechanism because the privacy metric, hence the probabilities depends on it */
	double cell_side;

	/* mu_sym[c] is mechanism diagonal element for a location in symmetry class c */
	DoubleColumn mu_sym;

	/* the coordinates (relative to the center of the grid) of the point representing class 0 in the upper-right quarter. */
	double x0, y0;

	/* the number of cells in the +ve x-direction, +ve y-directions starting from the center of the grid */
	int qh, qv;

};



#endif /* SRC_MECHANISMS_PLAN_TCRECTMECHANISM_H_ */
