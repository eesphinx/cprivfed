/*
 * TCRectMechanism.cpp
 *
 *  Created on: Feb 7, 2020
 *      Author: ehab
 */



/** Tight-constraints mechanism on a rectangular grid (using automorphism) **/

#include <cmath>
#include <iostream>
#include <sstream>
#include "GeoPriv.h"
#include "DoubleColumn.h"
#include "TCRectMechanism.h"
#include "PrivMetric.h"
#include "GSettings.h"


/** this constructor creates the TC mechanism specified by the size of the grid,
 * the cell_side (in km), and a pointer to privacy parameter */
TCRectMechanism::TCRectMechanism(int p_hsize, int p_vsize, double p_cell_side, const PrivMetric* p_priv)
:GridOutputLocalMechanism(p_hsize, p_vsize), cell_side(p_cell_side){
	priv = p_priv->clone();
	initialize();
}

TCRectMechanism::TCRectMechanism(const std::string p_info)
:GridOutputLocalMechanism(0,0) {   //these variables are set later
	std::string mech_type_str = AbstractOutputLocalMechanism::mechTypeToString(TC);
	std::istringstream iss(p_info);
	std::string p_mech_type_str; std::string p_priv_type_str; double p_eps;

	if (iss >> p_mech_type_str >> hsize >> vsize >> cell_side >> p_priv_type_str >> p_eps) {
		if (p_mech_type_str != mech_type_str) {
			std::cerr << "unable to construct " << mech_type_str << " mechanism from info string:\n "
					<< "info string must start with " << mech_type_str;
			exit(0);
		}
		if (p_priv_type_str != PrivMetric::privTypeToString(PRIV_GEO)) {
			std::cerr << "unable to construct " << mech_type_str << " mechanism from info string:\n "
					<< "specified privacy ("<< p_priv_type_str <<") metric is not supported";
			exit(0);
		}

		// construction work
		priv = new GeoPriv(p_eps);
		initialize();
	}
	else {
		std::cerr << "unable to parse the mechanism info to construct "<< mech_type_str <<" mechanism";
		exit(0);
	}

}


TCRectMechanism::~TCRectMechanism() {
	delete priv;
	priv = nullptr;
}


/**  get the index of the class that contains the given point. The cells of the upper-right
 * quarter describe the individual different classes. (p_x,p_y) are continuous coordinates
 * relative to the center of the grid. For example p = (0.5,0.5) */
int TCRectMechanism::get_sym_class(Point p) const {
	// the class is identified by going through the cells of the upper-right quarter of the grid in zigzag manner horizontally first.
	return (int) ( (abs(p.y) - y0)*qh + (abs(p.x) - x0) );
}


/* get the point (in continuous coordinates relative to the center of the grid)
 * that represents the class having the given index */
Point TCRectMechanism::get_sym_class_rep(int p_class_index) const {
	// the point is located by going through the cells of the upper-right quarter of the grid in zigzag manner horizontally first.
	double d_retx = x0 + p_class_index%qh;
	double d_rety = y0 + p_class_index/qh;
	return Point(d_retx,d_rety);
}



// get the number of classes in the grid (size of mu)
int TCRectMechanism::get_sym_classes_count() const {
	return qh*qv;
}


/* the coordinates of in and out points are in the original grid in the linear order */
double TCRectMechanism::getProb(int p_in, int p_out) const {

	int in_x=p_in%hsize; int in_y = p_in/hsize;
	int out_x=p_out%hsize; int out_y = p_out/hsize;

	double ret;
	double mid_x = (hsize-1)*0.5;
	double mid_y = (vsize-1)*0.5;

	// get the coordinates of the input and output points relative to the center.
	double rel_in_x = in_x - mid_x;  double rel_in_y = in_y - mid_y;
	double rel_out_x = out_x - mid_x;  double rel_out_y = out_y - mid_y;

	int class_ind = get_sym_class(Point(rel_out_x,rel_out_y));

	// required probability is the product of the entry in mu that corresponds
	// to the output point and the factor between output and input.
	ret = mu_sym.get(class_ind)*
			exp(-1*priv->get_priv_distance(rel_in_x*cell_side, rel_in_y*cell_side,
					rel_out_x*cell_side, rel_out_y*cell_side));
	return ret;
}

void TCRectMechanism::initialize() {

	// initialize data related to symmetric classes
	if (hsize%2==0) { // if hsize is even
		x0 = 0.5; qh = hsize/2;
	}
	else {
		x0 =0; qh = (hsize-1)/2 + 1;
	}

	if (vsize%2==0) { // if vsize is even
		y0 = 0.5; qv = vsize/2;
	}
	else {
		y0 =0; qv = (vsize-1)/2 + 1;
	}

	//construct the matrix phi_sym;
	int n = get_sym_classes_count();
	DoubleMatrix phi_sym(n,n);
	for (int i=0;i<n;i++)
		for (int j=0;j<n;j++)
			phi_sym.put(i,j, get_sym_phi_entry(i,j));
	//std::cout << "phi = " << phi_sym.toString() << "\n";


	//construct the column vector mu_sym;
	DoubleColumn b(n);
	b.fillWith(1.0);
	mu_sym = phi_sym.solve(b); // get mu_sym that satisfies phi_sym * mu_sym = 1

	// check the resulting mu. If it is nonzero, the mechanism exists, otherwise it does not.
	if (mu_sym.min()<0) {
		std::cerr << "TC-mechanism does not exist. vector entry is -ve.";
		std::cerr << "mu (transposed)= " << mu_sym.transpose().toString();
		exit(0);
	}

}

// get the entry of sym_phi for two classes p_i, p_j
double TCRectMechanism::get_sym_phi_entry(int p_i, int p_j) const {
	double sum=0;
	Point i_rep_point = get_sym_class_rep(p_i);
	Point j_rep_point = get_sym_class_rep(p_j);

	if (j_rep_point.x==0 && j_rep_point.y==0)     // j_rep_x = j_rep_y = 0
		sum = sum + exp(-1*priv->get_priv_distance(
				i_rep_point.x*cell_side, i_rep_point.y*cell_side,
				0, 0));
	else if (j_rep_point.y==0) {  // j_rep_x>0 , j_rep_y=0
		sum = sum + exp(-1*priv->get_priv_distance(
				i_rep_point.x*cell_side, i_rep_point.y*cell_side,
				j_rep_point.x*cell_side, 0));
		sum = sum + exp(-1*priv->get_priv_distance(
				i_rep_point.x*cell_side, i_rep_point.y*cell_side,
				-j_rep_point.x*cell_side, 0));
	}
	else if (j_rep_point.x==0) {  // j_rep_x=0 , j_rep_y>0
		sum = sum + exp(-1*priv->get_priv_distance(
				i_rep_point.x*cell_side, i_rep_point.y*cell_side,
				0, j_rep_point.y*cell_side));
		sum = sum + exp(-1*priv->get_priv_distance(
				i_rep_point.x*cell_side, i_rep_point.y*cell_side,
				0, -j_rep_point.y*cell_side));
	}
	else {        // j_rep_x>0 , j_rep_y >0
		sum = sum + exp(-1*priv->get_priv_distance(
				i_rep_point.x*cell_side, i_rep_point.y*cell_side,
				j_rep_point.x*cell_side, j_rep_point.y*cell_side));
		sum = sum + exp(-1*priv->get_priv_distance(
				i_rep_point.x*cell_side, i_rep_point.y*cell_side,
				j_rep_point.x*cell_side, -j_rep_point.y*cell_side));
		sum = sum + exp(-1*priv->get_priv_distance(
				i_rep_point.x*cell_side, i_rep_point.y*cell_side,
				-j_rep_point.x*cell_side, j_rep_point.y*cell_side));
		sum = sum + exp(-1*priv->get_priv_distance(
				i_rep_point.x*cell_side, i_rep_point.y*cell_side,
				-j_rep_point.x*cell_side, -j_rep_point.y*cell_side));
	}
	return sum;
}

/** get the basic info string of the mechanism */
std::string TCRectMechanism::get_info() const {
	std::string retStr = Mechanism::mechTypeToString(TC)+ " "+ std::to_string(hsize)+" "+
			std::to_string(vsize) + " "+ std::to_string(cell_side) + " "+
			priv->get_name() +" " + std::to_string(priv->eps);
	return retStr;
}

std::string TCRectMechanism::getTypeAndPriv(const int p_precision) const {
	std::string retStr = Mechanism::mechTypeToString(TC)+"_" 
	+ priv->get_name() +"_" + GSettings::to_string_with_precision(priv->eps, p_precision);
	return retStr;
}








