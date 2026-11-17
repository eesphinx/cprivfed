
#include <iostream>
#include <sstream>
#include <string>
#include <cmath>
#include "PrivMetric.h"
#include "GeoPriv.h"
#include "FullGridSpace.h"
#include "ExponentialMechanism.h"
#include "GSettings.h"


/**
 * This mechanism is generic with respect to privacy metric. Therefore we require to store
 * the privMetric object. Since this class is abstract, we store a pointer to the object.
 */
ExponentialMechanism::ExponentialMechanism(int p_hsize, int p_vsize, double p_cell_side, PrivMetric* p_privObjPtr)
:GridOutputLocalMechanism(p_hsize,p_vsize), cell_side(p_cell_side), mechanism_constants(p_hsize*p_vsize) {
	privObjPtr = p_privObjPtr->clone();
	construct_norm_constants();
}


ExponentialMechanism::ExponentialMechanism(const FullGridSpace& p_space, PrivMetric* p_privObjPtr)
:GridOutputLocalMechanism(p_space.hsize, p_space.vsize), cell_side(p_space.cell_side),
 mechanism_constants(p_space.hsize*p_space.vsize) {
	privObjPtr = p_privObjPtr->clone();
	construct_norm_constants();
}

ExponentialMechanism::ExponentialMechanism(const std::string p_info)
:GridOutputLocalMechanism(0,0), mechanism_constants(0) {   //these variables are set later
	std::string mech_type_str = AbstractOutputLocalMechanism::mechTypeToString(EXP);
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
		mechanism_constants = DoubleColumn(hsize * vsize);
		privObjPtr = new GeoPriv(p_eps);
		construct_norm_constants();
	}
	else {
		std::cerr << "unable to parse the given mechanism info ("
		<< p_info << ") to construct a "<< mech_type_str <<" mechanism\n";
		exit(0);
	}
}


/** release the memory pointed to by the stored pointer */
ExponentialMechanism::~ExponentialMechanism() {
	delete privObjPtr;
	privObjPtr = nullptr;
}

double ExponentialMechanism::getProb(int p_in, int p_out) const {
	int p_in_x = p_in%hsize; int p_in_y = p_in/hsize;
	int p_out_x = p_out%hsize; int p_out_y = p_out/hsize;
	return mechanism_constants.get(p_in)*
			exp(-0.5 * privObjPtr->get_priv_distance(
					p_in_x*cell_side, p_in_y*cell_side,
					p_out_x*cell_side, p_out_y*cell_side));
}

/** sets the normalization constant for every input location (x,y) the mechanism in the linear order y*hsize+x */
void ExponentialMechanism::construct_norm_constants() {
	int in_x, in_y, out_x, out_y;
	for (int in=0; in < hsize*vsize; in++) {  // iterate on inputs
		double sum =0;
		for (int out=0; out< hsize*vsize; out++) {  // iterate on outputs
			in_x = in % hsize; in_y = in / hsize;
			out_x = out % hsize; out_y = out / hsize;
			sum = sum + exp(-0.5 * privObjPtr->get_priv_distance(
					in_x*cell_side, in_y*cell_side,
					out_x*cell_side, out_y*cell_side));
		}
		mechanism_constants.put(in, (double)1.0/sum);
	}
}


/** get the basic info string of the mechanism */
std::string ExponentialMechanism::get_info() const {
	std::string retStr = AbstractOutputLocalMechanism::mechTypeToString(EXP)+ " "+ std::to_string(hsize)+" "+
			std::to_string(vsize) + " "+ std::to_string(cell_side) + " "+
			privObjPtr->get_name() +" " + std::to_string(privObjPtr->eps);
	return retStr;
}

std::string ExponentialMechanism::getTypeAndPriv(const int p_precision) const {
	std::string retStr = AbstractOutputLocalMechanism::mechTypeToString(EXP)+ "_"+
			privObjPtr->get_name() +"_" + GSettings::to_string_with_precision(privObjPtr->eps, p_precision);
	return retStr;
}
