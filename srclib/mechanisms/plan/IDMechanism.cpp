/*
 * IDMechanism.cpp
 *
 *  Created on: Sep 28, 2019
 *      Author: ehab
 */

#include <iostream>
#include <sstream>
#include "IDMechanism.h"

/* This mechanism is defined only by the space dimensions	 */
IDMechanism::IDMechanism(int p_hsize, int p_vsize)
:GridOutputLocalMechanism(p_hsize,p_vsize){
}

std::string IDMechanism::get_info() const {
	return AbstractOutputLocalMechanism::mechTypeToString(IDENT)+ " "+ std::to_string(hsize)+" "+
			std::to_string(vsize);
}

std::string IDMechanism::getTypeAndPriv(const int p_precision) const {
	return AbstractOutputLocalMechanism::mechTypeToString(IDENT);
}


IDMechanism::IDMechanism(const std::string p_info)
:GridOutputLocalMechanism(0,0) {
	std::string mech_type_str = AbstractOutputLocalMechanism::mechTypeToString(IDENT);
	std::istringstream iss(p_info);
	std::string p_mech_type_str;
	if (iss >> p_mech_type_str >> hsize >> vsize) {
		if (p_mech_type_str != mech_type_str) {
			std::cerr << "unable to construct " << mech_type_str << " mechanism from info string:\n "
					<< "info string must start with " << mech_type_str;
			exit(0);
		}
	}
	else {
		std::cerr << "unable to parse the mechanism info to construct "<< mech_type_str <<" mechanism";
	}
}

double IDMechanism::getProb(int p_in, int p_out) const {
	if (p_out == p_in)
		return 1;
	else
		return 0;
}



