/*
 * Mechanism.cpp
 *
 *  Created on: August 25, 2025
 *      Author: ehab
 */

#include <TTaxicabMechanism.h>
#include <iostream>
#include <sstream>
#include "Mechanism.h"
#include "IDMechanism.h"
#include "ExponentialMechanism.h"
#include "TCRectMechanism.h"
#include "PRRespMechanism.h"
#include "TPGeomMechanism.h"
#include "TPLapMechanism.h"
#include "TVLGeomMechanism.h"
#include "TVLLapMechanism.h"
#include "RMechanism.h"
#include "ScaledMechanism.h"
#include "BasicOneTimeRappor.h"
#include "PMatrixMechanism.h"
#include "MaxRiskGivenMaxLossMechanism.h"
#include "CLaplaceMechanism.h"
#include "CGeomMechanism.h"
#include "CGaussMechanism.h"
#include "FullGridSpace.h"

/* constructor */
Mechanism::Mechanism(int p_hsize, int p_vsize)
: hsize(p_hsize), vsize(p_vsize){
}

/* constructs a mechanism with a Full grid input space. However the mechanism output space is generic. */
Mechanism::Mechanism(const FullGridSpace& p_space)
: hsize(p_space.hsize) , vsize(p_space.vsize) {

}

/** returns the mechanism type as a string*/
std::string Mechanism::getTypeStr() const {
	std::istringstream iss(get_info());
	std::string mech_type_str;
	if (iss >> mech_type_str)
		return mech_type_str;
	else {
		std::cerr << "unable to return the mechanism type";
		exit(0);
	}
}


/** get a string representation for a mechanism type */
std::string Mechanism::mechTypeToString(MechanismType p_mech_type) {
	switch(p_mech_type) {
	case IDENT : return "IDENT";
	case VLLAP : return "VLLAP";
	case PLAP : return "PLAP";
	case VLGEOM : return "VLGEOM";
	case PGEOM : return "PGEOM";
	case TAXICAB : return "TAXICAB";
	case EXP : return "EXP";
	case TC : return "TC";
	case KRR : return "KRR";
	case PMAT : return "PMAT";
	case RMECH : return "RMECH";
	case SMECH : return "SMECH";
	case BOTRAPPOR : return "BOTRAPPOR";
	case MAXADVRISK : return "MAXADVRISK";  // shokri's mechanism
	case CLAP : return "CLAP";  // central Laplace mechanism
	case CGEOM : return "CGEOM";  // central geometric mechanism
	case CGAUSS : return "CGAUSS";  // central Gauss mechanism
	default :
		std::cerr << "invalid mechanism type \n";
		exit(0);
	}
//	return "";

}

/** create a mechanism from an info string, and return a pointer.
 * Important: delete the object in the calling procedure after use */
Mechanism* Mechanism::createMechanism(std::string p_info) {
	std::istringstream iss(p_info);
	std::string p_mech_type_str;
	iss >> p_mech_type_str;

	if (p_mech_type_str == Mechanism::mechTypeToString(IDENT))
		return new IDMechanism(p_info);
	if (p_mech_type_str == Mechanism::mechTypeToString(EXP))
		return new ExponentialMechanism(p_info);
	if (p_mech_type_str == Mechanism::mechTypeToString(TC))
		return new TCRectMechanism(p_info);
	if (p_mech_type_str == Mechanism::mechTypeToString(KRR))
		return new PRRespMechanism(p_info);
	if (p_mech_type_str == Mechanism::mechTypeToString(VLGEOM))
		return new TVLGeomMechanism(p_info);
	if (p_mech_type_str == Mechanism::mechTypeToString(PGEOM))
		return new TPGeomMechanism(p_info);
	if (p_mech_type_str == Mechanism::mechTypeToString(TAXICAB))
		return new TTaxicabMechanism(p_info);
	if (p_mech_type_str == Mechanism::mechTypeToString(VLLAP))
		return new TVLLapMechanism(p_info);
	if (p_mech_type_str == Mechanism::mechTypeToString(PLAP))
		return new TPLapMechanism(p_info);
	if (p_mech_type_str == Mechanism::mechTypeToString(BOTRAPPOR))
		return new BasicOneTimeRappor(p_info);
	if (p_mech_type_str == Mechanism::mechTypeToString(MAXADVRISK))
		return new MaxRiskGivenMaxLossMechanism(p_info);
	if (p_mech_type_str == Mechanism::mechTypeToString(PMAT))
		return new PMatrixMechanism(p_info);
	if (p_mech_type_str == Mechanism::mechTypeToString(RMECH))
		return new RMechanism(p_info);
	if (p_mech_type_str == Mechanism::mechTypeToString(SMECH))
		return new ScaledMechanism(p_info);
	if (p_mech_type_str == Mechanism::mechTypeToString(CLAP))
		return new CLaplaceMechanism(p_info);
	if (p_mech_type_str == Mechanism::mechTypeToString(CGEOM))
		return new CGeomMechanism(p_info);
	if (p_mech_type_str == Mechanism::mechTypeToString(CGAUSS))
		return new CGaussMechanism(p_info);
	std::cerr << "unable to construct obfuscating mechanism from info:" << p_info <<"\n";
	exit(0);
}


