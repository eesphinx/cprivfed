
#include <iostream>
#include <sstream>
#include <string>
#include <cmath>
#include "DoubleRow.h"
#include "GridObservation.h"
#include "Distributions.h"
#include "CGeomMechanism.h"
#include "GSettings.h"

/*
 * Constructs a CGeometric mechanism that satisfies eps differential privacy.
 */
CGeomMechanism::CGeomMechanism(int p_hsize, int p_vsize, double p_dp_eps)
:CentralMechanism(p_hsize, p_vsize) {
    dp_eps = p_dp_eps;
}


/* constructs a CGeometric mechanism from the info string: CGEOM hsize vsize eps*/
CGeomMechanism::CGeomMechanism(const std::string p_info)
: CentralMechanism(0,0){
	std::string mechTypeStr = Mechanism::mechTypeToString(CGEOM);
	std::istringstream iss(p_info);
	std::string p_mechTypeStr;
	if (iss >> p_mechTypeStr >> hsize >> vsize >> dp_eps) {
		if (p_mechTypeStr != mechTypeStr) {
			std::cerr << "unable to construct " << mechTypeStr << " mechanism from info string:\n "
					<< "info string must start with " << mechTypeStr;
			exit(1);
		}
	}
	else {
		std::cerr << "unable to parse the mechanism info to construct "<< mechTypeStr <<" mechanism";
        exit(1);
	}
}

/** Obfuscates the original dataset histogram p_origdata by adding discrete geometric 
 * noise to every bin and returns the resulting noisy histogram 
 */
DoubleRow CGeomMechanism::obfuscateDataset(const EmpData& p_raw_data) const {
	// get attributes of p_origData
	std::istringstream  p_origdata_iss(p_raw_data.getObfMechInfo());
	std::string p_origdata_mech_type; p_origdata_iss >> p_origdata_mech_type;
	if (p_origdata_mech_type != Mechanism::mechTypeToString(IDENT)) {
		std::cerr << "unable to obfuscate. Data are already obfuscated by a non-identity mechanism: " << p_origdata_mech_type << "\n";
		exit(1);
	}

	int p_origdata_hsize; p_origdata_iss >> p_origdata_hsize;
	int p_origdata_vsize; p_origdata_iss >> p_origdata_vsize;

	// validate the obfuscating mechanism with p_empdata
	if ((p_origdata_hsize!= hsize)||(p_origdata_vsize != vsize)) {
		std::cerr << "unable to obfuscate dataset. The grid space of the given " << 
		"data is different from that of the mechanism";
		exit(1);
	}
    
    // a noisy histogram produced by adding Laplace noise to every bin of the histogram of the raw data
    DoubleRow noisy(hsize*vsize);   
	noisy.fillWith(0);
	
    // fill the noisy vector with the original frequencies of p_raw_data
	for (size_t j=0; j<p_raw_data.getSize(); j++)
		if (GridObservation* obs_ptr = dynamic_cast<GridObservation*>(p_raw_data.data[j].obs_ptr) )
				noisy.put(obs_ptr->location, 
                    (double) noisy.get(obs_ptr->location) + p_raw_data.data[j].frequency );

	for (int k=0; k<hsize*vsize; k++)  {
		int d = noisy.get(k)+ Distributions::sampleSymGeometricInt(std::exp(-dp_eps));
		noisy.put(k,d);
	}
	return noisy;
}

double CGeomMechanism::getVariance() const {
	return 2.0 * std::exp(-dp_eps)/std::pow(1-std::exp(-dp_eps),2);
}


double CGeomMechanism::get_eps() const {
	return dp_eps;
}

/** get the basic info string of the mechanism */
std::string CGeomMechanism::get_info() const {
	return Mechanism::mechTypeToString(CGEOM)+ " "+ std::to_string(hsize)+" "+
			std::to_string(vsize) +" "+std::to_string(dp_eps);
}

/** get the type and privacy of the mechanism */
std::string CGeomMechanism::getTypeAndPriv(const int p_precision) const {
	return Mechanism::mechTypeToString(CGEOM)+
	GSettings::to_string_with_precision(dp_eps, p_precision);
}
