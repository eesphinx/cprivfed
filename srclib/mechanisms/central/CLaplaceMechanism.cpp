
#include <iostream>
#include <sstream>
#include <string>
#include <cmath>
#include "DoubleRow.h"
#include "GridObservation.h"
#include "Distributions.h"
#include "CLaplaceMechanism.h"
#include "GSettings.h"

/*
 * Constructs a CLaplace mechanism that satisfies eps differential privacy.
 */
CLaplaceMechanism::CLaplaceMechanism(int p_hsize, int p_vsize, double p_dp_eps)
:CentralMechanism(p_hsize, p_vsize) {
    dp_eps = p_dp_eps;
}


/* constructs a CLaplace mechanism from the info string: CLAP hsize vsize eps*/
CLaplaceMechanism::CLaplaceMechanism(const std::string p_info)
: CentralMechanism(0,0){
	std::string mechTypeStr = Mechanism::mechTypeToString(CLAP);
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

/** Obfuscates the original dataset (histogram) p_origdata by adding continous Laplace 
 * noise to every bin and returns a noisy histogram 
 */
DoubleRow CLaplaceMechanism::obfuscateDataset(const EmpData& p_raw_data) const {
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
	double lap_scale = 1.0/dp_eps; // laplace scale parameter (b) = sensitivity/eps_dp;
	
    // fill noisy_lap vector with the original frequencies of p_raw_data
	for (size_t j=0; j<p_raw_data.getSize(); j++)
		if (GridObservation* obs_ptr = dynamic_cast<GridObservation*>(p_raw_data.data[j].obs_ptr) )
				noisy.put(obs_ptr->location, 
                    (double) noisy.get(obs_ptr->location) + p_raw_data.data[j].frequency );
	
    // add Laplace noise to every entry of noisy_lap
	for (int k=0; k<hsize*vsize; k++)  {
		double d = noisy.get(k)+ Distributions::sampleLaplaceDouble(0, lap_scale);
		noisy.put(k,d);
	}
	return noisy;
}

double CLaplaceMechanism::getVariance() const {
	return 2.0/(dp_eps*dp_eps);
}


double CLaplaceMechanism::get_eps() const {
	return dp_eps;
}

/** get the basic info string of the mechanism */
std::string CLaplaceMechanism::get_info() const {
	return Mechanism::mechTypeToString(CLAP)+ " "+ std::to_string(hsize)+" "+
			std::to_string(vsize) +" "+std::to_string(dp_eps);
}

/** get the type and privacy of the mechanism */
std::string CLaplaceMechanism::getTypeAndPriv(const int p_precision) const {
	return Mechanism::mechTypeToString(CLAP)+
	GSettings::to_string_with_precision(dp_eps, p_precision);
}
