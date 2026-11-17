
#include <iostream>
#include <sstream>
#include <string>
#include <cmath>
#include "DoubleRow.h"
#include "GridObservation.h"
#include "Distributions.h"
#include "CGaussMechanism.h"
#include "GSettings.h"

/*
 * Constructs an optimal CGauss mechanism that satisfies (eps,delta) differential privacy.
 */
CGaussMechanism::CGaussMechanism(int p_hsize, int p_vsize, double p_dp_eps, double p_delta)
:CentralMechanism(p_hsize, p_vsize) {
    dp_eps = p_dp_eps;
    delta = p_delta;

}


/* constructs a CGausse mechanism from the info string: CGAUSS hsize vsize eps delta*/
CGaussMechanism::CGaussMechanism(const std::string p_info)
: CentralMechanism(0,0){
	std::string mechTypeStr = Mechanism::mechTypeToString(CGAUSS);
	std::istringstream iss(p_info);
	std::string p_mechTypeStr;
	if (iss >> p_mechTypeStr >> hsize >> vsize >> dp_eps >> delta) {
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

/** 
* computes the histogram of the given raw dataset (associate with identity mechanism) 
* and returns a new histogram after adding Gaussian noise to every bin. 
*/
DoubleRow CGaussMechanism::obfuscateDataset(const EmpData& p_raw_data) const {
	// get attributes of p_raw_data
	std::istringstream  p_raw_data_iss(p_raw_data.getObfMechInfo());
	std::string p_raw_data_mech_type; p_raw_data_iss >> p_raw_data_mech_type;
	if (p_raw_data_mech_type != Mechanism::mechTypeToString(IDENT)) {
		std::cerr << "unable to obfuscate. Data are already obfuscated by a non-identity mechanism: " << p_raw_data_mech_type << "\n";
		exit(1);
	}

	int p_raw_data_hsize; p_raw_data_iss >> p_raw_data_hsize;
	int p_raw_data_vsize; p_raw_data_iss >> p_raw_data_vsize;

	// validate the obfuscating mechanism with p_raw_data
	if ((p_raw_data_hsize!= hsize)||(p_raw_data_vsize != vsize)) {
		std::cerr << "unable to obfuscate dataset. The grid space of the given " << 
		"data is different from that of the mechanism";
		exit(1);
	}

    // a noisy histogram produced by adding Laplace noise to every bin of the original histogram
    DoubleRow noisy(hsize*vsize);   
	noisy.fillWith(0);
    double gauss_var = get_gaussian_sigma(1.0, dp_eps, delta);
	
    // fill noisy vector with the frequencies of p_raw_data
	for (size_t j=0; j<p_raw_data.getSize(); j++)
		if (GridObservation* obs_ptr = dynamic_cast<GridObservation*>(p_raw_data.data[j].obs_ptr) )
				noisy.put(obs_ptr->location, 
                    (double) noisy.get(obs_ptr->location) + p_raw_data.data[j].frequency );
	
    // add Laplace noise to every entry of noisy_lap
	for (int k=0; k<hsize*vsize; k++)  {
		double d = noisy.get(k)+ Distributions::sampleNormalDouble(0, gauss_var);
		noisy.put(k,d);
	}
	return noisy;
}

double CGaussMechanism::getVariance() const {
	return std::pow(get_gaussian_sigma(1.0, dp_eps, delta),2);
}


double CGaussMechanism::get_eps() const {
	return dp_eps;
}
double CGaussMechanism::get_delta() const {
	return delta;
}

/** get the info string of the mechanism */
std::string CGaussMechanism::get_info() const {
	return Mechanism::mechTypeToString(CGAUSS)+ " "+ std::to_string(hsize)+" "+
			std::to_string(vsize) +" "+std::to_string(dp_eps)+" "+std::to_string(delta);
}

/** get the type and privacy of the mechanism */
std::string CGaussMechanism::getTypeAndPriv(const int p_precision) const {
	return Mechanism::mechTypeToString(CGAUSS)+
	GSettings::to_string_with_precision(dp_eps, p_precision);
}

// returns the value of the Gaussian CDF with mean=0, and sigma=1, i.e. P( N(0,1) <= t )
// for Gaussian N(mu,sigma^2), we have 
// P[N(mu,sigma^2) <=t] = P[ N(0,1) <= (t -mu)/sigma ] = getStdNormal( (t - mu)/sigma )
double getStdNormalCDF(double t) {
	return 0.5 + 0.5*std::erf(t/sqrt(2.0));
}


// The inequality to satisfy:
// Φ(s/(2σ) - (εσ/s)) - e^ε Φ(-s/(2σ) - (εσ/s)) ≤ δ
double inequality_value(double sigma, double s, double epsilon, double delta) {
    if (sigma <= 0) return std::numeric_limits<double>::infinity();
    
    double term1 = s / (2.0 * sigma) - (epsilon * sigma / s);
    double term2 = -s / (2.0 * sigma) - (epsilon * sigma / s);
    
    double phi1 = getStdNormalCDF(term1);
    double phi2 = getStdNormalCDF(term2);
    
    return phi1 - std::exp(epsilon) * phi2 - delta;
}

// Find smallest σ using binary search (guaranteed convergence)
double find_minimal_sigma_binary(double s, double epsilon, double delta, 
                                double sigma_min, double sigma_max,
                                double tol = 1e-12,
                                int max_iter = 1000) {
    
    // Check if σ_max already satisfies the inequality
    double f_max = inequality_value(sigma_max, s, epsilon, delta);
    if (f_max <= 0) {
        // σ_max satisfies, need to search lower
        // Keep reducing σ_max until it doesn't satisfy
        while (sigma_max > sigma_min && f_max <= 0) {
            sigma_max /= 2.0;
            f_max = inequality_value(sigma_max, s, epsilon, delta);
        }
        sigma_max *= 2.0;  // Last one that satisfied
        f_max = inequality_value(sigma_max, s, epsilon, delta);
    }
    
    // Check σ_min
    double f_min = inequality_value(sigma_min, s, epsilon, delta);
    
    // If σ_min already satisfies, it's the answer
    if (f_min <= 0) {
        return sigma_min;
    }
    
    // Ensure we bracket the solution: f_min > 0, f_max < 0
    // If not, expand sigma_max
    while (f_min * f_max > 0 && sigma_max < 1e100) {
        sigma_max *= 2.0;
        f_max = inequality_value(sigma_max, s, epsilon, delta);
    }
    
    if (f_min * f_max > 0) {
        std::cerr << "Error: Could not bracket solution. f_min=" << f_min 
                  << ", f_max=" << f_max << std::endl;
        return std::numeric_limits<double>::quiet_NaN();
    }
    
    // Binary search
    for (int iter = 0; iter < max_iter; iter++) {
        double sigma_mid = (sigma_min + sigma_max) / 2.0;
        double f_mid = inequality_value(sigma_mid, s, epsilon, delta);
        
        if (std::abs(f_mid) < tol || (sigma_max - sigma_min) < tol) {
            return sigma_mid;
        }
        
        if (f_mid > 0) {
            // Inequality not satisfied (value > 0), need larger σ
            sigma_min = sigma_mid;
            f_min = f_mid;
        } else {
            // Inequality satisfied (value ≤ 0), can try smaller σ
            sigma_max = sigma_mid;
            f_max = f_mid;
        }
    }
    
    return (sigma_min + sigma_max) / 2.0;
}

// returns sigma (square root of the variance) of gaussian noise needed 
// to satisfy eps-delta differential privacy using (Theorem A.1 in 
// Dwark&Roth book)
double CGaussMechanism::get_gaussian_sigma(double p_sensitivity, double p_eps, double p_delta) {

	double min_error = 1e-7; // error tolerance of delta

	double sigma_l = min_error;
    // initial upper bound for sigma (obtained from Theorem A.1 in Dwark&Roth book)
    // Note that Sensitivity for a single histogram bin is D_1 = D_2. 
    // for eps<1 any sigma > p_sensitivity* sqrt(2*log(1.25/p_delta)) / p_eps  
    // satisfies (eps, delta), hence the above value is upper bound for sigma
    // for eps>1 sigma must be smaller than its value with eps = 1, which is 
    // p_sensitivity* sqrt(2*log(1.25/p_delta)) / 1. Therefore we set the upper 
    // bound as follows
    double sigma_h;
    if (p_eps < 1)  
	    sigma_h = p_sensitivity* sqrt(2*log(1.25/p_delta)) / p_eps;
    else
        sigma_h = p_sensitivity* sqrt(2*log(1.25/p_delta));
    // variance = sigma^2. 
    // Return this value if we don't want to optimize the Gaussian noise
    //	return sigma_h;  
	// std::cout << "sigma_h = " << sigma_h << "\n";

    // Trying to optimize the variance by binary search
	
    // u is the smallest delta that a given sigma satisfies
    // evaluated using Eq (6) in (Theorem 8 in [Balle,Wang 2018])
    // As we decrease sigma, u increases. So we try to find the 
    // smallest sigma such that u does not exceed delta, i.e. u = delta. 
    
    return find_minimal_sigma_binary(p_sensitivity,p_eps,p_delta,sigma_l,sigma_h);
	
	
	// double u; 
	// double sigma;
	// do {
	// 	sigma = (sigma_l + sigma_h)/2;
	// 	// compute the minimal delta satisfied for the given variance
	// 	u = getStdNormalCDF(p_sensitivity/(2*sigma) - (p_eps*sigma)/p_sensitivity ) -
	// 			exp(p_eps)*getStdNormalCDF((-1)*p_sensitivity/(2*sigma)- (p_eps*sigma)/p_sensitivity);
	// 	// if u is bigger than p_delta then we have to increase the variance
	// 	// otherwise we can make sigma smaller.
	// 	if (u > p_delta) sigma_l = sigma; else sigma_h = sigma;
	// }
	// while (abs(p_delta - u) > min_error);
	// std::cout << "sigma = " << sigma << "\n";
	// return sigma;
}




