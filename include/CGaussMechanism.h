/*
 * CGaussMechanism.h
 *
 *  Created on: Aug 28, 2025
 *      
 */

#ifndef SRC_MECHANISMS_CENTRAL_CGAUSSMECHANISM_H_
#define SRC_MECHANISMS_CENTRAL_CGAUSSMECHANISM_H_


#include "CentralMechanism.h"

class CGaussMechanism: public CentralMechanism {
public:
	/**
	 * This mechanism is parameterized by differential privacy parameter dp_eps and delta	 
     * Constructs the optimal mechanism (least variance) that satisfies (eps, delta)-DP */
	CGaussMechanism(int p_hsize, int p_vsize, double p_dp_eps, double p_delta);

	/* construct the mechanism from the info_string */
    CGaussMechanism(const std::string p_info);

	/** 
     * computes the histogram of the given raw dataset (associate with identity mechanism) 
     * and returns a new histogram after adding Gaussian noise to every bin. 
    */
    DoubleRow obfuscateDataset(const EmpData& p_raw_data) const override;

    double getVariance() const override; 

	double get_eps() const;
	double get_delta() const;

	std::string get_info() const override;
	std::string getTypeAndPriv(const int p_precision) const override;
private:
    double dp_eps;
    double delta;

    // returns sigma (square root of the variance) of gaussian noise needed 
    // to satisfy eps-delta differential privacy using (Theorem A.1 in 
    // Dwark&Roth book) or optimized using (Theorem 8 in [Balle,Wang 2018]).
    static double get_gaussian_sigma(double p_sensitivity, double p_eps, double p_delta);

};


#endif /* SRC_MECHANISMS_CENTRAL_CGAUSSMECHANISM_H_ */
