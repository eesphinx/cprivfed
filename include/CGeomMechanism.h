/*
 * CGeomMechanism.h
 *
 *  Created on: Jan 8, 2026
 *      
 */

#ifndef SRC_MECHANISMS_CENTRAL_CGEOMMECHANISM_H_
#define SRC_MECHANISMS_CENTRAL_CGEOMMECHANISM_H_


#include "CentralMechanism.h"

class CGeomMechanism: public CentralMechanism {
public:
	/**
	 * This mechanism is parameterized by differential privacy parameter dp_eps	 */
	CGeomMechanism(int p_hsize, int p_vsize, double p_dp_eps);

	/* construct the mechanism from the info_string */
    CGeomMechanism(const std::string p_info);
	
    /** Obfuscates the original dataset (histogram) p_origdata by adding discrete Geometric 
    * noise to every bin and returns a noisy histogram */
    DoubleRow obfuscateDataset(const EmpData& p_raw_data) const override;

	double getVariance() const override;

	double get_eps() const;

	std::string get_info() const override;
	std::string getTypeAndPriv(const int p_precision) const override;
private:
    double dp_eps;

};


#endif /* SRC_MECHANISMS_CENTRAL_CGEOMMECHANISM_H_ */
