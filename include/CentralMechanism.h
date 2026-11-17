/*
 * LocalMechanism.h
 *
 *  Created on: Aug 26, 2025
 *      
 */

#ifndef SRC_MECHANISMS_CENTRAL_CENTRALMECHANISM_H_
#define SRC_MECHANISMS_CENTRAL_CENTRALMECHANISM_H_

/**
 * This is a base class describing any central mechanism that works on a datasets 
 * of points in full grid input space. Therefore it is a subclass of Mechanism which 
 * has this basic information. This is a useful generalization of CLAP and CGAUSS 
 * mechanisms because they are manipulated in the same way by Estimator_EMP which 
 * combines the underlying histograms. 
 * */

#include <string>
#include "Mechanism.h"
#include "GridEmpData.h"

class DoubleColumn;
class Observation;

class FullGridSpace;
class CentralMechanism : public Mechanism {
public:
	/** basic mechanism attributes */
	//int hsize,vsize;      // The grid for which the mechanism is designed


	// constructor
	CentralMechanism(int p_hsize, int p_vsize)
	   : Mechanism(p_hsize, p_vsize) {};   
	   // constructors are not inherited automatically, we have to forward it explicitly to allow 
	   // subclasses to invoke it (they cannot invoke their grandparent constructors directly)

    /** Computes the histogram of the given raw empirical data (associated with the identity mechanism) 
	 * and returns a noisy histogram after adding random noise to every bin. 
	 * This function must be implemented by every sub mechanism. */
    virtual DoubleRow obfuscateDataset(const EmpData& p_raw_data) const = 0;

    virtual double getVariance() const = 0;

};

#endif /* SRC_MECHANISMS_CENTRAL_CENTRALMECHANISM_H_ */
