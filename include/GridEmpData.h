/*
 * GridEmpData.h  
 *  Created on: Sep 29, 2019
 *      
 */

#ifndef SRC_STAT_EMPIRICAL_GRIDEMPDATA_H_
#define SRC_STAT_EMPIRICAL_GRIDEMPDATA_H_

/** This is a subclass of EmpData. The data points in an instance of this class are grid 
 * points. This is a special case of data point because for example in Rappor protocol, 
 * the data points are bit-vectors. This is why we have a constructor that receives 
 * an vector of locations and constructs a GridEmpData from this vector. It is a subclass
 * of EmpData to inherit useful functions in EmpData, e.g. constructing the distribution 
 */

#include "EmpData.h"

class GridEmpData : public EmpData {
public:
	/** default constructor with empty data and mechanism_info (automatically calls the supe constructor) */
	GridEmpData(){}

	/** constructs a  GridEmpData object from an array in which the (int) value 
	 * of each entry is the frequency of the location corresponding to the 
	 * index of that entry. It also receives the dimensions of the grid, and 
	 * associates the data to an identity mechanism (with these dimensions) that 
	 * works on this grid. */ 
	GridEmpData(const std::vector<int>& p_freq_array, const int p_hsize, const int p_vsize);
};



#endif /* SRC_STAT_EMPIRICAL_GRIDEMPDATA_H_ */
