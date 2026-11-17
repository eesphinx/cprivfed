/*
 * GridEmpData.cpp
 *
 *  Created on: Sep 29, 2019
 *      Author: ehab
 */

#include "GridEmpData.h"
#include "GridObservation.h"
#include "Mechanism.h"

	/** constructs a  GridEmpData object from an array in which the (int) value 
	 * of each entry is the frequency of the location corresponding to the 
	 * index of that entry. It also receives the dimensions of the grid, and 
	 * associates the data to an identity mechanism (with these dimensions) that 
	 * works on this grid. The main purpose is to convert an array of frequencies 
	 * into a GridEmpData object */
	GridEmpData::GridEmpData(const std::vector<int>& p_freq_array, const int p_hsize, const int p_vsize) {
		mechanism_info = Mechanism::mechTypeToString(IDENT)
		+" "+std::to_string(p_hsize)+" "+std::to_string(p_vsize);
		// note: zero frequency data are ignored
		for (size_t loc=0;loc<p_freq_array.size();loc++)
			if (p_freq_array[loc]>0) {
				GridObservation loc_obs(loc);
				EmpDatum datum(&loc_obs,p_freq_array[loc]);
				data.push_back(datum);   // this uses the copy constructor of EmpDatum datum to resize the vector
			}
	}


