/*
 * CarlosAgeBinEmpData.h
 *
 *  Created on: Nov 15, 2024
 *      
 */

#ifndef SRC_DATASETS_ADULTDSET_CARLOSAGEBINEMPDATA_H_
#define SRC_DATASETS_ADULTDSET_CARLOSAGEBINEMPDATA_H_

#include "EmpData.h"

// Relative ages empirical data. Age data are relative to min_age. That is relative age = real age - min_age

class CarlosAgeBinEmpData : public EmpData {
public:

	/** construct the empirical data from Carlos's dataset file. Each line contains Age, frequency.
	 * Age data are selected between min_age, max_age.
	 * Age values are grouped into bins, where age_bin = (age - min_age)/age_bin_size. Secrets are the the
	 * set of age bins.
	 * duplicationFactor is used to duplicate the frequency of every observation. Needed
	 * when the available data is too small. Its default value is 1.
	 */
	CarlosAgeBinEmpData(std::string p_dataset_file, int p_min_age, int p_max_age, int p_age_bin_size, int duplicationFactor=1);

private:
	/**
	 * get real age and its frequency from the line.
	 * The real age and frequency are the first and second items in the line and assuming that line data are comma separated.
	 */
	void getAgeAndFreq(std::string p_line, int& age, int& freq) const;
	int min_age;
	int max_age;
	int age_bin_size;   // allows to have a range of ages in one bin (slot). Useful for experiments on Rappor
};




#endif /* SRC_DATASETS_ADULTDSET_CARLOSAGEBINEMPDATA_H_ */
