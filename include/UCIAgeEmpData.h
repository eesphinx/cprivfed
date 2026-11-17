/*
 * UCIAgeEmpData.h
 *
 *  Created on: Aug 29, 2020
 *      
 */

#ifndef SRC_DATASETS_ADULTDSET_UCIAGEEMPDATA_H_
#define SRC_DATASETS_ADULTDSET_UCIAGEEMPDATA_H_

#include "EmpData.h"

// Relative ages empirical data. Age data are relative to min_age. That is relative age = real age - min_age
// Data are organized in bins using the parameter bin_size given to construct the object

class UCIAgeEmpData : public EmpData {
public:

	/** construct the empirical data from a dataset file. Age data are selected between min_age, max_age.
	 * Age values are grouped into bins, where age_bin = (age - min_age)/age_bin_size. Secrets are the the
	 * set of age bins.
	 */
	UCIAgeEmpData(std::string p_dataset_file, int p_min_age, int p_max_age, int p_age_bin_size);

private:
	/**
	 * get real age value from the line.
	 * The real age is the first item in the line and assuming that line data are comma separated.
	 */
	int getAge(std::string p_line) const;
	int min_age;
	int max_age;
	int age_bin_size;   // allows to have a range of ages in one bin (slot). Useful for experiments on Rappor
};




#endif /* SRC_DATASETS_ADULTDSET_UCIAGEEMPDATA_H_ */
