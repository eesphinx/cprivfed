/*
 * UCILabTestsEmpData.h
 *
 *  Created on: July 6, 2026
 *      
 */

#ifndef SRC_DATASETS_DIABETES_UCILABTESTSEMPDATA_H_
#define SRC_DATASETS_DIABETES_UCILABTESTSEMPDATA_H_

#include "EmpData.h"

// lab tests performed for every patient during their stay in a diabetes hospital. 
// https://archive.ics.uci.edu/dataset/296/diabetes+130-us+hospitals+for+years+1999-2008
// The values are relative in the senste that a patient value = labTests_ - min_labTests
// where the min and max values are given to the construct the object.


class UCILabTestsEmpData : public EmpData {
public:

	/** construct the empirical data from a dataset file. LabTests data are selected between min_labTests, max_labTests.
	 * labTests values are grouped into bins, where labTests_bin = (labTests - min_labTests)/labTests_bin_size. Secrets are the
	 * set of labTests bins.
	 */
	UCILabTestsEmpData(std::string p_dataset_file, int p_min_labTests, int p_max_labTests, int p_labTests_bin_size);

private:
	/**
	 * get real labTests value from the line.
	 * The real labTests is the 13th item in the line and assuming that line data are comma separated.
	 */
	int getLabTests(std::string p_line) const;
	int min_labTests;
	int max_labTests;
	int labTests_bin_size;   // allows to have a range of labTests in one bin (slot). Useful for experiments on Rappor
};




#endif /* SSRC_DATASETS_DIABETES_UCILABTESTSEMPDATA_H_ */
