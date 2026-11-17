/*
 * UCILabTestsEmpData.cpp
 *
 *  Created on: Jul 6, 2026
 *      Author: ehab
 *
 * create an (original) EmpData object from a given dataset file.
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include "Mechanism.h"
#include "UCILabTestsEmpData.h"
#include "GridObservation.h"     // we reuse grid observation to represent a relative age


/** construct the empirical data from a dataset file. LabTests data are selected between min_age, max_age.
 * labTests values are grouped into bins, where labTests_bin = (labTests - min_labTests)/labTests_bin_size. Secrets are the the
 * set of labTests bins.
 */
UCILabTestsEmpData::UCILabTestsEmpData(std::string p_dataset_file, int p_min_labTests, int p_max_labTests, int p_labTests_bin_size)
:min_labTests(p_min_labTests), max_labTests(p_max_labTests), labTests_bin_size(p_labTests_bin_size){
	int read_row_count=0;     // number of valid rows read from in_data_file
	std::ifstream pfile(p_dataset_file);
	if (!pfile.is_open()) {
		std::cerr << "Unable to open dataset file " << p_dataset_file;
		exit(0);
	}

	// In the following (max_labTests-min_labTests+1) is the number of different labTests values between min,max
	// Casted to Double, and divided by labTests_bin_size to get the number of bins (with fractions).
	// We take the "ceil" because if the number of bins has fractions, this fraction should be considered a final bin.
	int labTests_bins = (int) ceil( (double)(max_labTests-min_labTests+1)/labTests_bin_size );

	std::vector<int> labTests_bin_freq_array(labTests_bins, 0); // create array and initialize to 0
	std::string line;
	// process the contents of the file
	while (getline(pfile, line)) {
		if (line.length()==0 || line.front()=='#')   continue;
		read_row_count++;
		int labTests = getLabTests(line);
		if (labTests < min_labTests || labTests > max_labTests) continue;
		int labTests_bin = (labTests-min_labTests)/labTests_bin_size;
		labTests_bin_freq_array[labTests_bin]++;
	}
	pfile.close();

	// data vector (note zero frequency data are ignored):
	for (int bin=0; bin < labTests_bins; bin++)
		if (labTests_bin_freq_array[bin]>0) {
			GridObservation labTestsBinObs(bin);
			EmpDatum datum(&labTestsBinObs,labTests_bin_freq_array[bin]);
			data.push_back(datum);   // this uses the copy constructor of EmpDatum datum to resize the vector
		}

	// mechanism info
	mechanism_info = Mechanism::mechTypeToString(IDENT)+" "+std::to_string(1)+" "+std::to_string(labTests_bins);
}


/**
 * get real labTests value from the line.
 * The real labTests is the 13th item in the line and assuming that the fields in the line are comma separated.
 */
int UCILabTestsEmpData::getLabTests(std::string p_line) const {
	
	std::stringstream ss(p_line);
    std::string field;
    // Read the first 12 fields
    for (int i = 0; i < 12; ++i)
    {
        if (!std::getline(ss, field, ','))
            throw std::runtime_error("Line has fewer than 12 fields.");
    }
    // Read the 13th field
    if (!std::getline(ss, field, ','))
        throw std::runtime_error("Line has fewer than 13 fields.");
	return std::stoi(field);
}



