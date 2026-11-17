/*
 * UCIAgeEmpData.cpp
 *
 *  Created on: Aug 29, 2020
 *      Author: ehab
 *
 * create an (original) EmpData object from a given dataset file.
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include "Mechanism.h"
#include "UCIAgeEmpData.h"
#include "GridObservation.h"     // we reuse grid observation to represent a relative age


/** construct the empirical data from a dataset file. Age data are selected between min_age, max_age.
 * Age values are grouped into bins, where age_bin = (age - min_age)/age_bin_size. Secrets are the the
 * set of age bins.
 */
UCIAgeEmpData::UCIAgeEmpData(std::string p_dataset_file, int p_min_age, int p_max_age, int p_age_bin_size)
:min_age(p_min_age), max_age(p_max_age), age_bin_size(p_age_bin_size){
	int read_row_count=0;     // number of valid rows read from in_data_file
	std::ifstream pfile(p_dataset_file);
	if (!pfile.is_open()) {
		std::cerr << "Unable to open dataset file " << p_dataset_file;
		exit(0);
	}

	// In the following (max_age-min_age+1) is the number of different ages between min,max
	// Casted to Double, and divided by age_bin_size to get the number of bins (with fractions).
	// We take the "ceil" because if the number of bins has fractions, this fraction should be considered a final bin.
	int age_bins = (int) ceil( (double)(max_age-min_age+1)/age_bin_size );

	std::vector<int> age_bin_freq_array(age_bins, 0); // create array and initialize to 0
	std::string line;
	// process the contents of the file
	while (getline(pfile, line)) {
		if (line.length()==0 || line.front()=='#')   continue;
		read_row_count++;
		int age = getAge(line);
		if (age < min_age || age > max_age) continue;
		int age_bin = (age-min_age)/age_bin_size;
		age_bin_freq_array[age_bin]++;
	}
	pfile.close();

	// data vector (note zero frequency data are ignored):
	for (int bin=0; bin < age_bins; bin++)
		if (age_bin_freq_array[bin]>0) {
			GridObservation AgeBinObs(bin);
			EmpDatum datum(&AgeBinObs,age_bin_freq_array[bin]);
			data.push_back(datum);   // this uses the copy constructor of EmpDatum datum to resize the vector
		}

	// mechanism info
	mechanism_info = Mechanism::mechTypeToString(IDENT)+" "+std::to_string(1)+" "+std::to_string(age_bins);
}


/**
 * get real age value from the line.
 * The real age is the first item in the line and assuming that line data are comma separated.
 */
int UCIAgeEmpData::getAge(std::string p_line) const {
	std::istringstream iss(p_line);
	std::string str_token="";
	char char_token;
	while (iss.get(char_token) && char_token != ',' && char_token != ' ')
		str_token = str_token + char_token;
	return std::stoi(str_token);
}



