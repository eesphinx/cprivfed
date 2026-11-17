/*
 * UCIAgeBinEmpData.cpp
 *
 *  Created on: Nov 14, 2024
 *      Author: ehab
 *
 * create an (original) EmpData object from a ages dataset file.
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include "AbstractOutputLocalMechanism.h"
#include "CarlosAgeBinEmpData.h"
#include "GridObservation.h"     // we reuse grid observation to represent an age


/** construct the empirical data from a dataset file. Age data are selected between min_age, max_age.
 * Age values are grouped into bins, where age_bin = (age - min_age)/age_bin_size. Secrets are the the
 * set of age bins.
 * duplicationFactor is used to duplicate the frequency of every observation.
 * Needed when the available data is too small. Its default value is 1.
 *
 */
CarlosAgeBinEmpData::CarlosAgeBinEmpData(std::string p_dataset_file, int p_min_age, int p_max_age, int p_age_bin_size, int duplicationFactor)
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


	std::vector<int> age_bin_freq_array(age_bins,0);  // create and initilaize to 0
	std::string line;
	// process the contents of the file
	while (getline(pfile, line)) {
		if (line.length()==0 || line.front()=='#')   continue;
		read_row_count++;
		int age,freq;
		getAgeAndFreq(line, age,freq);
		if (age < min_age || age > max_age) continue;
		int age_bin = (age-min_age)/age_bin_size;
		age_bin_freq_array[age_bin]=duplicationFactor*freq;
		//std::cout << "age:" << age << ", freq: " << freq << "\n";
	}
	std::cout << "Finished after reading " << read_row_count << " lines\n";
	pfile.close();

	// data vector (note zero frequency data are ignored):
	for (int bin=0; bin < age_bins; bin++)
		if (age_bin_freq_array[bin]>0) {
			GridObservation AgeBinObs(bin);
			EmpDatum datum(&AgeBinObs,age_bin_freq_array[bin]);
			data.push_back(datum);   // this uses the copy constructor of EmpDatum datum to resize the vector
		}

	// mechanism info
	mechanism_info = AbstractOutputLocalMechanism::mechTypeToString(IDENT)+" "+std::to_string(1)+" "+std::to_string(age_bins);
}


/**
 * get real age value from the line.
 * The real age is the first item in the line and assuming that line data are comma separated.
 */
void CarlosAgeBinEmpData::getAgeAndFreq(std::string p_line, int& age, int& freq) const {
	std::istringstream iss(p_line);
	char char_token;
	// find the age value
	std::string str_token="";
	while (iss.get(char_token) && char_token != ',' && char_token != ' ')
		str_token = str_token + char_token;
	age = std::stof(str_token);

	// find the frequency
	str_token = "";
	while (iss.get(char_token) && char_token != ',' && char_token != ' ')
		str_token = str_token + char_token;
	freq = std::stoi(str_token);
}



