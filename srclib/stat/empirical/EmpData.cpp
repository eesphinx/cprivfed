/*
 * EmpData.cpp
 *
 *  Created on: Sep 23, 2019
 *      Author: ehab
 */
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include "GSettings.h"
#include "Distributions.h"
#include "DoubleRow.h"
#include "GridObservation.h"
#include "BitVectorObservation.h"
#include "AbstractOutputLocalMechanism.h"
#include "GridEmpData.h"
#include "EmpData.h"


/** constructs an EmpData by obfuscating the given raw  
 * empirical data using a given AbstractOutputLocalMechanism. 
 * Note: Only local mechanism are allowed to generate EmpData */
EmpData::EmpData(const EmpData& p_raw_empdata, const AbstractOutputLocalMechanism* p_local_mech_ptr) {
	mechanism_info = p_local_mech_ptr -> get_info();  // arbitrary local mechanism
	obfuscateAndAdd(p_raw_empdata);  // This checks if the data are raw (associated with the identity mechanism IDENT)

}

/** create an Empdata object from an empirical data file. */
EmpData::EmpData(std::string& p_data_file)
:mechanism_info(""){
	std::ifstream pfile(p_data_file);
	if (!pfile.is_open()) {
		std::cerr << "unable to open empirical data file: " << p_data_file;
		exit(0);
	}
	// process the contents of the file
	std::string info_tag = "#info: ";
	int info_tag_len = info_tag.length();
	bool mechanism_info_found = false;

	std::string line;
	//find the mechanism info
	while (!mechanism_info_found && getline(pfile,line)) {
		if (line.substr(0, info_tag_len) == info_tag) {
			mechanism_info_found = true;
			mechanism_info = line.substr(info_tag_len);
			//std::istringstream iss(mechanism_info);
			//iss >> mech_type_str;
		}
	}

	if (mechanism_info == "") {
		std::cerr << "invalid EmpData: mechanism information is not specified";
		exit(0);
	}

	// std::unique_ptr allows us to forget about deleting the pointer.
	// use it if not needed to pass the pointer somewhere outside the scope. 
	std::unique_ptr<Mechanism> mech_ptr(Mechanism::createMechanism(mechanism_info));
	auto local_mech_ptr = dynamic_cast<AbstractOutputLocalMechanism*>(mech_ptr.get());
	if (!local_mech_ptr) {
    	std::cerr << "Mechanism is not an AbstractOutputLocalMechanism \n";
		exit(1);
	}
	while (getline(pfile,line)) {
		if (line.length()==0 || line.front()=='#')
			continue;
		EmpDatum datum(line, local_mech_ptr);     // construct a generic datum from line (given the mechanism)
		data.push_back(datum);
	}
	pfile.close();
}


/** write the data to the given empdata file */
void EmpData::write_data_to_file(const std::string& p_file) const{
	std::ofstream pfile (p_file);
	if (pfile.is_open()) {
		// write data
		pfile << "#info: " << mechanism_info << "\n";
		for (size_t k=0; k< getSize() ; k++)
			pfile << data[k].toString() << "\n";
		pfile.close();
	}
	else {
		std::cerr << "Unable to open file for writing the empirical data";
		exit(0);
	}
}

/** add the data of p_empdata (must be associated with the same mechanism) */
void EmpData::add(const EmpData& p_empdata) {
	// validate to have the same mechanism
	if ( p_empdata.getObfMechInfo() != mechanism_info) {
		std::cerr << "unable to add empirical data. The mechanism or space are not the same: ";
		exit(0);
	}
	// add the data of p_empdata 
	for (size_t l=0; l<p_empdata.getSize();l++) {
		size_t i=0; while (i< data.size() && *data[i].obs_ptr != *p_empdata.data[l].obs_ptr) i++;
		if (i < getSize()) {  // obsevation seen before?
			data[i].frequency = data[i].frequency + p_empdata.data[l].frequency;
		}
		else {
			EmpDatum obfDatum(p_empdata.data[l].obs_ptr ,1); // construct a new datum with the new observation and frequency 1
			data.push_back(obfDatum);
		}
	}
}

/** Obfuscate (original) p_empdata using existing mech_info and add the noisy data */ 
void EmpData::obfuscateAndAdd(const EmpData& p_raw_empdata) {

	std::unique_ptr<Mechanism> mech_ptr(Mechanism::createMechanism(mechanism_info));
	// Mechanism* mech_ptr = Mechanism::createMechanism(mechanism_info);

	// get attributes of p_empdata. The main purpose is to get hsize, vsize
	std::istringstream  p_raw_empdata_iss(p_raw_empdata.getObfMechInfo());
	std::string p_raw_empdata_mech_type; p_raw_empdata_iss >> p_raw_empdata_mech_type;
	if (p_raw_empdata_mech_type != Mechanism::mechTypeToString(IDENT)) {
		std::cerr << "The given empirical data are not raw. They have been produced by: " 
			<< p_raw_empdata_mech_type << "\n";
		exit(1);
	}

	int p_grid_empdata_hsize; p_raw_empdata_iss >> p_grid_empdata_hsize;
	int p_grid_empdata_vsize; p_raw_empdata_iss >> p_grid_empdata_vsize;

	// validate the obfuscating mechanism with p_empdata
	if ((p_grid_empdata_hsize!= mech_ptr->hsize)||(p_grid_empdata_vsize != mech_ptr->vsize)) {
		std::cerr << "unable to obfuscate and add data. The grid space of the given " << 
		"data (" << p_grid_empdata_hsize << "," << p_grid_empdata_vsize << 
		") is different from that of the mechanism that generated the existing data (" <<
		mech_ptr->hsize << "," << mech_ptr->vsize << ")";
		// delete(mech_ptr); mech_ptr = nullptr; Object deleted automatically since we use unique_ptr().
		exit(0);
	}
	// cast the mechanism of the current data into an AbstractOutputLocalMechanism to allow obfuscation
	AbstractOutputLocalMechanism* local_mech_ptr = 
			dynamic_cast<AbstractOutputLocalMechanism*>(mech_ptr.get()); //mech_ptr.get() returns the needed raw pointer
	if (!local_mech_ptr) {
		std::cerr << "unable to cast the mechanism of exiting data into a local mechanismto abfuscate the new data";
		// delete(mech_ptr); mech_ptr = nullptr;  // not needed since we use unique_ptr
		exit(0);
	}
	// Now everything should be fine
	for (size_t l=0; l<p_raw_empdata.getSize();l++)  // obfuscate and add
		if (GridObservation* original_obs_ptr = dynamic_cast<GridObservation*>(p_raw_empdata.data[l].obs_ptr)) {
			for (int u=0; u<p_raw_empdata.data[l].frequency; u++) {
				// a new observation is created and a pointer is returned
				Observation* new_obs_ptr = local_mech_ptr->obfuscateLoc(original_obs_ptr->location); 
				// if the new observation exists in the data, then update frequency, 
				// otherwise add a new datum for it with frequency 1
				size_t i=0; while (i< data.size() && *data[i].obs_ptr != *new_obs_ptr) i++;
				if (i < getSize()) {  // obsevation seen before?
					data[i].frequency++;
				}
				else {
					// construct a new datum with the new observation and frequency 1
					// Note that the cinstructor of EmpDatum clones new_obs_ptr
					EmpDatum obfDatum(new_obs_ptr ,1); 
					data.push_back(obfDatum);
				}
				delete(new_obs_ptr); new_obs_ptr = nullptr; // important
			}
		}
		else {
			std::cerr << "unable to cast the raw data into grid data";
			// delete(mech_ptr); mech_ptr = nullptr; // Not  needed by the use of unique_ptr
			exit(1);
		}
}


/** split the data into multiple empirical data objects according to the given probability distribution */
std::vector<EmpData> EmpData::randSplit(const DoubleRow& p_split_distribution) const {
	// check that p_split_distribution is a valid distribution
	int prec = std::stoi(GSettings::get_setting_value(precision));
	double error = pow(10, -1*prec);
	if (abs(p_split_distribution.sum() - 1) > error ) {
		std::cerr << "unable to split the grid data into multiple data objects: invalid split distribution:\n"
				<< "sum of probabilities = " << p_split_distribution.sum();
		exit(0);
	}
	std::vector<EmpData> emp_data_objects(p_split_distribution.getLength());
	for (int i=0; i< p_split_distribution.getLength(); i++)
		emp_data_objects[i].mechanism_info = mechanism_info;  // grid_emp_data_objects

	for (size_t l=0; l<getSize();l++)
		for (int u=0; u<data[l].frequency; u++) {
			int index = Distributions::sampleDiscrete(p_split_distribution,1)[0];
			if (emp_data_objects[index].getSize()>0) {
				int last = emp_data_objects[index].getSize()-1; // index of the last inserted item
				if (*emp_data_objects[index].data[last].obs_ptr == *data[l].obs_ptr) {
					emp_data_objects[index].data[last].frequency++;
					continue;
				}
			}
			EmpDatum datum(data[l].obs_ptr ,1); // construct a new datum the current observation and frequency 1
			emp_data_objects[index].data.push_back(datum); // add datum
		}

	return emp_data_objects;
}


/**
 * returns the empirical distribution on observables as a row vector. The j-th element corresponds
 * to the j-th datum (observation , frequency) in data. P.S. Since it is constructed from data, it
 * contains no zero entries (i.e. excludes observations that do not appear in the data).
 * It is used in IBU since it requires no zeros in the empirical distribution.
 */
DoubleRow EmpData::getDistributionOnObservables() const {
	DoubleRow emp_pr_z(getSize());
	for (size_t j=0; j< getSize(); j++)
		emp_pr_z.put(j, data[j].frequency);
	// Normalize
	double sum = emp_pr_z.sum();
	emp_pr_z.divi(sum);
	return emp_pr_z;
}

/** returns the number of observed outcomes with nonzero frequencies **/
size_t EmpData::getSize() const {
	return data.size();
}

/** return the sum of all frequencies in data */
int EmpData::get_totalObservations() const {
	int ret = 0;
	for (size_t i=0; i< getSize(); i++) {
		ret = ret + data[i].frequency;
	}
	return ret;
}

/** returns the mechanism info string */
std::string EmpData::getObfMechInfo() const {
	return mechanism_info;
}

/**  specific to BitVector observations data */

/**
 * This returns a row vector with length equal to the number of bits in one bit vector. The
 * j-th element is Tj/n where Tj is the number of bitVectors that have bit j equal to 1, and
 * n is the total number of observed bit vectors in the data (totalObservations).
 * This function is used for Rappor estimators BOTRappor_INVN and BOTRappor_INVP.  */
DoubleRow EmpData::get_relative_frequency_of_bits() const {
	std::istringstream iss(mechanism_info);
	std::string mechType;
	iss >> mechType;
	if (mechType != Mechanism::mechTypeToString(BOTRAPPOR)) {
		std::cerr << "unable to get empirical bit probabilities from data";
		exit(0);
	}
	int hsize,vsize;
	iss >> hsize >> vsize;

	DoubleRow retRow(hsize*vsize); retRow.fillWith(0);
	for (size_t i=0; i< getSize(); i++)
		if (BitVectorObservation* obs_ptr = dynamic_cast<BitVectorObservation*>(data[i].obs_ptr)) {
			for (int j=0; j<hsize*vsize; j++)
				if (obs_ptr->bitVector[j])
					retRow.put(j, retRow.get(j)+ data[i].frequency);
		}
		else {
			std::cerr << "unable to get empirical bit probabilities from data.\n"
					<< "unable to cast the empirical data observations to bitVectors";
		}

	int totalObservedBitVectors = get_totalObservations();
	retRow.divi(totalObservedBitVectors);
	return retRow;
}

/**
 * Takes a list of EmpData objects, combines the elements that are associated 
 * with the same mechanism, and returns the resulting list
 */
std::vector<EmpData> EmpData::compress(const std::vector<EmpData>& p_list_empdata ) {
	std::vector<EmpData> new_list_empdata; // construct a new list
	for (size_t l=0; l < p_list_empdata.size(); l++) {
		// search for a matching empdata element in the new_list_empdata
		size_t i=0; while (i< new_list_empdata.size() && new_list_empdata[i].getObfMechInfo() != p_list_empdata[l].getObfMechInfo()) i++;
		if (i < new_list_empdata.size()) {  // a similar empdata found in the new list?
			new_list_empdata[i].add(p_list_empdata[l]);  // add the empdata object to the existing one
		}
		else {
			EmpData e(p_list_empdata[l]); // make a copy 
			new_list_empdata.push_back(e); // create a new element in the new list with the new copy.
		}
	}
	return new_list_empdata;
}




