
#include <string>
#include <sstream>    // for reading from string as a file
#include <fstream>    // for reading and writing to files
#include <iostream>   // for cerr, cout
#include <iomanip>
#include <ios>
#include <vector>
#include "DoubleRow.h"
#include "GSettings.h"
#include "Prior.h"



/** 
 * The prior describes a probability distribution over a discrete points
 *  in the space. These points are the centers of the grid that contains them
 * (as cell centers).
 */ 

/** constructs a prior on hsize*vsize grid using a given row vector. This is useful to construct
 * a prior based on a distribution, which is represented by DoubleRow vector */
Prior::Prior(const DoubleRow& p_vect, int p_hsize, int p_vsize)
:hsize(p_hsize), vsize(p_vsize), prior_vect(p_hsize*p_vsize), prior_info("") {
	prior_vect = p_vect;
}


/** For creating a prior from a file. It should include the matrix */
Prior::Prior(std::string p_file)
:prior_vect(1), prior_info("") {  // prior_vector must be explicitly initialized here since it has no default constructor.
	std::vector<double> prob_list;
	std::string line;         // to hold a line string from the input file

	std::ifstream pfile(p_file); // construct an ifstream associated with p_file. This also opens the file (see refence on ifstream)
	if (pfile.is_open())
	{
		int curr_y=0; int curr_x;
		double prob;
		while (getline(pfile,line)) {
			if (line.length() == 0 || line.front()=='#')     // skip empty or commented line
				continue;
			// if the above conditions are not true for the current line
			curr_x =0;
			std::stringstream ss(line);
			while (ss >> prob) {
				prob_list.push_back(prob);
				curr_x++;
			}
			if (curr_x >0) {
				// already found elements with coordinate curr_y
				hsize=curr_x;     // set the number of different x coordinates
				curr_y=curr_y+1;   // increase the curr_y if we have we have
			}
		}
		vsize=curr_y;            // set the number of rows
		prior_vect.resize(prob_list.size());
		for (size_t k=0;k<prob_list.size();k++)
			prior_vect.put(k, prob_list.at(k));

		pfile.close();
	}
	else {
		std::cerr << "Unable to open prior file " << p_file << std::endl;
		exit(0);
	}
}

/**
 * Creates a prior object from a file having the probabilities in the linear order
 * i.e., the first probability in p_str corresponds to P(y=0,x=0)
 */
Prior::Prior(std::string p_file, int p_hsize, int p_vsize)
:hsize(p_hsize), vsize(p_vsize), prior_vect(p_hsize*p_vsize), prior_info("") {
	std::string line;         // to hold a line string from the input file
	double prob;
	std::ifstream pfile(p_file);
	if (pfile.is_open()) {
		if (getline(pfile,line)) {
			std::stringstream ss(line);
			int i=0;
			while (ss >> prob && i<prior_vect.getLength()) {
				prior_vect.put(i,prob);
				i++;
			}
			if (i<prior_vect.getLength()) {
				std::cerr << "The data in the file is insufficient for the specified dimensions";
				exit(0);
			}
		}
		pfile.close();
	}
	else {
		std::cerr << "Unable to open prior file " << p_file;
		exit(0);
	}
}

/** return a uniform prior on the grid */
Prior Prior::uniform(int p_hsize, int p_vsize) {
	double u = 1.0/(p_hsize*p_vsize);
	return Prior(  DoubleRow(p_hsize*p_vsize).fillWith(u), p_hsize, p_vsize );
}


/** return a random prior on the grid */
Prior Prior::random(int p_hsize, int p_vsize) {
	DoubleRow vec(p_hsize*p_vsize);
	vec.fillRand();
	double sum = vec.sum();
	return Prior(vec.divi(sum), p_hsize, p_vsize  );
}

///** fill the prior matrix with random probabilities */
//void Prior::fill_prior_random() {
//	prior_vect.fillRand();
//	double sum = prior_vect.sum();
//	prior_vect.divi(sum); // normalize
//}


/** write the prior_matrix to the given file. If p_with_info=true, the header information is written too. */
void Prior::write_prior_as_matrix(std::string p_file, bool p_with_info) const {
	std::ofstream pfile (p_file); // construct an ofstream associated with p_file. This also opens the file (see refence on ofstream)
	if (pfile.is_open()) {
		if (p_with_info)
			pfile << prior_info << "\n\n";
		int index =0;
		int default_precision = std::stoi(GSettings::get_setting_value(precision)); // or use std::atoi(GSettings::get_setting_value(precision).c_str())
		pfile << std::fixed << std::setprecision(default_precision);
		for (int y=0; y<vsize; y++) {
			for (int x=0; x< hsize; x++) {
				pfile << prior_vect.get(index) << " ";
				index++;
			}
			pfile << "\n";
		}
		//pfile.flush();
		pfile.close();
	}
	else {
		std::cerr << "Unable to open file";
		exit(0);
	}
}


/** returns the prior data in a single string format as required by Kostas first row is written first, etc. */
std::string Prior::toString () const {
	std::ostringstream oss;
	for (int i = 0; i< hsize*vsize; i++)
		oss << prior_vect.get(i) << " ";
	std::string s = oss.str();
    return s;
}

	// write the prior as one vector in S format
void Prior::write_prior_as_vector(std::string p_file) const {
	std::ofstream pfile (p_file);
	if (pfile.is_open()) {
		for (int i = 0; i< hsize*vsize; i++)
			pfile << prior_vect.get(i) << " ";
		pfile.flush();
		pfile.close();
	}
	else std::cout << "Unable to open file" << std::endl;
}



/** get the probability of cell p specified in the linear order (row first) */
double Prior::getProb(int p) const {
	return prior_vect.get(p);
}

/**  returns the row vector representation of the prior */
DoubleRow Prior::get_prior_vector() const {
	return prior_vect;
}


void Prior::add_info(std::string p_extra_info) {
	prior_info = prior_info + "\n" + p_extra_info;
}

std::string Prior::get_info() const {
	return prior_info;
}

