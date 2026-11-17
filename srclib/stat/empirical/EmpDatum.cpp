/*
 * EmpDatum.cpp
 *
 *  Created on: Sep 23, 2019
 *      Author: ehab
 */
#include <iostream>
#include <sstream>
#include "EmpDatum.h"
#include "AbstractOutputLocalMechanism.h"
#include "Observation.h"

/** construct observation-frequency datum */
EmpDatum::EmpDatum(Observation* p_obs_ptr, int p_frequency)
:frequency(p_frequency) {
	obs_ptr = p_obs_ptr->clone();
}

/** construct a datum from a given data line in the format: Observation | frequency (as written in files)
 * Note that the mechanism p_mech_ptr is necessary to parse the observation. This mechanism may be local or central */
EmpDatum::EmpDatum(std::string p_data_line, const AbstractOutputLocalMechanism* p_local_mech_ptr)
:frequency(0) {
    std::istringstream iss(p_data_line);
    std::string obs_str;
    std::string token;
    while (iss >> token && token != "|")
    	obs_str = obs_str + token + ' ';
    if (token != "|") {
    	std::cerr << "unable to construct datum from data line: invalid format";
    	exit(0);
    }
    if (iss >> frequency && frequency != 0)  // set the frequency field of the datum
    	// set the observation pointer field to a pointer to newly created observation
    	obs_ptr = p_local_mech_ptr->readObservation(obs_str);
    else {
    	std::cerr << "unable to construct datum from data line: invalid frequency";
    	exit(0);
    }
}

/** copy constructor */
EmpDatum::EmpDatum(const EmpDatum& other)
: frequency(other.frequency) {
	//std::cout << "EmpDatum copy constructor for" << other.toString() << "\n";
	obs_ptr = other.obs_ptr->clone();
}

/** move constructor: used by the compiler to move rvalue EmpDatum (other) to this */
EmpDatum::EmpDatum(EmpDatum&& other)
: frequency(other.frequency) {
	//std::cout << "EmpDatum move constructor for" << other.toString() << "\n";
	obs_ptr = other.obs_ptr;
	other.obs_ptr = nullptr;
}

/** The assignment operator, used for assigning objects. Works for both copy and move assignments.
 * Important:
 * Unlike usual implementations, in this function we pass other by value:
 * The true genius of this method is what happens when the passed-by-value object to be assigned from
 * (the source object) is created as a temporary object in the assignment operator function. There are
 * three possibilities:
 *  i.  The source object was an lvalue: the source object’s copy constructor will be called to make a copy
 *      and the overload will behave as a copy assignment operator.
 *  ii. The source object was an rvalue: nothing is called; the compiler elides the temporary copy and passes
 *      the object by value; the overload will behave as a move assignment operator.
 *  iii.The source object was an rvalue: The source object’s move constructor will be called and the overload
 *      will behave as a move assignment operator. */
EmpDatum& EmpDatum::operator=(EmpDatum other) {
	//todo:debug
	// todo: do we need to check if the underlying observations have the same type ?
	//std::cout << "assignment operation\n";
	frequency = other.frequency;
	std::swap(obs_ptr, other.obs_ptr);   // safe swap of the value of pointers obs_ptr and other.obs_ptr, hence no data copy is done.
	return *this;
}


/** destructor: to delete the observation */
EmpDatum::~EmpDatum() {
	delete(obs_ptr);
	obs_ptr = nullptr;
}

/** return a string representation for an EmpDatum in the form [observation | frequency] */
std::string EmpDatum::toString() const {
	return obs_ptr->toString() + " | " + std::to_string(frequency);
}



