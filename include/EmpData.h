/*
 * EmpData.h
 *
 *  Created on: Sep 23, 2019
 *      
 */

#ifndef SRC_STAT_EMPIRICAL_EMPDATA_H_
#define SRC_STAT_EMPIRICAL_EMPDATA_H_

#include <vector>
#include "EmpDatum.h"

class DoubleRow;
class AbstractOutputLocalMechanism;
class GridEmpData;

class EmpData {
public:
	/** a vector of Observation/frequency where frequencies must be nonzero **/
	std::vector<EmpDatum> data;

	/** default constructor invoked automatically when an object is declared */
	EmpData(){};

	/** constructs an EmpData by obfuscating (raw) empirical data using 
	 * a given AbstractOutputLocalMechanism. <p_raw_empdata> is of the general type 
	 * EmpData to allow all subclasses of Empdata including for instances GridEmpData 
	 * and UCIBinAgeData, etc. p_raw_empdata is checked to be raw by inspecting its
	 * associated mechanism: it must be the identity mechanism (IDENT) 
	 * Note: Only local mechanism are allowed to generate EmpData */
	EmpData(const EmpData& p_raw_empdata, const AbstractOutputLocalMechanism* p_local_mech_ptr);

	/** constructs an Empdata object from a data file. */
	EmpData(std::string& p_data_file);

	/** add noisy empirical data p_empdata (they must have been obfuscated by the same 
	 * mechanism associated to the existing data) */
	void add(const EmpData& p_empdata);

	/** Obfuscate raw data p_raw_empdata using the mechanism of the current 
	 *  EmpData object, and add the noisy data to the existing data 
	 * given data are of type EmpData to be general to any raw data
	 * Raw data are distinguished by being associated with the identity mechanism (IDENT)*/
	void obfuscateAndAdd(const EmpData& p_raw_empdata);

	/** write the data to the given file */
	void write_data_to_file(const std::string& p_file) const;

	/** split the data into multiple empirical data objects according to the given probability distribution */
	std::vector<EmpData> randSplit(const DoubleRow& p_split_distribution) const;

	static std::vector<EmpData> compress(const std::vector<EmpData>& p_list_empdata );
	
	/**
	 * returns the empirical distribution on observables as a row vector. The j-th element corresponds
	 * to the j-th datum (observation , frequency) in data. P.S. Since it is constructed from data, it
	 * contains no zero entries (i.e. excludes observations that do not appear in the data).
	 * It is used in the IBU estimation since it requires no zeros in the empirical distribution.
	 */
	DoubleRow getDistributionOnObservables() const;

	/** returns the number of observed outcomes with nonzero frequencies **/
	size_t getSize() const;

	/** return the sum of all frequencies in data */
	int get_totalObservations() const;


	/** returns the mechanism info string */
	std::string getObfMechInfo() const;

	/** returns a pointer to a clone of the mechanism that generated the data: Important delete
	 * the object in the calling procedure after use */
	//AbstractOutputMechanism* getObfMech() const;

	/**  specific to BitVector observations data */

	/**
	 * This returns a row vector with length equal to the number of bits in one bit vector. The
	 * j-th element is Tj/n where Tj is the number of bitVectors that have bit j equal to 1, and
	 * n is the total number of observed bit vectors in the data (totalObservations). Note the
	 * resulting vector is not a distribution.
	 * This function is used for Rappor estimators BOTRappor_INVN and BOTRappor_INVP.  */
	DoubleRow get_relative_frequency_of_bits() const;


protected: // accessed only by derived classes
	/** information about the mechanism that generated the data */
	std::string mechanism_info;

};



#endif /* SRC_STAT_EMPIRICAL_EMPDATA_H_ */
