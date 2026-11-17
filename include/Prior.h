/*
 * Prior.h
 *
 *  Created on: Jun 21, 2019
 *      
 */
#ifndef PRIV_PLAN_PRIORS_PRIOR_H_
#define PRIV_PLAN_PRIORS_PRIOR_H_

#include <string>
#include "DoubleRow.h"    // has to be included because Prior has a member of type DoubleMatrix
//(compiler has to know complete information about it)

class Prior {
public:
	int hsize;
	int vsize;


	/** constructs a prior on hsize*vsize grid using a given row vector (elements are in the linear order). */
	Prior(const DoubleRow& p_vect, int p_hsize, int p_vsize);

	/** For creating a prior from a file. It should include the matrix */
	Prior(std::string p_file);

	/**
	 * Creates a prior object from a file having the probabilities in the linear order
	 * i.e., the first probability in p_str corresponds to P(y=0,x=0)
	 */
	Prior(std::string p_file, int p_hsize, int p_vsize);


	/** return a uniform prior on the grid */
	static Prior uniform(int p_hsize, int p_vsize);

	/** return a random prior on the grid */
	static Prior random(int p_hsize, int p_vsize);


//	/** fill the prior matrix with random probabilities */
//	void fill_prior_random();

	/** write the prior_matrix to the given file. If p_with_info=true, the header information is written too. */
	void write_prior_as_matrix(std::string p_file, bool p_with_info) const;

	/** returns the prior data in a single string format as required by Kostas first row is written first, etc. */
	std::string toString () const;

	/** write the prior as one vector in the linear (row major) format */
	void write_prior_as_vector(std::string p_file) const;

	/** get the probability of cell p specified in the linear order (row first) */
	double getProb(int p) const;

	/** returns the row vector representation of the prior */
	DoubleRow get_prior_vector() const;

	void add_info(std::string p_extra_info);

	std::string get_info() const;

private:
	DoubleRow prior_vect;  // a matrix of one row
	std::string prior_info;


};




#endif /* PRIV_PLAN_PRIORS_PRIOR_H_ */
