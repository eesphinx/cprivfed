/*
 * GridSecretsSpace.h
 *
 *  Created on: Jul 1, 2019
 *      
 */

#ifndef SPACES_FULLGRIDSPACE_H_
#define SPACES_FULLGRIDSPACE_H_

#include <cmath>
#include "AbstractGridSpace.h"

class EmpData;
class BitVectorEmpData;

/**
 * This class describes the space X of secrets as the entire grid of locations specified by hsize,vsize,
 */
class FullGridSpace: public AbstractGridSpace {
public:
	FullGridSpace(int p_hsize, int p_vsize, double p_cell_side);

	/**
	 * get the row vector describing a uniform prior on the set of secrets.
	 * the size of this vector and the mapping from the grid cells to its entries is
	 * maintained by the specific class.
	 */
	DoubleRow get_uniform_distribution() const;


	/** return the linear indexes (in S order) of grid cells in the space. This depends on the specifics of the space */
	std::vector<int> getSecretsIndexes() const;

	int getSize() const;

	/** return max distance between two cells: used to get the max loss for EMD  */
	double getMaxDistance() const;



	/**
	 * returns a row vector, that is the distribution on the secrets.
	 * It takes a GridEmpData object and produces from it a distribution on the full grid.
	 * Note that it may contain zero elements. The distribution respects the semantics
	 * of the space secrets: i-th item is location on the hsize*vsize grid represented by
	 * its linear order. This method is used to obtain the original distribution on secrets
	 * from empirical data. It is also used in the matrix inversion method
	 **/
	DoubleRow get_distribution_from_emp_data(const EmpData& p_empdata) const;

	DoubleRow get_histogram_from_emp_data(const EmpData& p_empdata) const;
	
	
	/**
	 * returns a row vector, the distribution on the secrets as a row vector.
	 * takes a BitVectorEmpData object (resulting from a Rappor mechanism) and produces
	 * from it a distribution on the secrets. Note that it may contain zero elements.
	 * This method is used to obtain the empirical distribution on secrets from a bitVector
	 * obfuscated data.
	 **/
	//  This does not reflect noisy data.
	//DoubleRow get_distribution_from_bitVectEmpData(const BitVectorEmpData& p_bitVectorEmpData);


};




#endif /* SPACES_FULLGRIDSPACE_H_ */
