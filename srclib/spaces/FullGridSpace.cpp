#include <iostream>
#include <sstream>
#include "FullGridSpace.h"
#include "EmpData.h"
#include "GridObservation.h"
#include "DoubleRow.h"


FullGridSpace::FullGridSpace(int p_hsize, int p_vsize, double p_cell_side)
: AbstractGridSpace(p_hsize, p_vsize, p_cell_side) {
}

/**
 * get the uniform distribution (as row vector) on the set of secrets of the space.
 * the size of this vector and the mapping from the grid cells to its entries is
 * maintained by the specific class. Remember the the caller deletes he pointer when
 * it is done with it
 */
DoubleRow FullGridSpace::get_uniform_distribution() const {
	double p = 1.0/(hsize*vsize);
	return DoubleRow(hsize*vsize).fillWith(p);
}


/** return the linear indexes (in S order) of grid cells in the space. This depends on the specifics of the space */
std::vector<int> FullGridSpace::getSecretsIndexes() const {
	std::vector<int> result;
	for (int i=0; i<hsize*vsize; i++)
		result.push_back(i);
	return result;
}

/** returns the number of locations (cells) in the space */
int FullGridSpace::getSize() const {
	return hsize*vsize;
}


DoubleRow FullGridSpace::get_histogram_from_emp_data(const EmpData& p_empdata) const {
	DoubleRow retHistogram(hsize*vsize);
	retHistogram.fillWith(0);
	for (size_t i=0; i<p_empdata.getSize();i++)
		if (GridObservation* obs_ptr = dynamic_cast<GridObservation*>(p_empdata.data[i].obs_ptr)) {
			// we colud also check if the location in the space or not.
			retHistogram.put(obs_ptr->location, p_empdata.data[i].frequency);
		}
		else {
			std::cerr << "unable to construct histogram from empirical data. Invalid format";
			exit(1);
		}
	return retHistogram;
}


/**
 * the distribution on the secrets (in this case it is a full grid).
 * It takes a GridEmpData object and produces from it a distribution on the secrets.
 * Note that it may contain zero elements. The distribution respects the semantics
 * of the space secrets: i-th item is location on the hsize*vsize grid represented by
 * its linear order. This method is used to obtain the original distribution on secrets
 * from empirical data. It is also used in the matrix inversion method to obtain the
 * (full) noisy distribution.
 **/
DoubleRow FullGridSpace::get_distribution_from_emp_data(const EmpData& p_empdata) const {
	DoubleRow hist =  get_histogram_from_emp_data(p_empdata);
	double sum = hist.sum();
	return hist.div(sum);
}

/** return max distance between two cells: used to get the max loss for EMD  */
double FullGridSpace::getMaxDistance() const {
	return cell_side * sqrt(hsize*hsize + vsize*vsize);
}



/**
 * returns a row vector, the distribution on the secrets as a row vector.
 * takes a BitVectorEmpData object (resulting from a Rappor mechanism) and produces
 * from it a distribution on the secrets. Note that it may contain zero elements.
 * This method is used to obtain the empirical distribution on secrets from a bitVector
 * obfuscated data.
 **/
/*
 *   not applicable since noisy data are not defined by bit vectors.
DoubleRow FullGridSpace::get_distribution_from_bitVectEmpData(const BitVectorEmpData& p_bitVectorEmpData) {
	DoubleRow retDistribution(hsize*vsize);
	retDistribution.fillWith(0);
	int numObservedVectors = p_bitVectorEmpData.data.size();
	for (int i=0; i< numObservedVectors; i++)
		for (int j=0; j<hsize*vsize; j++)
			if (p_bitVectorEmpData.data.at(i).bitVector[j])
				retDistribution.put(j, retDistribution.get(j) + p_bitVectorEmpData.data.at(i).frequency);
	double sum = retDistribution.sum();
	retDistribution.divi(sum);
	return retDistribution;
}
*/

