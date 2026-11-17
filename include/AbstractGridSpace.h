/*
 * SecretsSpace.h
 *
 *  Created on: Jul 1, 2019
 *      
 */

#ifndef SPACES_ABSTRACTGRIDSPACE_H_
#define SPACES_ABSTRACTGRIDSPACE_H_

#include <vector>

class DoubleRow;
class EmpData;

/**
 * Describes an abstract space X of secrets, as a general subset of a grid of locations
 * The grid is specified by hsize,vsize, cell_side. The subset may be the entire grid as
 * in (FullGridSecretsSpace) or a discrete set of cells, or
 * a convex hull of locations. Every case of these is implemented by a
 * subclass of AbstractGridSpace.
 */
class AbstractGridSpace {
public:
	// The following are essential information about the grid. The space is defined as the centers of the cells.
	// The cell side is important when applying d-privacy mechanisms since the privacy on the space
	// depends on the distance between centers of the cell. However, the cell side is not relevant when applying
	// local differential privacy mechanisms, e.g. kRR and RAPPOR. The cell_side is also important to determine
	// the Euclidean distance between two points which is used in computing the EMD distance implemented in the
	// StatDistance class.

	int hsize;
	int vsize;
	double cell_side;

	AbstractGridSpace(int p_hsize, int p_vsize, double p_cell_side);

	/** make a virtual destructor in case we delete a object using a pointer,
	 * so the most derived destructor is called */
	virtual ~AbstractGridSpace(){};

	/**
	 * get the row vector describing a uniform prior on the set of secrets.
	 * the size of this vector and the mapping from the grid cells to its entries is
	 * maintained by the specific class.
	 */
	virtual DoubleRow get_uniform_distribution() const = 0;  // pure virtual to force overriding

	/**
	 * returns the distribution on the secrets, given empirical data (the empirical data should include only observables from
	 * the secret space). Note that it may contain zero elements. The distribution respects the semantics of the space secrets.
	 * Therefore it has to be implemented by every subclass. This method is used to obtain the original distribution on secrets
	 * from empirical data. It is also used in the matrix inversion method.
	 **/
	virtual DoubleRow get_distribution_from_emp_data(const EmpData& p_empdata) const=0;

	/** return the linear indexes (in S order) of grid cells in the space. This depends on the specifics of the space */
	virtual std::vector<int> getSecretsIndexes() const = 0;

	/** returns the number of locations in the soace */
	virtual int getSize() const=0;
	/**
	 * for future implementation:
	 * returns a row vector, that is a distribution on the secrets. It takes a GridEmpData object
	 * that contains data drawn from (and hence is subset of) the secret space, and produces from
	 * it a distribution on the secrets. Note that it may contain zero elements. The mapping between
	 * the components in this distribution and the secrets depends on the definition of the secret
	 * space, i.e. is obtained via the subclass. For example in the subclass FullGridSecretsSpace,
	 * this method is used to obtain the original distribution on secrets from empirical data constructed
	 * from check-ins. Also in this case i-th element in the resulting distribution corresponds to the
	 * i-th cell of the full grid in  the linear order. It is also used in the matrix inversion
	 * method which is defined only for the subclass FullGridSecretsSpace.
	 * We remove it in the current implementation because the above purposes are needed only for the fullGridSpace
	 * Hence we implement it only for that subclass
	 **/
	//virtual DoubleRow get_distribution_from_emp_data(const GridEmpData& p_empdata) const = 0;


	/** return the continuous (double) x coordinate for a given location in the linear order (linear_index)
	 * This is used for computing the loss between two points, hence computing the EMD distance */
	double getx(int linear_index) const;

	/** return the constinuous (double) y coordinate for a given location in the linear order (linear_index) */
	double gety(int linear_index) const;

	/** return max distance between two cells: used to get the max loss for EMD
	 * may be implemented by subclasses */
	virtual double getMaxDistance() const =0 ;


};




#endif /* SPACES_ABSTRACTGRIDSPACE_H_ */
