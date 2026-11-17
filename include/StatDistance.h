/*
 * StatDistance.h
 *
 *  Created on: Jul 1, 2019
 *      
 */

#ifndef STAT_DISTRIBUTIONS_STATDISTANCE_H_
#define STAT_DISTRIBUTIONS_STATDISTANCE_H_

#include <vector>

class DoubleRow;
class AbstractGridSpace;
class Loss;

// statDistance types
enum StatDistanceType {TV, KL, EMD, SQD, RSQD};

class StatDistance {
public:

	/**
	 * returns a string representation for the given statistical distance type
	 */
	static std::string distTypeToString(StatDistanceType p_type);

	/** return the distance type corresponding to the given string */
	static StatDistanceType stringToDistType(std::string p_typeStr);

	/**
	 * returns L1-norm distance between two distributions.
	 * Note that the total variation distance TV is half of the returned distance
	 */
	static double getNorm1Dist(const DoubleRow& p_distribution1, const DoubleRow& p_distribution2);

	/**
	 * returns total variation distance between two distributions.
	 * Note that the total variation distance TV is half of the L1-norm distance (for discrete distributions).
	 */
	static double getTotalVariationDist(const DoubleRow& p_distribution1, const DoubleRow& p_distribution2);


	/** returns squared distance between two distributions */
	static double getSquaredDist(const DoubleRow& p_distribution1, const DoubleRow& p_distribution2);

	/** returns root squared distance between two distributions */
	static double getRootSquaredDist(const DoubleRow& p_distribution1, const DoubleRow& p_distribution2);

	/**
	 * returns the EMD distance between two distributions given a generic gridSecretSpace
	 * on which the distributions are defined and the loss object used to measure
	 * the loss between elements of p_xspace.
	 */
	static double getEMDDist(const DoubleRow& p_distribution1, const DoubleRow& p_distribution2, AbstractGridSpace* p_xspace_ptr, const Loss& p_LossObj);


	/**
	 * returns the EMD distance between two distributions given a generic gridSecretSpace.
	 * The loss in this case is assumed to be LinGeoLoss defined according to the size of
	 * the given secrets space. This procedure tries to use a cashed copy of the loss matrix
	 * if it exists, otherwise constructs it and store a copy.
	 *
	 * It is also effecient since
	 * it uses the cashed loss two dimensional vector.
	 */
	static double getGeoEMDDist(const DoubleRow& p_distribution1, const DoubleRow& p_distribution2,
			AbstractGridSpace* p_xspace_ptr);


	/**
	 * Returns the Symmetric Kullback-Liebler measure from the 1st distribution to the 2nd one.
	 * 	zero probabilities in the 1st distribution (the original) are not an issue because pi log pi
	 *  tends to 0, so they are handled by adding 0 for these components. However the zero components
	 *  in the second distributions make pi log qi diverge. In privacy this situation happens when
	 *  the cell size is small (e.g. 0.2km) or the number of samples is small (many cells do not include
	 *  noisy checkins). In this case, we smooth the KL measure when qi=0 by setting q = a.U + (1-a).q, where
	 *  U is the uniform distribution and a is a small probability such that a << 1/grid_size to keep
	 *  the resulting distribution far from the uniform distribution in which qi=1/grid_size
	 */
	static double getSymKLDiff(const DoubleRow& p_distribution1, const DoubleRow& p_distribution2);

	/** return the distance the original distribution depending on the distance type */
	static double getDist(const DoubleRow& p_distribution1, const DoubleRow& p_distribution2, StatDistanceType p_distType, AbstractGridSpace* xspace_ptr);

private:

	/*
	 * returns a two dimensional cost vector that corresponds to the given loss object and generic
	 * grid space.
	 */
	static std::vector<std::vector<double> > getLossTwoDimVectorForSpace(AbstractGridSpace* p_xspace_ptr, const Loss& p_LossObj);

	/*
	 * returns a two dimensional cost vector that corresponds to LinGeoLoss object constructed on the
	 * given generic space. This function tries to return a cashed copy that is compatible with the
	 * given space. If not found this two dimensional vector is constructed, cashed, and then returned.
	 */
	static std::vector<std::vector<double> > getLinGeoLossTwoDimVectorForSpace(AbstractGridSpace* p_xspace_ptr);


	/*
	 * returns a two dimensional cost vector that corresponds to LinGeoLoss on a grid.
	 * It checks if it exists it returns it otherwise creates it
	 */
	static void setLossTwoDimVectorForSpace(AbstractGridSpace* p_xspace);

	static std::vector<std::vector<double> > cashedLinGeoLossTwoDimVector;



};


#endif /* STAT_DISTRIBUTIONS_STATDISTANCE_H_ */
