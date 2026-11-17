#include <AbstractGridSpace.h>
#include <vector>
#include <iostream>
#include <cmath>
#include "DoubleRow.h"
#include "StatDistance.h"
#include "emd_hat.hpp"
#include "Loss.h"
#include "LinGeoLoss.h"

std::vector< std::vector<double> > StatDistance::cashedLinGeoLossTwoDimVector {{-1}};  // this definition is necessary for every static variable

/**
 * returns a string representation for the given statistical distance type
 */
std::string StatDistance::distTypeToString(StatDistanceType p_type) {
	switch(p_type) {
	case TV : return "TV";
	case KL : return "KL";
	case EMD : return "EMD";
	case SQD : return "SQD";    // squared distance
	case RSQD : return "RSQD";    // root squared distance
	default :
		std::cerr << "invalid statistical distance type \n";
		exit(1);
	}
	return "";
}

StatDistanceType StatDistance::stringToDistType(std::string p_typeStr) {
	if (p_typeStr == "TV") return TV;
	if (p_typeStr == "KL") return KL;
	if (p_typeStr == "EMD") return EMD;
	if (p_typeStr == "SQD") return SQD;
	if (p_typeStr == "RSQD") return RSQD;
	std::cerr << "invalid statistical distance type string\n";
	exit(1);
}


/**
 * returns L1-norm distance between two distributions. Note that the total variation distance TV is half of the returned distance
 */
double StatDistance::getNorm1Dist(const DoubleRow& p_distribution1, const DoubleRow& p_distribution2) {
	return p_distribution1.distance1(p_distribution2);
}

/**
 * returns total variation distance between two distributions.
 * Note that the total variation distance TV is half of the L1-norm distance (for discrete distributions).
 */
double StatDistance::getTotalVariationDist(const DoubleRow& p_distribution1, const DoubleRow& p_distribution2) {
	return  0.5*p_distribution1.distance1(p_distribution2);
}

/** returns squared distance between two distributions */
double StatDistance::getSquaredDist(const DoubleRow& p_distribution1, const DoubleRow& p_distribution2) {
	return  p_distribution1.squaredDistance(p_distribution2);
}

/** returns root squared distance between two distributions */
double StatDistance::getRootSquaredDist(const DoubleRow& p_distribution1, const DoubleRow& p_distribution2) {
	return  sqrt(p_distribution1.squaredDistance(p_distribution2));
}



/*
 * returns a two dimensional cost vector that corresponds to the given loss object and generic
 * grid space.
 */
std::vector<std::vector<double> > StatDistance::getLossTwoDimVectorForSpace(AbstractGridSpace* p_xspace_ptr, const Loss& p_LossObj) {
	int xspace_size = p_xspace_ptr->getSize();
	// construct the loss twoDimVector
	std::vector< std::vector<double> > loss_twoDimVector;
	// declare and construct a row vector of size (length)
	std::vector<double> loss_row(xspace_size);
	// construct the twoDimVector by adding row by row
	for (int i=0; i<xspace_size; ++i) loss_twoDimVector.push_back(loss_row);
	// fill the values of the matrix by setting loss_twoDimVector[i][j].
	std::vector<int> xspace_indexes = p_xspace_ptr->getSecretsIndexes();
	for (int i=0;i<xspace_size;i++)
		for (int j=0;j<xspace_size;j++)
			loss_twoDimVector[i][j] = p_LossObj.get_loss(p_xspace_ptr->getx(xspace_indexes[i]), p_xspace_ptr->gety(xspace_indexes[i]),
					p_xspace_ptr->getx(xspace_indexes[j]), p_xspace_ptr->gety(xspace_indexes[j]));
	return loss_twoDimVector;
}

/*
 * returns a two dimensional cost vector that corresponds to LinGeoLoss object constructed on the
 * given generic space. This function tries to return a cashed copy that is compatible with the
 * given space. If not found this two dimensional vector is constructed, cashed, and then returned.
 */
std::vector<std::vector<double> > StatDistance::getLinGeoLossTwoDimVectorForSpace(AbstractGridSpace* p_xspace_ptr) {
	//	int xspace_size = p_xspace_ptr->getSize();
	if (cashedLinGeoLossTwoDimVector[0][0] > -1
			&& cashedLinGeoLossTwoDimVector.size() ==p_xspace_ptr -> getSize())
		return cashedLinGeoLossTwoDimVector;

	double max_dist = p_xspace_ptr->getMaxDistance();
	LinGeoLoss loss(max_dist,max_dist);
	std::vector<std::vector<double> > linGeoLossTwoDimVector = getLossTwoDimVectorForSpace(p_xspace_ptr, loss);

	cashedLinGeoLossTwoDimVector = linGeoLossTwoDimVector;
	return linGeoLossTwoDimVector;
}



/**
 * returns the EMD distance between two distributions given a generic gridSecretSpace
 * on which the distributions are defined and the loss object used to measure
 * the loss between elements of p_xspace.
 * This is the most flexible function.
 */
double StatDistance::getEMDDist(const DoubleRow& p_distribution1, const DoubleRow& p_distribution2,
		AbstractGridSpace* p_xspace_ptr, const Loss& p_LossObj) {
	if (p_distribution1.getLength() != p_distribution2.getLength()) {
		std::cerr << "unable to compute emd distance: the two distributions must have the same size\n";
		exit(0);
	}
	if (p_xspace_ptr ->getSize() != p_distribution1.getLength()) {
		std::cerr << "unable to compute emd distance: the distributions must have the same size as the given space\n";
		exit(0);
	}

	// construct the two distributions in the format vector<double>
	std::vector<double> distribution_vector1 = p_distribution1.toVector();
	std::vector<double> distribution_vector2 = p_distribution2.toVector();
	std::vector<std::vector<double> > lossTwoDimVector = StatDistance::getLossTwoDimVectorForSpace(p_xspace_ptr, p_LossObj);

	double ret_distance = emd_hat_gd_metric<double>()(distribution_vector1,distribution_vector2, lossTwoDimVector);
	return ret_distance;
}


/**
 * returns the EMD distance between two distributions given a generic gridSecretSpace.
 * The loss in this case is assumed to be LinGeoLoss. It is also effecient since
 * it uses the cashed loss two dimensional vector.
 */
double StatDistance::getGeoEMDDist(const DoubleRow& p_distribution1, const DoubleRow& p_distribution2,
		AbstractGridSpace* p_xspace_ptr) {

	if (p_distribution1.getLength() != p_distribution2.getLength()) {
		std::cerr << "unable to compute emd distance: the two distributions must have the same size\n";
		exit(0);
	}
	if (p_xspace_ptr ->getSize() != p_distribution1.getLength()) {
		std::cerr << "unable to compute emd distance: the distributions must have the same size as the given space\n";
		exit(0);
	}

	// construct the two distributions in the format vector<double>
	std::vector<double> distribution_vector1 = p_distribution1.toVector();
	std::vector<double> distribution_vector2 = p_distribution2.toVector();
	std::vector<std::vector<double> > lossTwoDimVector = StatDistance::getLinGeoLossTwoDimVectorForSpace(p_xspace_ptr);

	double ret_distance = emd_hat_gd_metric<double>()(distribution_vector1,distribution_vector2, lossTwoDimVector);
	return ret_distance;
}


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
double StatDistance::getSymKLDiff(const DoubleRow& p_distribution1, const DoubleRow& p_distribution2) {
	double sum=0;
	double p; double q;
	double a = 1e-8;  // smoothing factor. It should be much smaller than uniform probability, i.e. 1/grid_size
	//		System.out.println("distribution sums = "+ p_distribution1.sum() +", "+ p_distribution2.sum());
	if (p_distribution1.getLength() != p_distribution2.getLength()) {
		std::cerr << "the given distributions must have the same length\n";
		exit(0);
	}
	int l = p_distribution1.getLength();
	DoubleRow distribution1 = p_distribution1;
	DoubleRow distribution2 = p_distribution2;
	// smooth distribution1
	if (distribution1.min()< a) {  // smoothing
		DoubleRow tempRow(l);
		tempRow.fillWith(1.0).divi(l).muli(a);
		distribution1.muli(1-a).addi(tempRow);
	}

	// smooth distribution2
	if (distribution2.min()< a) {  // smoothing
		DoubleRow tempRow(l);
		tempRow.fillWith(1.0).divi(l).muli(a);
		distribution2.muli(1-a).addi(tempRow);
	}

	for (int i=0;i<l;i++) {
		p = distribution1.get(i);
		q = distribution2.get(i);
		if (p>0)
			sum = sum + p * log(p/q)/log(2); //returns the result in bits;
		if (q>0)
			sum = sum + q * log(q/p)/log(2); //returns the result in bits;
	}
	return sum;
}


/** return the distance the original distribution depending on the distance type */
double StatDistance::getDist(const DoubleRow& p_distribution1, const DoubleRow& p_distribution2, StatDistanceType p_distType, AbstractGridSpace* xspace_ptr = NULL) {
	double dist=0;
	if (p_distType == TV)
		dist = StatDistance::getTotalVariationDist(p_distribution1,p_distribution2);
	else if (p_distType == KL)
		dist = StatDistance::getSymKLDiff(p_distribution1, p_distribution2);
	else if (p_distType == SQD)
		dist = StatDistance::getSquaredDist(p_distribution1, p_distribution2);
	else if (p_distType == RSQD)
		dist = StatDistance::getRootSquaredDist(p_distribution1, p_distribution2);
	else if (p_distType == EMD) {
		if (xspace_ptr == nullptr || xspace_ptr->hsize == 0 || xspace_ptr->vsize==0 || xspace_ptr->cell_side ==0) {
			std::cerr << "unable to return EMD: the given xspace is invalid";
			exit(0);
		}
		dist = StatDistance::getGeoEMDDist(p_distribution1, p_distribution2, xspace_ptr);
	}
	else {
		std::cerr << "invalid distance type";;
		exit(0);
	}
	return dist;
}

