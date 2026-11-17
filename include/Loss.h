/*
 * Loss.h
 *
 *  Created on: Jun 22, 2019
 *      
 */

#ifndef PRIV_PLAN_LOSS_LOSS_H_
#define PRIV_PLAN_LOSS_LOSS_H_


/** This class specifies the necessary function (get_loss) that computes
    the loss based on the coordinates of given cells. Every loss function
    is assumed to be increasing. This function behaves as follows.
    loss = get_loss_by_dist(distance) For all 0<=distance<= saturation_dist,
    loss = max_loss for all saturation_dist < dist.
 */

// loss types
enum lossType {LOSS_BIN, LOSS_LIN_GEO, LOSS_LIN_MAN};

class DoubleMatrix;     // forward declaration is enough because we use it as parameter
class Loss {

public:

	Loss(double p_saturation_dist, double p_max_loss);
	// make a virtual destructor in case we delete a Loss object using a pointer, so the proper derived destructor is called
	virtual ~Loss(){};

	/**
	 * Loss for two points i,j: gets the loss for input i=(i_x, i_y) and output j=(j_x,j_y)
	 **/
	virtual double get_loss(double i_x, double i_y, double j_x, double j_y) const = 0;


	/* Loss for a distance:
	 * This function takes the distance p_dist between two points
	 * and return the loss due to this distance.
	 * The semantics of the distance depends on the type of the loss */
	virtual double get_loss_by_dist(double p_dist) const = 0;    // pure virtual function: must be implemented by subclasses

	// get a descriptive name for the loss function
	virtual lossType getType() const = 0;

	double get_saturation_dist() const;

	double get_max_loss() const;

	/**
	 * for a grid with dimensions hsize,vsize, cell_side returns a matrix
	 * representing the loss for every pair of input and output cells.
	 **/
	DoubleMatrix get_loss_matrix(int p_hsize, int p_vsize, double p_cell_side) const;

protected:  // accessable by derived classes

	double saturation_dist;  // default max distance which get_loss() respects.
	double max_loss;    // loss for distances larger than saturation_dist

};



#endif /* PRIV_PLAN_LOSS_LOSS_H_ */
