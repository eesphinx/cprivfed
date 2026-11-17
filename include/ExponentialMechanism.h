/*
 * ExponentialMechanism.h
 *
 *  Created on: Jul 8, 2019
 *      
 */

#ifndef PRIV_PLAN_MECHANISMS_EXPONENTIALMECHANISM_H_
#define PRIV_PLAN_MECHANISMS_EXPONENTIALMECHANISM_H_

#include "GridOutputLocalMechanism.h"
#include "DoubleColumn.h"
class PrivMetric;

class ExponentialMechanism : public GridOutputLocalMechanism {
public:

	double cell_side;
	/**
	 * This mechanism is parameterize by a generic privacy metric. Therefore the constructor receives
	 * a pointer to the (abstract) privObject. The object it cloned and a pointer to the clone is stored.
	 * This is to avoid the destruction of the original privObj. The clone is deleted in the destructor.
	 * Since privacy metric depends on the distance, cell_side must be specified together with hsize, vsize.
	 * These information is obtained from the passed FullGridSpace.
	 */
	ExponentialMechanism(int p_hsize, int p_vsize, double p_cell_side, PrivMetric* p_privObjPtr);

	ExponentialMechanism(const FullGridSpace& p_space, PrivMetric* p_privObjPtr);

	ExponentialMechanism(const std::string p_info);

	/** destructor: to free the dynamically allocated the normailzation constants array */
	~ExponentialMechanism();

	std::string get_info() const override;
	std::string getTypeAndPriv(const int p_precision) const override;

	PrivMetric* privObjPtr;   // This is abstract class, so include only a pointer. This points to a clone that remains alive with the mechanism

private:
	double getProb(int p_in, int p_out) const override;

	/**
	 * A column of constants of size hsize*vsize. Every constant corresponds to a mechanism row.
	 **/
	DoubleColumn mechanism_constants;

	/** sets the normalization constant for every input location (x,y) the mechanism in the linear order y*hsize+x */
	void construct_norm_constants();
};




#endif /* PRIV_PLAN_MECHANISMS_EXPONENTIALMECHANISM_H_ */
