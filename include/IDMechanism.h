/*
 * IDMechanism.h
 *
 *  Created on: Sep 28, 2019
 *      
 */

/** identity local mechanism which maps every location to itself with probability */

#ifndef SRC_MECHANISMS_PLAN_IDMECHANISM_H_
#define SRC_MECHANISMS_PLAN_IDMECHANISM_H_

#include "GridOutputLocalMechanism.h"


class IDMechanism: public GridOutputLocalMechanism {
public:
	/* This mechanism is defined only by the space dimensions	 */
	IDMechanism(int p_hsize, int p_vsize);

	IDMechanism(const std::string p_info);

	std::string get_info() const override;
	std::string getTypeAndPriv(const int p_precision) const override;

private:
	double getProb(int p_in, int p_out) const override;
};



#endif /* SRC_MECHANISMS_PLAN_IDMECHANISM_H_ */
