/*
 * OptimizedMechanism.h
 *
 *  Created on: Aug 16, 2024
 *      
 */

#ifndef MECHANISMS_PLAN_OPTIMIZEDMECHANISM_H_
#define MECHANISMS_PLAN_OPTIMIZEDMECHANISM_H_

#include "glpk.h"
#include "PMatrixMechanism.h"

class OptimizedMechanism: public PMatrixMechanism {
public:
	OptimizedMechanism(int p_hsize, int p_vsize);
	~OptimizedMechanism();   // destructor: destroys the pointer to the glpk linear programming problem lp_ptr

protected:

	glp_prob * lp_ptr;  // pointer to a glpk linear programming problem.

	// some private functions for debugging

	/*  display the variables of the linear programming problem */
	void showVars();

	/*  display the objective function of the linear programming problem */
	void showObjectiveFunction ();

	/* show all constraints of the linear programming problem */
	void showConstraints ();

	/* show the values of the problem variables in the linear programming problem */
	void showVarValues ();
};

#endif /* MECHANISMS_PLAN_OPTIMIZEDMECHANISM_H_ */
