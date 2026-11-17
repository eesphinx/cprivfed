/*
 * OptimizedMechanism.cpp
 *
 *  Created on: Aug 16, 2024
 *      Author: ehab
 */

#include <iostream>
#include "OptimizedMechanism.h"

OptimizedMechanism::OptimizedMechanism(int p_hsize, int p_vsize)
:PMatrixMechanism(p_hsize, p_vsize){
	lp_ptr = glp_create_prob();
	glp_create_index(lp_ptr);    // important to find the index of a column (variable) from its name using glp_find_col()
}

OptimizedMechanism::~OptimizedMechanism() {
	glp_delete_prob(lp_ptr);
	lp_ptr = nullptr;
}

/*  display the variables of the linear programming problem */
void OptimizedMechanism::showVars() {
	// to check the added columns
	for (int t=1; t<=glp_get_num_cols(lp_ptr); t++)
		std::cout << glp_get_col_name(lp_ptr, t) << "\n";
}


/*  display the objective function of the linear programming problem */
void OptimizedMechanism::showObjectiveFunction () {
	// Debug: printing objective
	int num_cols = glp_get_num_cols(lp_ptr);
	std::string objectiveStr = std::string(glp_get_obj_name(lp_ptr))+ "= ";
	for (int i=1; i<= num_cols;i++) {
		double coef = glp_get_obj_coef(lp_ptr,i);
		if (coef != 0)
			objectiveStr = objectiveStr + std::to_string(coef)+"*"+ glp_get_col_name(lp_ptr, i)+" + ";
	}
	std::cout << objectiveStr << std::endl;
}

/* show all constraints of the linear programming problem */
void OptimizedMechanism::showConstraints () {
	for (int i=1; i <= glp_get_num_rows(lp_ptr);i++) {
		int ind[1000]; double val[1000];
		int len = glp_get_mat_row(lp_ptr, i, ind, val);
		std::string conStr ="";
		for (int j=1; j<len; j++)
			conStr = conStr + std::to_string(val[j]) +"*"+ glp_get_col_name(lp_ptr, ind[j])+ "+";
		conStr = conStr + std::to_string(val[len]) +"*"+ glp_get_col_name(lp_ptr, ind[len]);

		int row_type = glp_get_row_type(lp_ptr,i);
		if (row_type == GLP_UP)
			conStr = conStr + " <= " + std::to_string(glp_get_row_ub(lp_ptr, i));
		else if (row_type == GLP_LO)
			conStr = conStr + " >= " + std::to_string(glp_get_row_lb(lp_ptr, i));
		else if (row_type == GLP_FX)
			conStr = conStr + " == " + std::to_string(glp_get_row_ub(lp_ptr, i));
		else if (row_type == GLP_DB)
			conStr = std::to_string(glp_get_row_lb(lp_ptr, i)) + "<= " + conStr + " <= " + std::to_string(glp_get_row_ub(lp_ptr, i));

		std::cout << conStr << "\n" << std::endl;
	}
}


void OptimizedMechanism::showVarValues () {
	int num_cols = glp_get_num_cols(lp_ptr);
	for (int i=1; i< num_cols;i++)
		std::cout << glp_get_col_name(lp_ptr, i) << " = " << glp_get_col_prim(lp_ptr, i) << "\n";
	std::cout << std::endl;
}

