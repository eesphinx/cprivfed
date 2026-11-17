/*
 * Min_Loss_Given_dp_Mechanism.cpp
 *
 *  Created on: Feb 1, 2020
 *      Author: ehab
 */
#include <iostream>
#include <string>
#include <cmath>
#include "Prior.h"
#include "PrivMetric.h"
#include "Loss.h"
#include "FullGridSpace.h"
#include "Min_Loss_Given_dp_Mechanism.h"

Min_Loss_Given_dp_Mechanism::Min_Loss_Given_dp_Mechanism(const Prior& p_prior, const FullGridSpace& p_xspace,
		const PrivMetric* p_privMetric, const Loss* p_lossMetric, const std::string p_mechFile)
:OptimizedMechanism(p_xspace.hsize, p_xspace.vsize) {

	int ncells = hsize*vsize; // number of grid cells, also the number of columns in the linear programming problem
	const double MAX_PRIV = 1e+4;

	glp_set_obj_dir(lp_ptr, GLP_MIN);   // We want to minimize the objective function
	glp_set_obj_name(lp_ptr, "eloss");  // set the name of the objective function


	// Add the columns (variables) of the linear programming problem and construct the objective function
    std::cout << "Adding columns (variables) and constructing the objective function" << "\n" << std::endl;

	glp_add_cols(lp_ptr, ncells*ncells);
	std::cout << "No. of variables in the LP = " << glp_get_num_cols(lp_ptr) << "\n";
    int col_index=0;    //
    for (int i=0; i<ncells; i++)
    	for (int k=0; k<ncells; k++) {
    		col_index++;
    		std::string s = "m"+std::to_string(i) + "_" + std::to_string(k);
    		glp_set_col_name(lp_ptr, col_index, s.c_str());
    		glp_set_col_bnds(lp_ptr, col_index, GLP_DB, 0, 1);
    		double loss = p_lossMetric->get_loss(p_xspace.getx(i), p_xspace.gety(i), p_xspace.getx(k), p_xspace.gety(k));
    		double prob = p_prior.getProb(i);
    		glp_set_obj_coef(lp_ptr, col_index, prob*loss);
    	}


	  std::cout << "Adding inequality constraints (privacy) \n" << std::endl; //every one is:  m_ik - priv_factor * m_jk < 0

    for (int k=0; k<ncells; k++)   // for every output cell
    	for (int i=0; i<ncells; i++)  // an input cell
    		for (int j=0;j<ncells; j++) // another input cell
    			if (i != j) {
    				double priv_factor = exp(p_privMetric->get_priv_distance(p_xspace.getx(i), p_xspace.gety(i), p_xspace.getx(j), p_xspace.gety(j)));
    				if (priv_factor < MAX_PRIV)  {
    					std::string s = "m"+std::to_string(i)+"_"+std::to_string(k);
    					int left_var_ind = glp_find_col(lp_ptr, s.c_str());

    					s= "m"+std::to_string(j)+"_"+std::to_string(k);
    					int right_var_ind = glp_find_col(lp_ptr, s.c_str());

    					// specify the terms of the constraints at positions 1,2,...
    					int ind[3]; ind[1]=left_var_ind; ind[2]= right_var_ind; // never use position 0
    					double val[3]; val[1] = 1; val[2]= -1*priv_factor;  // never use position 0

    					int row_ind = glp_add_rows(lp_ptr, 1);  // add one new row (constraint)
    					glp_set_row_bnds(lp_ptr, row_ind, GLP_UP, 0, 0); // the expression in the constraint is upper bounded by 0
    					// add the constraint to the problem: the 3rd argument is the number of non-zero terms in the expression
    					glp_set_mat_row(lp_ptr, row_ind, 2, ind, val);

    				}
    			}


    std::cout << "Adding equality constraints \n" <<std::endl;

    for (int i=0;i<ncells;i++) {
			int row_ind = glp_add_rows(lp_ptr, 1);  // add a new constraint;
			glp_set_row_bnds(lp_ptr, row_ind, GLP_FX, 1, 1); // specify it is equality
			int ind[ncells+1]; double val[ncells+1];
			for (int k=0;k<ncells;k++) {  // for every output
				std::string s= "m"+ std::to_string(i)+ "_" +std::to_string(k) ;
				int var_ind = glp_find_col(lp_ptr, s.c_str());
				ind[k+1] = var_ind;
				val[k+1] =1;
			}
			glp_set_mat_row(lp_ptr, row_ind, ncells, ind, val);  // set the values of the new constraint
    }


//    showConstraints(lp);

    std::cout << "Solving linear programming ...\n";

    int ret= glp_simplex(lp_ptr, NULL);

    std:: cout << "Finished solving. \n" << std::endl;

	if (ret==0) {
		std::cout << "A solution is found. \n"<< std::endl;
		std::cout << "minimal expected loss is " << glp_get_obj_val(lp_ptr) << std::endl;
	    mechanism_matrix = DoubleMatrix(ncells, ncells);
		// fill the mechanism matrix
		for (int i=0; i<ncells; i++)
			for (int k=0; k<ncells; k++) {
				std::string s = "m"+std::to_string(i)+"_"+std::to_string(k);
				int var_ind = glp_find_col(lp_ptr, s.c_str());
				double val = glp_get_col_prim(lp_ptr, var_ind);
				mechanism_matrix.put(i,k, val);
				mechanism_file = p_mechFile;
				write_mechanism_to_file(p_mechFile);
			}
		showVarValues();
	} else
		std::cout << "The problem could not be solved!" << std::endl;

	//glp_delete_prob(lp);
}




