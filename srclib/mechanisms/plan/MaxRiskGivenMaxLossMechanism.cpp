/*
 * Max_l_Risk_Given_Max_Loss_Mechanism.cpp
 *
 *  Created on: Feb 1, 2020
 *  Updated on: Aug 20, 2025
 *      Author: ehab
 */
#include <iostream>
#include <sstream>
#include "GSettings.h"
#include "Prior.h"
#include "Loss.h"
#include "LinGeoLoss.h"
#include "FullGridSpace.h"
#include "MaxRiskGivenMaxLossMechanism.h"

/** Construct the mechanism and store the matrix in a file with the name in the format 
 * <MAXADVRISK_hsize_vsize_cellside_maxQualityLoss.txt>.
*/
MaxRiskGivenMaxLossMechanism::MaxRiskGivenMaxLossMechanism(const Prior p_prior, const FullGridSpace p_xspace,
		const Loss* p_advLoss, const Loss* p_loss, double p_maxLoss)
:OptimizedMechanism(p_xspace.hsize, p_xspace.vsize) {
	cell_side = p_xspace.cell_side;
	maxUtilityLoss = p_maxLoss;
	std::string mech_type_str = AbstractOutputLocalMechanism::mechTypeToString(MAXADVRISK);
	mechanism_file = mech_type_str + "_" 
		+ std::to_string(hsize) + "_" + std::to_string(vsize) + "_" + 
		GSettings::to_string_with_precision(cell_side,2)+ "_" + 
		GSettings::to_string_with_precision(maxUtilityLoss,2) + ".txt";
		write_mechanism_to_file(mechanism_file);
	
	construct(p_prior, p_advLoss, p_loss);
	if (mechanism_matrix.sum() == 0) {
		std::cerr << "Unable to construct "<< 
		AbstractOutputLocalMechanism::mechTypeToString(MAXADVRISK) << " mechanism";
		exit(1);
	}
	if (mechanism_file != "")
		write_mechanism_to_file(mechanism_file);
//	std::cout << "adversary's l_risk = " << get_l_risk_for_prior(p_prior, p_loss, p_xspace.cell_side) << " , ";
//	std::cout << "Expected loss = " << get_exp_loss_for_prior(p_prior, p_loss, p_xspace.cell_side) << "\n";
}

/**
 * construct the MaxAdvRiskGivenMaxLoss from a given info string in the format:
 * <MAXADVRISK_hsize_vsize_cellside_maxQualityloss> 
 * contains the type, hsize, vsize, cell_side, and maxLoss. The mechanism 
 * is constructed for the uniform prior and Euclidean metric for both adversary 
 * risk and utility loss */ 
MaxRiskGivenMaxLossMechanism::MaxRiskGivenMaxLossMechanism(const std::string p_info)
:OptimizedMechanism(0,0) {
	std::string mech_type_str = AbstractOutputLocalMechanism::mechTypeToString(MAXADVRISK);
	std::istringstream iss(p_info);
	std::string p_mech_type_str; int p_hsize, p_vsize; double p_cell_side; double p_maxUtilityLoss;
	if (iss >> p_mech_type_str >> p_hsize >> p_vsize >> p_cell_side >> p_maxUtilityLoss) {
		if (p_mech_type_str != mech_type_str) {
			std::cerr << "unable to construct " << mech_type_str 
			<< " mechanism from info string:\n " << p_info << "\n";
			exit(0);
		}
		// now we parse the rest of the string .. construction work
		hsize = p_hsize; vsize = p_vsize; cell_side = p_cell_side; 
		maxUtilityLoss = p_maxUtilityLoss;
		
		double Maxdistance = std::sqrt(hsize*hsize+vsize*vsize) * cell_side;
		LinGeoLoss lossObj(Maxdistance,Maxdistance); 
		
		construct(Prior::uniform(hsize,vsize),&lossObj, &lossObj);
		if (mechanism_matrix.sum() == 0) {
			std::cerr << "Unable to construct "<< 
			AbstractOutputLocalMechanism::mechTypeToString(MAXADVRISK) << " mechanism";
			exit(1);
		}
		mechanism_file = mech_type_str + "_" 
		+ std::to_string(hsize) + "_" + std::to_string(vsize) + "_" + 
		GSettings::to_string_with_precision(cell_side,2)+ "_" + 
		GSettings::to_string_with_precision(maxUtilityLoss,2) + ".txt";
		write_mechanism_to_file(mechanism_file);
	}  // the info string starts correctly with the right mechanism type but the rest is incorrect
	else {
		std::cerr << " can not parse the info string of " << mech_type_str;
		exit(1);
	}
}


/** returns the info string in the format: <PMAT_MAXADVRISK_hsize_vsize_cellside_maxQualityloss.txt>
 * Note that this info string identifies the mechanism as PMAT mechanism instead of MAXADVRISK. This
 * is useful to avoid invoking the optimization again when the mechanism is reconstructed from the 
 * info string, for example by the estimators when they process the Empdata objects (which store the 
 * mechanisms infos together with the data).  
 */
std::string MaxRiskGivenMaxLossMechanism::get_info() const {
		return PMatrixMechanism::get_info();
		// return AbstractOutputMechanism::mechTypeToString(MAXADVRISK) + " " + std::to_string(hsize) + " "+
		// 	std::to_string(vsize) + " "+ std::to_string(cell_side) + " "+ std::to_string(maxUtilityLoss);
}

std::string MaxRiskGivenMaxLossMechanism::getTypeAndPriv(const int p_precision) const {
	return AbstractOutputLocalMechanism::mechTypeToString(MAXADVRISK) + 
	GSettings::to_string_with_precision(maxUtilityLoss, p_precision);
}


/** construct the mechanism matrix */
void MaxRiskGivenMaxLossMechanism::construct(const Prior p_prior, 
		const Loss* p_advLoss, const Loss* p_loss) {
	int ncells = hsize*vsize; // number of grid cells, also the number of columns in the linear programming problem
	FullGridSpace xspace(hsize,vsize,cell_side);

	glp_set_obj_dir(lp_ptr, GLP_MAX);   // We want to maximize the objective function
	glp_set_obj_name(lp_ptr, "l_risk");  // set the name of the objective function

	// Add some columns (the x variables) of the linear programming problem 
	// and construct the objective function (Eq. (16) in the ccs paper)
	//std::cout << "Adding the x auxiliary variables and constructing the objective function" << "\n" << std::endl;
	for (int i=0; i<ncells; i++) {
		int col_index = glp_add_cols(lp_ptr, 1);  // add one column to the problem
		std::string col_name = "x"+std::to_string(i);
		glp_set_col_name(lp_ptr, col_index, col_name.c_str());
		glp_set_col_bnds(lp_ptr, col_index, GLP_LO, 0, 0); // xj >= 0
		glp_set_obj_coef(lp_ptr, col_index, 1);
	}

	// Add the rest of problem columns (the entries mi_j of the required mechanism)
	//std::cout << "Adding the mechanism variables M_ij to the problem" << "\n" << std::endl;
	for (int i=0; i<ncells; i++)
		for (int j=0; j<ncells; j++) {
			int col_index = glp_add_cols(lp_ptr, 1); // add one column to the problem
			std::string col_name = "m"+std::to_string(i) + "_" + std::to_string(j);
			glp_set_col_name(lp_ptr, col_index, col_name.c_str());
			glp_set_col_bnds(lp_ptr, col_index, GLP_DB, 0, 1);  // set bound for the column (Eq. (20) in the CCS12 paper)
		}

	//showObjectiveFunction();

	// add privacy constrains (Eq 17) in Sh CCS12 paper
	// i : input cell, denoted in the paper by r.
	// j : an output cell, denoted in the paper by r'.
	// k : remapped (guessed )cell, denoted in the paper by r^.
	//std::cout << "Adding an (ncells+1)-length privacy constraint for every (j,k) ..." << "\n" << std::endl;
	for (int j = 0; j < ncells; ++j)
    for (int k = 0; k < ncells; ++k) {
        std::vector<int> var_indices;
        std::vector<double> var_coefs;
        
		// scan the terms (p_i loss(i,k) m_ij) in the (j,k)-th constraint
		for (int i = 0; i < ncells; ++i) {
            double prob = p_prior.getProb(i);
            double loss = p_advLoss->get_loss(
                xspace.getx(i), xspace.gety(i),
                xspace.getx(k), xspace.gety(k)
            );
            std::string col_name = "m" + std::to_string(i) + "_" + std::to_string(j);
            int col_index = glp_find_col(lp_ptr, col_name.c_str());
            if (col_index == 0) {
                std::cerr << "Error: variable " << col_name << " not found\n";
                exit(0);
            }
            var_indices.push_back(col_index);
            var_coefs.push_back(prob * loss);
        }

        // Add the (-1 x_j) term
        std::string x_name = "x" + std::to_string(j);
        int x_index = glp_find_col(lp_ptr, x_name.c_str());
        if (x_index == 0) {
            std::cerr << "Error: variable " << x_name << " not found\n";
            exit(0);
        }
        var_indices.push_back(x_index);
        var_coefs.push_back(-1.0);

        // Fill GLPK arrays. GLPK arrays are 1-based, so prepend dummy element
        std::vector<int> glpk_indices(var_indices.size() + 1);
        std::vector<double> glpk_coefs(var_coefs.size() + 1);
        for (size_t idx = 0; idx < var_indices.size(); ++idx) {
            glpk_indices[idx + 1] = var_indices[idx];
            glpk_coefs[idx + 1] = var_coefs[idx];
        }
        int row_ind = glp_add_rows(lp_ptr, 1);
        glp_set_row_bnds(lp_ptr, row_ind, GLP_LO, 0.0, 0.0); // the expression in the constraint is lower bounded by 0
        glp_set_mat_row(lp_ptr, row_ind, var_indices.size(), glpk_indices.data(), glpk_coefs.data());
    }

	// add one long (n*n)-length constraint for the utility (Eq 18) in Sh CCS12 paper.
	//std::cout << "Adding one utility constraint ..." << "\n" << std::endl;
	// Use vectors to store the indices and coefficients
	std::vector<int> constraint_var_indices;
	std::vector<double> constraint_var_coefs;

	// Fill the vectors with one entry for every (i,j)
	for (int i = 0; i < ncells; i++) {
    	for (int j = 0; j < ncells; j++) {
        	double prob = p_prior.getProb(i);
        	double loss = p_loss->get_loss(xspace.getx(i), xspace.gety(i),
                                       	xspace.getx(j), xspace.gety(j));
        	std::string col_name = "m" + std::to_string(i) + "_" + std::to_string(j);
        	int col_index = glp_find_col(lp_ptr, col_name.c_str());

        	if (col_index == 0) {
            	std::cerr << "Error: variable " << col_name << " not found\n";
            	exit(0);
        	}

        	constraint_var_indices.push_back(col_index);
        	constraint_var_coefs.push_back(prob * loss);
    	}
	}

	// GLPK arrays are 1-based, so we need to create temporary arrays starting from index 1
	std::vector<int> glpk_indices(constraint_var_indices.size() + 1);
	std::vector<double> glpk_coefs(constraint_var_coefs.size() + 1);

	for (size_t k = 0; k < constraint_var_indices.size(); k++) {
    	glpk_indices[k + 1] = constraint_var_indices[k];
    	glpk_coefs[k + 1] = constraint_var_coefs[k];
	}
	// Add a new row for the constraint
	int row_ind = glp_add_rows(lp_ptr, 1);
	glp_set_row_bnds(lp_ptr, row_ind, GLP_UP, 0, maxUtilityLoss); // upperbounded by maxUtilityLoss
	// Set the constraint matrix row
	glp_set_mat_row(lp_ptr, row_ind, constraint_var_indices.size(),
                		glpk_indices.data(), glpk_coefs.data());	
	
	
	// equality constraints (Eq. 19 in the CCS12 paper)
	//std::cout << "Adding (n)-length equality constraint for every i \n" <<std::endl;
	for (int i = 0; i < ncells; ++i) {
        std::vector<int> var_indices;
		// we dont need var_coefs because they are all 1.0
		// scan the terms (m_ij) in the constraint
		for (int j = 0; j < ncells; ++j) {
            std::string col_name = "m" + std::to_string(i) + "_" + std::to_string(j);
            int col_index = glp_find_col(lp_ptr, col_name.c_str());
            if (col_index == 0) {
                std::cerr << "Error: variable " << col_name << " not found\n";
                exit(0);
            }
            var_indices.push_back(col_index);
        }

        // Fill GLPK arrays. GLPK arrays are 1-based, so prepend dummy element
        std::vector<int> glpk_indices(var_indices.size() + 1);
        std::vector<double> glpk_coefs(var_indices.size()+ 1);
        for (size_t idx = 0; idx < var_indices.size(); ++idx) {
            glpk_indices[idx + 1] = var_indices[idx];
            glpk_coefs[idx + 1] = (1.0);
        }
        int row_ind = glp_add_rows(lp_ptr, 1);
        glp_set_row_bnds(lp_ptr, row_ind, GLP_FX, 1.0, 1.0); // the expression in the constraint is equal to 1
        glp_set_mat_row(lp_ptr, row_ind, var_indices.size(), glpk_indices.data(), glpk_coefs.data());
    }

    //showConstraints(lp);

	//std::cout << "Solving linear programming ...\n";

	glp_term_out(GLP_OFF);  // turn off all GLPK terminal output
	int ret= glp_simplex(lp_ptr, NULL);

	//std:: cout << "Finished solving. \n" << std::endl;

	if (ret==0) {
		std::cout << "A solution is found. \n";
		std::cout << "maximal expected adversary's risk " << glp_get_obj_val(lp_ptr) << std::endl;
		mechanism_matrix = DoubleMatrix(ncells, ncells);
		// fill the mechanism matrix
		for (int i=0; i < ncells; i++) {
			for (int j=0; j < ncells; j++) {
				std::string s = "m"+std::to_string(i)+"_"+std::to_string(j);
				int var_ind = glp_find_col(lp_ptr, s.c_str());
				if (var_ind ==0 ) {
					std::cout << "variable " << s <<" not found\n";
					exit(0);
				}
				double val = glp_get_col_prim(lp_ptr, var_ind);
				mechanism_matrix.put(i,j, val);
			}
		}
		//showVarValues();
	} else {
		std::cerr << "The problem could not be solved!" << std::endl;
		exit(1);
	}
}



