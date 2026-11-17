/*
 * RMechanism.cpp
 *
 *  Created on: Jan 26, 2020
 *      Author: ehab
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include "DoubleColumn.h"
#include "FullGridSpace.h"
#include "GridObservation.h"
#include "Prior.h"
#include "Loss.h"
#include "RMechanism.h"


/** constructs a mechanism by remapping an original mechanism (p_mech_ptr) using a remap p_remap.
 * Important: in the remap file every output has to be mapped to a target even itself
 * Note we pass a pointer because GridOutputMechanism is abstract (has pure virtual functions)*/
RMechanism::RMechanism(GridOutputLocalMechanism* p_mech_ptr, const std::string p_remap_file, bool clone_original)
: GridOutputLocalMechanism(p_mech_ptr->hsize, p_mech_ptr->vsize),
  remap_file(p_remap_file), original_mechanism_constructed(clone_original)  {

	if (original_mechanism_constructed)
		if (GridOutputLocalMechanism* ptr =
				dynamic_cast<GridOutputLocalMechanism*>(AbstractOutputLocalMechanism::createMechanism(p_mech_ptr->get_info()))) {
			original_mechanism_ptr = ptr;
			ptr = nullptr;
		}
		else {
			std::cerr << "unable to construct RMechanism: unable to clone the original mechanism";
			exit(0);
		}
	else
		original_mechanism_ptr = p_mech_ptr;

	remap = read_remap_file(remap_file);
}


/** construct a remapped version of p_mech using the optimal remap computed using the prior p_prior and
 * an abstract loss p_loss_ptr. p_cell_side is also passed because it is used to compute and therefore
 * to optimize the loss */
RMechanism::RMechanism(GridOutputLocalMechanism* p_mech_ptr, const Prior& p_prior,
		const Loss* p_loss_ptr, double p_cell_side, const std::string p_remap_file, bool clone_original)
: GridOutputLocalMechanism(p_mech_ptr->hsize, p_mech_ptr->vsize),
  remap_file(p_remap_file), original_mechanism_constructed(clone_original) {

	if (hsize != p_prior.hsize || vsize != p_prior.vsize) {
		std::cerr << "the mechanism and the prior are incompatible";
		exit(0);
	}

	if (original_mechanism_constructed)
		if (GridOutputLocalMechanism* ptr =
				dynamic_cast<GridOutputLocalMechanism*>(AbstractOutputLocalMechanism::createMechanism(p_mech_ptr->get_info()))) {
			original_mechanism_ptr = ptr;
			ptr = nullptr;
		}
		else {
			std::cerr << "unable to construct RMechanism: unable to clone the original mechanism";
			exit(0);
		}
	else
		original_mechanism_ptr = p_mech_ptr;


	//std::cout << "Loss before remapping = " << p_mech_ptr->get_exp_loss_for_prior(xspace,p_prior, p_loss_ptr) << "\n";

	DoubleColumn prior_as_column = p_prior.get_prior_vector().transpose();   // a column vector
	DoubleMatrix loss_mat = p_loss_ptr->get_loss_matrix(hsize, vsize, p_cell_side);

	// each element of remap corresponds to a target and is a vector containing the outputs mapped to this target
	remap.resize(hsize*vsize);

	// construct the optimal remap: inspecting every output to determine which target it should be remapped to.
	for (int out=0;out<hsize*vsize;out++) {

		// notes: In general we have x--perturbed--> y and we want to map y --> z
		// using the formula z(y) = argmin_z Sum_x P(x|y).loss(x,z)
		//                        = argmin_z Sum_x P(x).P(y|x).loss(x,z)
		// The above formula is computed by element-wise multiplying column vectors P(x)
		// and P(y|x) to yield column vector C of dim(X*1) and then matrix-multiplying
		// the loss_matrix L of dim(Z*X) and the vector C.
		// This yield a column vector (losses_of_targets) with dim(Z*1) in which the minimum
		// entry is the minimum loss, and the corresponding z is the best target for
		// y=out. Note in our case Z (guesses) is exactly the domain X, however if they are different
		// the above loss_matrix has to be in the form Z*X, i.e. |Z| rows and |X| columns.
		// Note also the minimum entry in losses_of_targets is exactly the conditional
		// l-risk (the dual of g-vulnerability) with l=loss when the observable out is y.
		// performing this operation as a linear algebraic expression instead of using
		// loops has significantly improved the computation.
		// Note the overall loss of the remapped mechanism = Sum_y (Sum_x P(x).P(y|x).loss(x, z(y))).
		// The inner sum in the above expression is exactly the minimum value in losses_of_targets
		// of output y. Therefore the overall expected loss is obtained by summing these values over y.

		//System.out.println("processing out = "+ out);
		// target_losses = loss_mat * (prior_vect % out_col)
		// out_col is the out'th column in the mechanism matrix.
		// * is matrix product, % is element-wise product

		GridObservation out_obs(out);
		DoubleColumn out_col = original_mechanism_ptr->get_mech_column(&out_obs);
		DoubleColumn losses_of_targets = loss_mat.mmul(prior_as_column.mul(out_col));
		int best_target = losses_of_targets.argmin();
		remap[best_target].push_back(out);

		//opt_loss=opt_loss + losses_of_targets.get(best_target);
		//remapped_mechanism_mat.putColumn(best_target, remapped_mechanism_mat.selectColumn(best_target).addi(out_col)); // col(t) = col(t) + col(out)
	}

	//std::cout << "Loss after remapping = " << opt_loss << "\n";
	write_remap_file(remap, remap_file);
}



/** construct the mechanism from an info string
 * info has the syntax : RMECH ( original_mechanism_info ) R> remap_file_name */
RMechanism::RMechanism(const std::string p_info)
: GridOutputLocalMechanism(0,0)
{
	original_mechanism_constructed = true;
	std::istringstream iss(p_info);

	// these are used to parse the info string
	std::string str_token;
	char char_token;
	std::string p_original_mech_info;

	//std::cout << "parsing: " << p_info << "\n";

	if (!(iss >> str_token && str_token == AbstractOutputLocalMechanism::mechTypeToString(RMECH))) {
		std::cerr<< "unable to construct RMechanism: type must be specified to "
				<< AbstractOutputLocalMechanism::mechTypeToString(RMECH);
		exit(0);
	}  // mechanism type is Ok


	// find the original mechanism info: find "(", then orig_mech_info, then find ")"
	// note we need to parse the string char by char using iss.get(char_token)
	while (iss.get(char_token) && char_token == ' '); // skip spaces until first non space char
	int lbrackets = 0; // start counting opening brackets
	if (char_token != '(') { // first non empty char
		std::cerr<< "unable to construct RMechanism: found " << char_token << " instead of (, lbrackets = " << lbrackets;
		exit(0);
	}  // found '('
	lbrackets++;
	while (lbrackets>0 && iss.get(char_token)) {
		if (char_token == '(')
			lbrackets++;
		else if (char_token == ')')
			lbrackets--;
		if (lbrackets>0)
			//add the char to mech_info (including brackets) as long as lbrackets>0
			p_original_mech_info = p_original_mech_info + char_token;
	}
	if (lbrackets!=0) { // i.e. the parsed string ended while there are opening brackets
		std::cerr<< "unable to construct RMechanism: failed to fine end of the original mechanism info";
		exit(0);
	}  // found ")" and the original mechanism string obtained
	//std::cout << "original_mechanism info = " << p_original_mech_info << "\n";

	// Construct the original mechanism.
	if (GridOutputLocalMechanism* ptr =
			dynamic_cast<GridOutputLocalMechanism*>(AbstractOutputLocalMechanism::createMechanism(p_original_mech_info))) {
		original_mechanism_ptr = ptr;
		original_mechanism_constructed = true;
		ptr = nullptr;
	}
	else {
		std::cerr << "unable to construct RMechanism: failed to construct the original mechanism";
		exit(0);
	}

	if (!(iss >> str_token && str_token == "R>")) {
		std::cerr<< "unable to construct RMechanism: failed to find the map operator R> ";
		exit(0);
	}  // found "R>" operator

	// find the map_file_name
	if (!(iss>> str_token)) {
		std::cerr << "unable to construct RMechanism: the remap file is unspecified";
		exit(0);
	} // found remap file name

	remap_file = str_token;
	hsize = original_mechanism_ptr->hsize; vsize = original_mechanism_ptr->vsize;
	remap = read_remap_file(remap_file);
}

RMechanism::~RMechanism() {
	if (original_mechanism_constructed ) {
		delete original_mechanism_ptr;
		original_mechanism_ptr = nullptr;
	}
}


/** returns the mechanism matrix: override the super class method to improve performance */
DoubleMatrix RMechanism::get_mech_matrix() const {
	DoubleMatrix mat(hsize*vsize, hsize*vsize);
	for (int t=0; t< hsize*vsize; t++) {
		mat.putColumn(t, get_mech_column(t));
	}
	return mat;
}

/** returns the mechanism matrix column: override the super class method to improve performance
 * p_out is the out location in linear order, i.e. out_y*hsize +out_x. */
DoubleColumn RMechanism::get_mech_column(int p_out) const {
	DoubleColumn col(hsize*vsize);
	col.fillWith(0);
	for (size_t j=0; j<remap[p_out].size(); j++) {
		GridObservation obs(remap[p_out][j]);
		col.addi(original_mechanism_ptr->get_mech_column(&obs));
	}
	return col;
}

/** returns the entry of the mechanism for input in and output out (in the grid linear order) */
double RMechanism::getProb(int p_in, int p_out) const {
	double prop =0;
	for (size_t j=0; j<remap[p_out].size(); j++) {
		GridObservation obs(remap[p_out][j]);
		prop = prop + original_mechanism_ptr->getProb(p_in, &obs);
	}
	return prop;
}

std::string RMechanism::get_info() const {
	return AbstractOutputLocalMechanism::mechTypeToString(RMECH) + " (" +
			original_mechanism_ptr->get_info() + ") R> " + remap_file;
}

std::string RMechanism::getTypeAndPriv(const int p_precision) const {
	return AbstractOutputLocalMechanism::mechTypeToString(RMECH);
}


