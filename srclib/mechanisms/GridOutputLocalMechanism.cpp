/*
 * GridOutputLocalMechanism.cpp
 *
 *  Created on: Sep 24, 2019
 *      Author: ehab
 */
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include "DoubleColumn.h"
#include "GridOutputLocalMechanism.h"
#include "GridObservation.h"
#include "Distributions.h"
#include "Prior.h"
#include "Loss.h"
#include "GeoFullGridSpace.h"
#include "GDataset.h"
#include "UserStatistic.h"
#include "Distance.h"


/** basic constructor */
GridOutputLocalMechanism::GridOutputLocalMechanism(int p_hsize, int p_vsize) 
// note we have to call the direct parent constructor not the grandparent one. So parent must have it.
: AbstractOutputLocalMechanism(p_hsize,  p_vsize) {} 


/** returns a mechanism column for a given observation. The given observation must be GridObservation */
double GridOutputLocalMechanism::getProb(int p_in, const Observation* p_out_ptr) const {
	if(const GridObservation* grid_out_ptr = dynamic_cast<const GridObservation*>(p_out_ptr))
	{
		int out = grid_out_ptr->location; // safe to call
		return getProb(p_in, out);
	}
	else {
		std::cerr << "unable to cast a generic observation to grid observation to get the probability";
		exit(0);
	}
}

/** returns a mechanism column for a given observation. The given observation must be GridObservation */
DoubleColumn GridOutputLocalMechanism::get_mech_column(const Observation* p_out_ptr) const {
	if(const GridObservation* grid_out_ptr = dynamic_cast<const GridObservation*>(p_out_ptr))
	{
		int out = grid_out_ptr->location; // safe to call
		return get_mech_column(out);
	}
	else {
		std::cerr << "unable to cast a generic observation to grid observation to get a mechanism column";
		exit(0);
	}
}


/** obfuscate the given location p_loc (in the linear order) into an observation (GridObservation) */
GridObservation* GridOutputLocalMechanism::obfuscateLoc(int p_loc) const {
	DoubleRow mech_row = get_mech_row(p_loc);
	std::vector<int> obf_locs = Distributions::sampleDiscrete(mech_row, 1);
	return new GridObservation(obf_locs[0]);  // create a new observation and return a pointer to it.
}

/** returns a pointer to newly created observation object constructed from the given string
 * Note that the string representation of the observation depends on the mechanism, hence
 * this function must be implemented by every subclass */
GridObservation* GridOutputLocalMechanism::readObservation(std::string p_observation_str) const {
	return new GridObservation(std::stoi(p_observation_str));
}


/** computes the expected loss of the current mechanism for a prior and loss function.
 * p_cell_side is passed to allow computing the loss between cells.
 * A pointer to Loss is passed because Loss is abstract (has pure virtual functions)*/
double GridOutputLocalMechanism::get_exp_loss_for_prior(const Prior& p_prior, const Loss* p_loss_ptr, double p_cell_side) const {
	if (p_prior.hsize!=hsize || p_prior.vsize!=vsize) {
		std::cerr << "error: can not compute expected loss - invalid prior dimensions for the mechanism\n";
		std::cerr << "mech: hsize = " << hsize << ", vsize = " << vsize << "\n";
		std::cerr << "prior: hsize = " << p_prior.hsize << ", vsize = " << p_prior.vsize << "\n";
		exit(0);
	}

	//std::cout << "computing loss ..." << "\n";

	double exp_loss=0;
	FullGridSpace space(hsize,vsize, p_cell_side);
	for (int in = 0; in < hsize*vsize; in++) {
		double loss_in = 0;
		for (int out =0; out < hsize*vsize; out++) {
			loss_in = loss_in + getProb(in,out)*
					p_loss_ptr -> get_loss(space.getx(in), space.gety(in),
							space.getx(out), space.gety(out));
		}
		exp_loss = exp_loss + loss_in* p_prior.getProb(in);
	}
	return exp_loss;
	// a single linear algebraic expression: does not improve performance and require large memory
}



/**
 * compute the expected losses of the given mechanism for the priors of the users in a given gdataset
 * object (containing a list of locations) on a region. The region and cell side are specified by
 * p_space (p_space.hsize, p_space.vsize must be equal to those of the mechanism).
 * The function returns a list of losses for these priors. Every loss is computed only to the users having
 * at least the given threshold of check-ins.
 * Important: This implementation assumes that checkin rows for the same user are contiguous in the dataset
 * */
std::vector<double> GridOutputLocalMechanism::get_exp_losses_for_gdata_priors(const GDataset& p_gdset, const GeoFullGridSpace p_space,
		const Loss* p_loss_ptr, int p_user_checkins_threshold) const {

	// validation
	if (p_space.hsize != hsize || p_space.vsize != vsize) {
		std::cerr << "unable to evaluate losses from geographic dataset. The given space must be identical to that of the mechanism";
		exit(0);
	}

	// get the list of users having above the thresholds.
	std::vector<UserStatistic> us = p_gdset.get_users_statistics(p_user_checkins_threshold);

	std::vector<double> users_priors_exp_losses;
	// find the first check-in for the first user in the list
	int k=0;
	while (p_gdset.get_checkin(k).user_id != us[0].user_id) k++;
	std::string old_user_id = p_gdset.get_checkin(k).user_id;
	int old_user_checkins=1;
	// add the loss of the cell of the new position to the total loss of the user loss
	double old_user_prior_loss = get_exp_loss_for_cell(p_space, std::stod(p_gdset.get_checkin(k).lat),
			std::stod(p_gdset.get_checkin(k).lng), p_loss_ptr);


	for (k=k+1;k<p_gdset.get_size();k++) {
		// check if us contains the user of the current record
		size_t i=0; while (us[i].user_id != p_gdset.get_checkin(k).user_id && i < us.size()) i++;
		if (i < us.size()) {
			// check if a new user is the same as before
			if (p_gdset.get_checkin(k).user_id == old_user_id) {
				old_user_checkins++;
				//update the user's loss;
				old_user_prior_loss = old_user_prior_loss + get_exp_loss_for_cell(p_space, std::stod(p_gdset.get_checkin(k).lat),
						std::stod(p_gdset.get_checkin(k).lng),p_loss_ptr);
			}
			else {
				// here we normalize the user's loss and add the result to the list
				users_priors_exp_losses.push_back((double) old_user_prior_loss/old_user_checkins);
				// now we have a new user: update old user info
				old_user_id = p_gdset.get_checkin(k).user_id;
				old_user_checkins=1;
				// add the loss of the new position
				old_user_prior_loss = get_exp_loss_for_cell(p_space, std::stod(p_gdset.get_checkin(k).lat),
						std::stod(p_gdset.get_checkin(k).lng), p_loss_ptr);
			}
		}
	}
	// process the last user
	// here we normalize the user's loss and add the result to the list
	users_priors_exp_losses.push_back( (double) old_user_prior_loss/old_user_checkins);
	//std::cout << "processed " << users_priors_exp_losses.size() << " users" << std::endl;
	return users_priors_exp_losses;
}


/**
 * compute the expected losses of the mechanism for the the users in a given gdataset object (containing a
 * list of locations), based on the distance from the input 'point' and the output cell. Since the loss
 * is measured from the input checkin to the output point, we specify the geographic region coverd by the
 * mechanism, p_space including the cell side (p_space.hsize, p_space.vsize must be equal to those of the
 * mechanism). Every loss is computed only to the users having at least the given threshold of check-ins.
 * Important: This implementation assumes that checkin rows for the same user are contiguous in the dataset
 * */
std::vector<double> GridOutputLocalMechanism::get_exp_losses_for_gdata_points(const GDataset& p_gdset, const GeoFullGridSpace p_space,
		const Loss* p_loss_ptr, int p_user_checkins_threshold) const {
	// validation
	if (p_space.hsize != hsize || p_space.vsize != vsize) {
		std::cerr << "unable to evaluate losses from geographic dataset. The given space must be identical to that of the mechanism";
		exit(0);
	}

	// get the list of users having above the thresholds.
	std::vector<UserStatistic> us = p_gdset.get_users_statistics(p_user_checkins_threshold);

	std::vector<double> users_exp_losses;

	// find the first check-in for the first user in the list
	int k=0;
	while (p_gdset.get_checkin(k).user_id != us[0].user_id) k++;
	std::string old_user_id = p_gdset.get_checkin(k).user_id;
	int old_user_checkins=1;
	// add the loss of the cell of the new position to the total loss of the user loss
	double old_user_prior_loss = get_exp_loss_for_point(p_space, std::stod(p_gdset.get_checkin(k).lat),
			std::stod(p_gdset.get_checkin(k).lng), p_loss_ptr);


	for (k=k+1;k<p_gdset.get_size();k++) {
		// check if us contains the user of the current record
		size_t i=0; while (us[i].user_id != p_gdset.get_checkin(k).user_id && i < us.size()) i++;
		if (i < us.size()) {
			// cout << "user " << p_gdset.get_checkin(k).user_id << "is in the list" << endl;
			// check if a new user is the same as before
			if (p_gdset.get_checkin(k).user_id == old_user_id) {
				//cout << "user " << p_gdset.get_checkin(k).user_id << "is the same as before" << endl;
				old_user_checkins++;
				//update the user's loss;
				old_user_prior_loss = old_user_prior_loss + get_exp_loss_for_point(p_space, std::stod(p_gdset.get_checkin(k).lat),
						std::stod(p_gdset.get_checkin(k).lng), p_loss_ptr);
			}
			else {
				// cout << "user " << p_gdset.get_checkin(k).user_id << "is a new user" << endl;
				// process old user
				// cout << "for user " << old_user_id << " processed " << old_user_checkins << " checkins" << endl;
				// here we normalize the user's loss and add the result to the list
				users_exp_losses.push_back((double) old_user_prior_loss/old_user_checkins);
				// now we have a new user: update old user info
				old_user_id = p_gdset.get_checkin(k).user_id;
				old_user_checkins=1;
				// add the loss of the new position
				old_user_prior_loss = get_exp_loss_for_point(p_space, std::stod(p_gdset.get_checkin(k).lat),
						std::stod(p_gdset.get_checkin(k).lng), p_loss_ptr);
			}
		}
	}
	// process the last user
	// cout << "for user " << old_user_id << " processed " << old_user_checkins << " checkins" << endl;
	// here we normalize the user's loss and add the result to the list
	users_exp_losses.push_back( (double) old_user_prior_loss/old_user_checkins);
	//std::cout << "processed " << users_exp_losses.size() << " users" << std::endl;
	return users_exp_losses;
}

/** returns the adversary's l-risk based on the prior p_prior, and adversary loss function p_loss_ptr.
 * p_cell_side is passed to allow computing the loss between cells
 * A pointer to Loss is passed because the Loss is abstract (has pure virtual functions) */
double GridOutputLocalMechanism::get_l_risk_for_prior(const Prior& p_prior, const Loss* p_advloss_ptr, double p_cell_side) const {
	if (p_prior.hsize!=hsize || p_prior.vsize!=vsize) {
		std::cerr << "unable to compute l-risk: the prior dimensions are incompatible with the mechanism\n";
		exit(0);
	}

	DoubleColumn prior_as_column = p_prior.get_prior_vector().transpose();   // a column vector
	DoubleMatrix advloss_mat = p_advloss_ptr->get_loss_matrix(hsize, vsize, p_cell_side);

	// we use: l-risk = sum_y min_w sum_x p(x).p(y|x).loss(w,x)
	double advloss=0;
	for (int out=0;out<hsize*vsize;out++)
		advloss = advloss + advloss_mat.mmul(prior_as_column.mul(get_mech_column(out))).min();

	return advloss;
}

/** returns the mechanism matrix. This is executed only if its subclass does not
 * already have a mechanism matrix. If it has a matrix, the subclass should
 * override this method for performance */
DoubleMatrix GridOutputLocalMechanism::get_mech_matrix() const {
	DoubleMatrix retMat(hsize*vsize, hsize*vsize);
	for (int in = 0; in < hsize*vsize; in ++)
		for (int out =0; out < hsize*vsize; out++) {
			retMat.put(in, out,getProb(in,out));
		}
	return retMat;
}

/* returns a mechanism row distribution (probabilities) at p_in:
 * it is executed only if its subclass does not
 * already have a mechanism matrix. If the subclass does, it should
 * override this method for performance. This method is used by EmpData class
 * in the estimation framework to generate obfuscated data.
 * p_in is specified in linear order (y*hsize + x) */
DoubleRow GridOutputLocalMechanism::get_mech_row(int p_in) const {
	if ((p_in < 0) || (p_in>=hsize*vsize)) {
		std::cerr << "unable to return mechanism row: invalid row index";
		exit(0);
	}
	DoubleRow retRow(hsize*vsize);
	for (int out=0; out<hsize*vsize;out++)
		retRow.put(out, getProb(p_in, out));
	return retRow;
}

/** Returns ideal empirical distribution on the observables (which is the same input grid):
 *  Given a prior p_prior, it returns p_prior.A where A is the mechanism matrix */
DoubleRow GridOutputLocalMechanism::get_ideal_noisy_distribution(const Prior& p_prior) const {
	if (p_prior.hsize!=hsize || p_prior.vsize!=vsize) {
		std::cerr << "error: can not compute expected loss - invalid prior dimensions for the mechanism\n";
		std::cerr << "mech: hsize = " << hsize << ", vsize = " << vsize << "\n";
		std::cerr << "prior: hsize = " << p_prior.hsize << ", vsize = " << p_prior.vsize << "\n";
		exit(0);
	}

	return p_prior.get_prior_vector().mmul(get_mech_matrix());
}



/** This method writes the current mechanism to a given file. */
void GridOutputLocalMechanism::write_mechanism_to_file(std::string p_file) const {
	std::ofstream pfile (p_file);
	if (pfile.is_open()) {
		pfile << "#info: " << get_info() << "\n\n";
		for (int in = 0; in < hsize*vsize; in++) {
			for (int out =0; out < hsize*vsize; out++)
				pfile << getProb(in, out) << " ";
			pfile << "\n";
		}
		pfile.close();
	}
	else {
		std::cerr << "Unable to open file to save the mechanism";
		exit(0);
	}

}

void GridOutputLocalMechanism::write_diagonal_to_file(std::string p_diag_file) const {
	std::ofstream pfile (p_diag_file);
	if (pfile.is_open()) {
		for (int loc = 0; loc < hsize*vsize; loc++)
			pfile << getProb(loc, loc) << " ";
		pfile.close();
	}
	else {
		std::cerr << "Unable to open file to write the mechanism diagonal";
		exit(0);
	}
}

/** write the remap data to remap_file (specified as an attribute for the RMechanism) */
void GridOutputLocalMechanism::write_remap_file(const std::vector<std::vector<int>>& p_remap, const std::string p_file) {
	std::ofstream pfile (p_file);
	if (pfile.is_open()) {
		// write data
		for (size_t t=0; t< p_remap.size() ; t++)
			if (p_remap[t].size()>0) {
				pfile << t << " < ";
				for (size_t out=0; out < p_remap[t].size(); out++)
					pfile << p_remap[t][out] << " ";
				pfile << "\n";
			}
		pfile.close();
	}
	else {
		std::cerr << "Unable to open file " << p_file << "to write the remap data ";
		exit(0);
	}
}

/** write a remap that maps every output in a grid of p_hsize * p_vsize to itself  */
void GridOutputLocalMechanism::write_identity_remap_file(int p_hsize, int p_vsize, const std::string p_file) {
	std::ofstream pfile (p_file);
	if (pfile.is_open()) {
		// write data
		for (int t=0; t< p_hsize*p_vsize ; t++)
			pfile << t << " < " << t << "\n";
		pfile.close();
	}
	else {
		std::cerr << "Unable to open file " << p_file << "to write the remap data ";
		exit(0);
	}

}

/** read the remap file into a vector of hsize*vsize vectors<int>, each vector corresponds to a target locations,
 * and includes the mapped output locations to the target */
std::vector<std::vector<int>> GridOutputLocalMechanism::read_remap_file(const std::string p_remap_file) const {
	std::ifstream pfile(p_remap_file);
	if (!pfile.is_open()) {
		std::cerr << "unable to open the remap file: " << p_remap_file;
		exit(0);
	}
	// create hsize*vsize empty vectors, each one corresponds to a target and will contain the outputs mapped to the target
	std::vector<std::vector<int>> remap(hsize*vsize);
	std::string line;
	int target;
	int output;
	std::string map_op;

	while (getline(pfile,line)) {
		if (line.length() == 0 || line.front()=='#')     // skip empty or commented line
			continue;
		// if the above conditions are not true proceed to read the current line
		std::stringstream ss(line);
		if (!(ss >> target && target < hsize*vsize)) {
			std::cerr << "unable to parse the remap file: " << p_remap_file << ": line must start with " << target << " which is not a valid location";
			exit(0);
		}
		if (!(ss >> map_op && map_op == "<")) {
			std::cerr << "unable to parse the remap file: " << p_remap_file << ": map operator is missing";
			exit(0);
		}
		while (ss >> output) {
			if (output >= hsize*vsize) {
				std::cerr << "unable to parse the remap file: " << p_remap_file << ": line contains " << output << " which is not a valid location";
				exit(0);
			}
			remap[target].push_back(output);
		}
	}
	pfile.close();
	return remap;
}

/** returns a mechanism column (out) as a column vector: executed only if its subclass does not
 * already have a mechanism matrix. If it has a matrix, it should override this method for performance.
 * p_out is specified in linear order (y*hsize + x) */
DoubleColumn GridOutputLocalMechanism::get_mech_column(int p_out) const {
	if ((p_out < 0) || (p_out>=hsize*vsize)) {
		std::cerr << "unable to return mechanism column: invalid column index";
		exit(0);
	}
	DoubleColumn retColVec(hsize*vsize);
	for (int in=0; in<hsize*vsize;in++)
		retColVec.put(in,getProb(in,p_out));
	return retColVec;

}

/**
 * returns the conditional expected loss for the center of the input cell that incloses the given point.
 * To determine the input cell, the region covered by the mechanism grid is specified by p_space
 * (p_space.hsize, p_space.vsize must be equal to those of the mechanism). p_space.cell_side is also
 * used to compute the loss from input to output.
 */
double GridOutputLocalMechanism::get_exp_loss_for_cell(const GeoFullGridSpace p_space, double p_lat, double p_lng,
		const Loss* p_loss_ptr) const {
	// validation
	if (p_space.hsize != hsize || p_space.vsize != vsize) {
		std::cerr << "unable to evaluate losses from geographic dataset. The given space must be identical to that of the mechanism";
		exit(0);
	}
	// get discrete (cell) coordinates relative to the space
	int in_x = p_space.get_cell_x_of_lng(p_lng);
	int in_y = p_space.get_cell_y_of_lat(p_lat);

	double loss =0;
	// consider all (out) cells of the mechanism grid
	for (int out=0; out<hsize*vsize; out++) {
		int out_x = out%hsize; int out_y = out/hsize;
		double in_to_out_dist = Distance::get_euclid_dist(in_x*p_space.cell_side, in_y*p_space.cell_side,
				out_x*p_space.cell_side, out_y*p_space.cell_side);
		loss = loss + getProb(in_y*hsize+in_x , out) *
				p_loss_ptr->get_loss_by_dist(in_to_out_dist);
	}
	return loss;
}

/** returns the conditional expected loss for a specific input point (lat,lng) measured from the
 * precise input point (lat, lng) to the center point of the output cell. To determine the lat/lng
 * of the center of the output cell, we specify the region covered by the mechanism grid as p_space.
 * (p_space.hsize, p_space.vsize must be equal to those of the mechanism).
 * */
double GridOutputLocalMechanism::get_exp_loss_for_point(const GeoFullGridSpace& p_space, double p_lat, double p_lng,
		const Loss* p_loss_ptr) const {
	// validation
	if (p_space.hsize != hsize || p_space.vsize != vsize) {
		std::cerr << "unable to evaluate losses from geographic dataset. The given space must be identical to that of the mechanism";
		exit(0);
	}

	// get discrete (cell) coordinates relative to the space (used to obtain the mechanism probabilities)
	int in_x = p_space.get_cell_x_of_lng(p_lng);
	int in_y = p_space.get_cell_y_of_lat(p_lat);

	double loss =0;
	for (int out=0; out < hsize*vsize; out++) {
		int out_x = out%hsize;
		int out_y= out/hsize;
		// dx, dy are the continuous coordinates of out from south-west corner
		double dy = (0.5+ out_y)*p_space.cell_side;
		double dx = (0.5+ out_x)*p_space.cell_side;
		// get the lat/lng of out
		double out_cell_center_lat = Distance::adjust_latitude(p_space.south, p_space.west, dy);
		double out_cell_center_lng = Distance::adjust_longitude(p_space.south, p_space.west, dx);
		// measure the distance from in to out
		double in_to_out_dist = Distance::get_distance_on_earth(p_lat, p_lng,
				out_cell_center_lat, out_cell_center_lng);
		loss = loss + getProb(in_y*hsize+in_x, out) *
				p_loss_ptr->get_loss_by_dist(in_to_out_dist);
	}
	return loss;
}



