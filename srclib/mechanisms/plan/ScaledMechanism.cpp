/*
 * ScaledMechanism.cpp
 *
 *  Created on: Feb 8, 2020
 *      Author: ehab
 */

#include <iostream>
#include <sstream>
#include "ScaledMechanism.h"

/** constructs a mechanism on the given fine-grained grid (p_fine_hsize, p_fine_vsize) from an original mechanism
 * (p_mech_ptr). The new fine dimensions must be integer multiples of those of the original mechanism.
 * clone_original determines if the original mechanism is cloned or not (default is true for safety) */
ScaledMechanism::ScaledMechanism(GridOutputLocalMechanism* p_mech_ptr, double p_fine_hsize, double p_fine_vsize,bool clone_original)
:GridOutputLocalMechanism(p_fine_hsize, p_fine_vsize), original_mechanism_constructed(clone_original) {

	validate_mechanism_to_be_scaled(p_mech_ptr);
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
}


/** construct the mechanism from info string
 * SMECH (original_mechanism_info) S> hsize vsize  */
ScaledMechanism::ScaledMechanism(const std::string p_info)
: GridOutputLocalMechanism(0,0){

	original_mechanism_constructed = true;
	std::istringstream iss(p_info);

	// these are used to parse the info string
	std::string str_token;
	char char_token;
	std::string p_original_mech_info;

	//std::cout << "parsing: " << p_info << "\n";

	if (!(iss >> str_token && str_token == AbstractOutputLocalMechanism::mechTypeToString(SMECH))) {
		std::cerr<< "unable to construct ScaledMechanism: type must be specified to "
				<< AbstractOutputLocalMechanism::mechTypeToString(SMECH);
		exit(0);
	}  // mechanism type is Ok


	// find the original mechanism info: find "(", then orig_mech_info, then find ")"
	// note we need to parse the string char by char using iss.get(char_token)
	while (iss.get(char_token) && char_token == ' '); // skip spaces until first non space char
	int lbrackets = 0; // start counting opening brackets
	if (char_token != '(') { // first non empty char
		std::cerr<< "unable to construct ScaledMechanism: found " << char_token << " instead of (, lbrackets = " << lbrackets;
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
		std::cerr<< "unable to construct ScaledMechanism: failed to fine end of the original mechanism info";
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
		std::cerr << "unable to construct ScaledMechanism: failed to construct the original mechanism";
		exit(0);
	}

	if (!(iss >> str_token && str_token == "S>")) {
		std::cerr<< "unable to construct RMechanism: failed to find the scaling operator S> ";
		exit(0);
	}  // found "S>" operator

	// find the map_file_name
	if (!(iss>> hsize >> vsize)) {
		std::cerr << "unable to construct ScaledMechanism: the fine grid is unspecified";
		exit(0);
	} // found hsize, vsize
}

ScaledMechanism::~ScaledMechanism() {
	if (original_mechanism_constructed ) {
		delete original_mechanism_ptr;
		original_mechanism_ptr = nullptr;
	}
}


void ScaledMechanism::validate_mechanism_to_be_scaled(const GridOutputLocalMechanism* p_mech_ptr) const {
	double hs = (double) hsize/p_mech_ptr->hsize;
	double vs = (double) vsize/p_mech_ptr->vsize;

	if (hs != (int) hs || vs != (int) vs) {
		std::cerr << "unable to scale mechanism:\n"
				<< "the resolution of p_fine_grid must be multiples of the mechanism grid";
		exit(0);
	}
}




/** returns the entry of the mechanism for input in and output out (in the grid linear order) */
double ScaledMechanism::getProb(int p_row, int p_col) const {
	// one cell in the coarse grid spans over a block of cells in the fine grid
	int block_width = hsize/original_mechanism_ptr->hsize;
	int block_height = vsize/original_mechanism_ptr->vsize;
	int block_size = block_width*block_height;

	//get x,y coordinates of the given row in in the fine grid.
	int row_x = p_row%hsize;
	int row_y = p_row/hsize;
	// get block coordinates in the original coarse grid
	int row_block_x = row_x/block_width;
	int row_block_y = row_y/block_height;

	//get x,y coordinates of col in the fine grid.
	int col_x = p_col%hsize;
	int col_y = p_col/hsize;

	// get block coordinates of col in the original coarse grid
	int col_block_x = col_x/block_width;
	int col_block_y = col_y/block_height;
	GridObservation obs(col_block_y * original_mechanism_ptr->hsize + col_block_x);

	return original_mechanism_ptr->getProb(row_block_y * original_mechanism_ptr->hsize + row_block_x, &obs)/block_size;
}


std::string ScaledMechanism::get_info() const {
	return AbstractOutputLocalMechanism::mechTypeToString(SMECH) + " (" +
			original_mechanism_ptr->get_info() + ") S> " + std::to_string(hsize) + " " + std::to_string(vsize);
}

std::string ScaledMechanism::getTypeAndPriv(const int p_precision) const {
	return AbstractOutputLocalMechanism::mechTypeToString(SMECH);
}




