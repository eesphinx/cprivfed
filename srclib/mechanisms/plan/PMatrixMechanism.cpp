#include <string>
#include <sstream>    // for reading from string as a file
#include <fstream>    // for reading and writing to files
#include <iostream>   // for cerr, cout
#include "DoubleRow.h"
#include "DoubleColumn.h"
#include "FullGridSpace.h"
#include "PMatrixMechanism.h"

PMatrixMechanism::PMatrixMechanism(int p_hsize, int p_vsize)
:GridOutputLocalMechanism(p_hsize,p_vsize) {
}

/** constructs a mechanism from a given file, in consistency with the given space. That is
 * More precisely, the dimensions of the constructed mechanisms have p_space hsize, vsize.
 *  Any line starting with # is skipped **/
PMatrixMechanism::PMatrixMechanism(int p_hsize, int p_vsize, std::string p_file)
:GridOutputLocalMechanism(p_hsize, p_vsize), mechanism_file(p_file), mechanism_matrix(p_hsize*p_vsize,p_hsize*p_vsize) {
	read_file(p_file);
}


PMatrixMechanism::PMatrixMechanism(const FullGridSpace& p_space, std::string p_file)
: GridOutputLocalMechanism(p_space.hsize, p_space.vsize),
  mechanism_matrix(p_space.hsize*p_space.vsize,p_space.hsize*p_space.vsize)  {
	read_file(p_file);
}


PMatrixMechanism::PMatrixMechanism(const std::string p_info)
: GridOutputLocalMechanism(0,0){
	std::string mech_type_str = AbstractOutputLocalMechanism::mechTypeToString(PMAT);
	std::istringstream iss(p_info);
	std::string p_mech_type_str;
	int p_hsize, p_vsize;
	std::string p_mech_file;
	if (iss >> p_mech_type_str >> p_hsize >> p_vsize >> p_mech_file) {
		if (p_mech_type_str != mech_type_str) {
			std::cerr << "unable to construct " << mech_type_str << " mechanism from info string:\n "
					<< "info string must start with " << mech_type_str;
			exit(0);
		}
		// construction work
		hsize = p_hsize; vsize = p_vsize;
		mechanism_matrix = DoubleMatrix(p_hsize*p_vsize, p_hsize*p_vsize);
		//matrix_file(p_mech_input_file)
		read_file(p_mech_file);
		mechanism_file = p_mech_file;
		//std::cout<< "constructed the mechansim using info: " << p_info << "\n";
	}
	else {
		std::cerr << "unable to parse the mechanism info to construct "<< mech_type_str <<" mechanism";
		exit(0);
	}

}


/** returns the mechanism matrix: override the super class method to improve performance */
DoubleMatrix PMatrixMechanism::get_mech_matrix() const {
	return mechanism_matrix;
}

/** returns a mechanism matrix row probabilities: override the super class method to improve performance
 * used For EM. p_in is the out location in linear order, i.e. out_y*hsize +out_x. */
DoubleRow PMatrixMechanism::get_mech_row(int p_in) const {
	return mechanism_matrix.selectRow(p_in);
}

/** returns the mechanism matrix column: override the super class method to improve performance
 * p_out is the out location in linear order, i.e. out_y*hsize +out_x. */
DoubleColumn PMatrixMechanism::get_mech_column(int p_out) const {
	return mechanism_matrix.selectColumn(p_out);
}

/** returns the entry of the mechanism for input in and output out (specified in the linear order) */
double PMatrixMechanism::getProb(int p_in,int p_out) const {
	return mechanism_matrix.get(p_in, p_out);
}


void PMatrixMechanism::read_file(std::string p_file) {
	std::string line;         // to hold a line string from the input file
	std::ifstream pfile(p_file);
	if (!(pfile.is_open())) {
		std::cerr << "mechanism input file " << p_file << " is not found";
		exit(0);
	}

	int row=0;
	double prob;
	while (getline(pfile,line)) {
		if (line.length() == 0)    // skip empty line
			continue;
		if (line.substr(0,1).compare("#")==0) // skip the commented line
			continue;

		// otherwise
		int column =0;
		std::stringstream ss(line);
		while (ss >> prob) {
			mechanism_matrix.put(row, column, prob);
			column++;
		}    // here column should be the number of seen columns
		if (column != hsize*vsize) {
			std::cerr << "invalid mechanism:too few or too many probabilities in one row";
			exit(0);
		}
		row++;
	}   // here row should be the number of encountered valid rows

	if (row != hsize*vsize) {
		std::cerr << "invalid mechanism:too few/many probability rows";
		exit(0);
	}
	pfile.close();
}

std::string PMatrixMechanism::get_info() const {
	return AbstractOutputLocalMechanism::mechTypeToString(PMAT) + " " +
			std::to_string(hsize) + " "+ std::to_string(vsize) + " " + mechanism_file;
}

std::string PMatrixMechanism::getTypeAndPriv(const int p_precision) const {
	return AbstractOutputLocalMechanism::mechTypeToString(PMAT)+"_"+mechanism_file;
}

