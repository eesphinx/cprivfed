/*
 * PMatrixMechanism.h
 *
 *  Created on: Jul 7, 2019
 *      
 */

#ifndef PRIV_PLAN_MECHANISMS_PMATRIXMECHANISM_H_
#define PRIV_PLAN_MECHANISMS_PMATRIXMECHANISM_H_

/**
 * A mechanism (on a grid) read from a matrix stored in a file. The grid parameters (hsize, vsize).
 * are passed in the constructor. rows/columns in the file are in S order.
 * Optimized mechanisms are subclasses of this mechanism (all of them are of type PMAT).
 * Useful to interface with Kostas qif library.
 * */

#include <string>
#include "GridOutputLocalMechanism.h"
#include "DoubleMatrix.h"

class PMatrixMechanism: public GridOutputLocalMechanism {
public:

	PMatrixMechanism(int p_hsize, int p_vsize);

	/** constructs a mechanism from a given file, in consistency with the given space. That is
     * More precisely, the dimensions of the constructed mechanisms have p_space hsize, vsize.
     *  Any line starting with # is skipped **/
	PMatrixMechanism(int p_hsize, int p_vsize, std::string p_mech_input_file);


	PMatrixMechanism(const FullGridSpace& p_space, std::string p_mech_input_file);

	/** construct the mechanism from info string */
	PMatrixMechanism(const std::string p_info);

    /** returns the mechanism matrix: override the super class method to improve performance */
    DoubleMatrix get_mech_matrix() const override;

    /** returns a mechanism matrix row probabilities: override the super class method to improve performance
     * used For EM. p_in is the out location in linear order, i.e. out_y*hsize +out_x. */
    DoubleRow get_mech_row(int p_in) const override;

	/** get the basic info string of the mechanism */
	std::string get_info() const override;
	std::string getTypeAndPriv(const int p_precision) const override;

protected:

	// important to store this mechanism matrix because mechanism info (in empdata) will refer to the file
    std::string mechanism_file;

    // mechanism_matrix
    DoubleMatrix mechanism_matrix;

private:


	/** returns the entry of the mechanism for input in and output out (in linear order). */
	double getProb(int p_in, int p_out) const override;

    /** returns the mechanism matrix column: override the super class method to improve performance
     * p_out is the out location in linear order, i.e. out_y*hsize +out_x. */
    DoubleColumn get_mech_column(int p_out) const override;

    void read_file(std::string p_file);
};




#endif /* PRIV_PLAN_MECHANISMS_PMATRIXMECHANISM_H_ */
