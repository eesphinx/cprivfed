/*
 * ScaledMechanism.h
 *
 *  Created on: Feb 8, 2020
 *      
 */

#ifndef SRC_MECHANISMS_PLAN_SCALEDMECHANISM_H_
#define SRC_MECHANISMS_PLAN_SCALEDMECHANISM_H_

/** constructing a mechanism on a fine grid from an original mechanism that operates on
 * a coarse grid. The new mechanism is defined by original mechanism info, the fine
 * grid dimensions hsize, vsize. The mechanism is virtually constructed by expanding the
 * original mechanism to the fine grid.
 */

#include "GridOutputLocalMechanism.h"

class ScaledMechanism : public GridOutputLocalMechanism {
public:

	/** constructs a mechanism on the given fine-grained grid (p_fine_hsize, p_fine_vsize) from an original mechanism
	 * (p_mech_ptr). The new fine dimensions must be integer multiples of those of the original mechanism.
	 * clone_original determines if the original mechanism is cloned or not (default is true for safety) */
	ScaledMechanism(GridOutputLocalMechanism* p_mech_ptr, double p_fine_hsize, double p_fine_vsize, bool clone_original=true);

    /** construct the mechanism from info string
     * SMECH (original_mechanism_info) S> hsize vsize  */
	ScaledMechanism(const std::string p_info);

	/** distroys the original mechanism if cloned */
	~ScaledMechanism();

    /** get the basic info string of the mechanism: has the syntax
     * SMECH (original_mechanism_info) S> hsize vsize  */
	std::string get_info() const override;
	std::string getTypeAndPriv(const int p_precision) const override;

private:


    /** returns the entry of the mechanism for input p_in and output p_out (in the linear order) */
	double getProb(int p_in, int p_out) const override;

    /** check if the original mechanism has dimensions that are multiple of this */
    void validate_mechanism_to_be_scaled(const GridOutputLocalMechanism* p_mech_ptr) const;

    /** return column (col) of the expanded version of the original mechanism matrix. col is the linear index of
     * the required column in the new (fine-grained) grid (hsize, vsize) */
    //DoubleColumn get_expanded_matrix_column(const GridOutputMechanism* p_original_mech_ptr, int p_col) const;

	GridOutputLocalMechanism* original_mechanism_ptr;

	// if the original mechanism is cloned, the destructor distroys it. By default it is set to true to be safe.
	bool original_mechanism_constructed;


};



#endif /* SRC_MECHANISMS_PLAN_SCALEDMECHANISM_H_ */
