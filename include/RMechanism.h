/*
 * RMechanism.h
 *
 *  Created on: Jan 26, 2020
 *      
 */

#ifndef MECHANISMS_PLAN_RMECHANISM_H_
#define MECHANISMS_PLAN_RMECHANISM_H_

#include "DoubleMatrix.h"
#include "GridOutputLocalMechanism.h"

/** A mechanism constructed by remapping a given GridOutputMechanism using a specific remap function,
 * or by constructing the optimal remap with respect a given prior. */

class RMechanism: public GridOutputLocalMechanism {
public:
	/** constructs a mechanism from an original mechanism (p_mech_ptr) and a remap p_remap
	 * Note we pass a pointer because GridOutputMechanism is abstract (has pure virtual functions)
	 * if 'clone_original is true (default) the original mechanism is cloned. Otherwise,
	 * the internal variable original_mechanism_ptr is set directly to the given p_mech_ptr */
	RMechanism(GridOutputLocalMechanism* p_mech_ptr, const std::string p_remap_file, bool clone_original=true);

	/** construct a remapped version of p_mech using the optimal remap computed using the prior p_prior
	 * and the loss p_loss and a cell side (required to compute the loss)
	 * the resulting remap is written to a file named remap_file */
	RMechanism(GridOutputLocalMechanism* p_mech_ptr, const Prior& p_prior,
			const Loss* p_loss_ptr, double p_cell_side, const std::string p_remap_file, bool clone_original=true);


	/** construct the mechanism from an info string
	 * info has the syntax : RMECH (original_mechanism_info) R> remap_file_name */
	RMechanism(const std::string p_info);

	/** distroy the original mechanism if it was constructed **/
	~RMechanism();


    /** returns the mechanism matrix: override the super class method to improve performance */
    DoubleMatrix get_mech_matrix() const override;

	/** get the basic info string of the mechanism:
	 * info has the syntax : RMECH (original_mechanism_info) R> remap_file_name */
	std::string get_info() const override;
	std::string getTypeAndPriv(const int p_precision) const override;


private:

    /** returns the mechanism matrix column: override the super class method to improve performance
     * p_out is the out location in linear order, i.e. out_y*hsize +out_x. */
    DoubleColumn get_mech_column(int p_out) const override;

	/** returns the entry of the mechanism for input p_in and output p_out (in the linear order) */
	double getProb(int p_in, int p_out) const override;

    std::string remap_file;

    /* the mechanism probabilities P(j|i) are all determined by the following members */
	std::vector<std::vector<int>>  remap;
	GridOutputLocalMechanism* original_mechanism_ptr;

	// if the original mechanism is cloned, the destructor distroys it. By default it is set to true to be safe.
	bool original_mechanism_constructed;


	// only for debug
    //double opt_loss=0;

};



#endif /* MECHANISMS_PLAN_RMECHANISM_H_ */
