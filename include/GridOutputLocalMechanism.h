/*
 * GridOutputMechanism.h
 *
 *  Created on: Sep 24, 2019
 *      
 */

#ifndef SRC_MECHANISMS_GRIDOUTPUTLOCALMECHANISM_H_
#define SRC_MECHANISMS_GRIDOUTPUTLOCALMECHANISM_H_


/**
 * This is a a subclass of AbstractOutputMechanism and therefore works also on a grid of cells,
 * specified by hsize,vsize. The input and output grids have to be identical. This is also abstract
 * because sume functions pure virtual functions e.g. getTypePriv() not implemented
 */

#include <vector>
#include "GridObservation.h"

class DoubleRow;
class Prior;
class Loss;
class DoubleMatrix;
class FullGridSpace;

class GDataset;
class GeoFullGridSpace;

#include "AbstractOutputLocalMechanism.h"

class GridOutputLocalMechanism : public AbstractOutputLocalMechanism {

public:

	/** basic constructor */
	GridOutputLocalMechanism(int p_hsize, int p_vsize);

	/** returns the probability for input cell (specified in linear order) and output observation */
	double getProb(int p_in, const Observation* p_out_ptr) const;

	/** returns a mechanism column for a given observation. The given observation must be GridObservation */
	DoubleColumn get_mech_column(const Observation* p_out_ptr) const;

	/** obfuscate the given location p_loc (in the linear order) into an observation (GridObservation) */
	GridObservation* obfuscateLoc(int p_loc) const; // GridObservation is co-variant for Observation

	/** returns a pointer to newly created observation object constructed from the given string
	 * Note that the string representation of the observation depends on the mechanism, hence
	 * this function must be implemented by every subclass */
	GridObservation* readObservation(const std::string p_observation_str) const;

	/** computes the expected loss of the current mechanism for a prior and loss function.
	 * p_cell_side is passed to allow computing the loss between cells.
	 * A pointer to Loss is passed because Loss is abstract (has pure virtual functions)*/
	double get_exp_loss_for_prior(const Prior& p_prior, const Loss* p_loss_ptr, double p_cell_side) const;

	/**
	 * compute the expected losses of the given mechanism for the priors of the users in a given gdataset
	 * object (containing a list of locations) on a region. The region and cell side are specified by
	 * p_space (p_hsize, p_vsize must be equal to those of the mechanism).
	 * The function returns a list of losses for these priors. Every loss is computed only to the users having
	 * at least the given threshold of check-ins.
	 * Important: This implementation assumes that checkin rows for the same user are contiguous in the dataset
	 * */
	std::vector<double> get_exp_losses_for_gdata_priors(const GDataset& p_gdset, const GeoFullGridSpace p_space,
			const Loss* p_loss_ptr, int p_user_checkins_threshold) const;


	/**
	 * compute the expected losses of the mechanism for the the users in a given gdataset object (containing a
	 * list of locations), based on the distance from the input 'point' and the output cell. Since the loss
	 * is measured from the input checkin to the output point, we specify the geographic region coverd by the
	 * mechanism, p_space including the cell side (p_space.hsize, p_space.vsize must be equal to those of the
	 * mechanism). Every loss is computed only to the users having at least the given threshold of check-ins.
	 * Important: This implementation assumes that checkin rows for the same user are contiguous in the dataset
	 * */
	std::vector<double> get_exp_losses_for_gdata_points(const GDataset& p_gdset, const GeoFullGridSpace p_space,
			const Loss* p_loss_ptr, int p_user_checkins_threshold) const;


	/** returns the adversary's l-risk based on the prior p_prior, and adversary loss function p_loss_ptr.
	 * p_cell_side is passed to allow computing the loss between cells
	 * A pointer to Loss is passed because the Loss is abstract (has pure virtual functions) */
	double get_l_risk_for_prior(const Prior& p_prior, const Loss* p_advloss_ptr, double p_cell_side) const;

	/** returns the mechanism matrix. This is executed only if its subclass does not
	 * already have a mechanism matrix. If it has a matrix, the subclass should
	 * override this method for performance */
	virtual DoubleMatrix get_mech_matrix() const;

	/* returns a mechanism row distribution (probabilities) at p_in:
	 * it is executed only if its subclass does not
	 * already have a mechanism matrix. If the subclass does, it should
	 * override this method for performance. This method is used by EmpData class
	 * in the estimation framework to generate obfuscated data.
	 * p_in is specified in the linear order (y*hsize + x) */
	virtual DoubleRow get_mech_row(int p_in) const;

	/** Returns ideal empirical distribution on the observables (which is the same input grid):
	 *  Given a prior p_prior, it returns p_prior.A where A is the mechanism matrix */
	DoubleRow get_ideal_noisy_distribution(const Prior& p_prior) const;


	/** This method writes the current mechanism to a given file. */
	void write_mechanism_to_file(const std::string p_file) const;

	void write_diagonal_to_file(const std::string p_diag_file) const;

	/** functions that allows adapting the mechanism to finer grids and losses with respect to check-ins */


	// remap functions

	/** write the remap data to remap_file (specified as an attribute for the RMechanism) */
	static void write_remap_file(const std::vector< std::vector<int> >& p_remap, const std::string p_file);

	/** write a remap that maps every output in a grid of p_hsize * p_vsize to itself  */
	static void write_identity_remap_file(int p_hsize, int p_vsize, const std::string p_file);




protected:

	/** read the remap file into a vector of hsize*vsize vectors<int>, each vector corresponds to a target locations,
	 * and includes the mapped output locations to the target */
	std::vector< std::vector<int> > read_remap_file(const std::string p_remap_file) const;


private:
	/** returns the probability for input and output cells: must be implemented by every subclass
	 * 	It is made private because it is invoked only from the public getProb(in, Observation*) */
	virtual double getProb(int p_in, int p_out) const = 0;

	/** returns a mechanism column (out) as a column vector: executed only if its subclass does not
	 * already have a mechanism matrix. If it has a matrix, it should override this method for performance.
	 * p_out is specified in linear order (y*hsize + x).
	 * This method is invoked from the public get_mech_column(Observation* obs_ptr). This is why it shadows
	 * p_out is specified in the linear order */
	virtual DoubleColumn get_mech_column(int p_out) const;


	/**
	 * returns the conditional expected loss for the center of the input cell that incloses the given point.
	 * To determine the input cell, the region covered by the mechanism grid is specified by p_space
	 * (p_space.hsize, p_space.vsize must be equal to those of the mechanism). p_space.cell_side is also
	 * used to compute the loss from input to output.
	 */
	double get_exp_loss_for_cell(const GeoFullGridSpace p_space, double p_lat, double p_lng,
			const Loss* p_loss_ptr) const;

	/** returns the conditional expected loss for a specific input point (lat,lng) measured from the
	 * precise input point (lat, lng) to the center point of the output cell. To determine the lat/lng
	 * of the center of the output cell, we specify the region covered by the mechanism grid as p_space.
	 * (p_space.hsize, p_space.vsize must be equal to those of the mechanism).
	 * */
	double get_exp_loss_for_point(const GeoFullGridSpace& p_space, double p_lat, double p_lng,
			const Loss* p_loss_ptr) const;


};


#endif /* SRC_MECHANISMS_GRIDOUTPUTLOCALMECHANISM_H_ */
