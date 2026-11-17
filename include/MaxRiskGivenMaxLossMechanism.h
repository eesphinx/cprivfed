/*
 * MaxRiskGivenMaxLossMechanism.h
 *
 *  Created on: Feb 1, 2020
 *      
 */

#ifndef SRC_MECHANISMS_PLAN_MAXRISKGIVENMAXLOSSMECHANISM_H_
#define SRC_MECHANISMS_PLAN_MAXRISKGIVENMAXLOSSMECHANISM_H_

#include "OptimizedMechanism.h"
#include "Loss.h"

class MaxRiskGivenMaxLossMechanism: public OptimizedMechanism {
public:

	/** Constructs the optimal mechanisms that maximizes the adversary's l-risk (based on generic aversary loss function
	 * p_advLoss) under the constraint that the user's quality loss (with respect to a quality loss function p_loss) is
	 * at most p_maxLoss (to satisfy a reasonable quality loss). p_mechFile is used to store the mechanism
	 * [Shokri et al, CCS2012].
	 * Important: The user quality loss does not take into account the optimal remap that the user may apply */
	MaxRiskGivenMaxLossMechanism(const Prior p_prior, const FullGridSpace p_xspace,
				const Loss* p_advLoss, const Loss* p_loss, double p_maxLoss);

	/** Constructs the mechanism from the given info string 
	 * It uses the uniform prior and the Euclidean distance 
	 * for both the adversary and utility losses.  
	*/
	MaxRiskGivenMaxLossMechanism(const std::string p_info);

	/** returns the info string in the format: <PMAT_MAXADVRISK_hsize_vsize_cellside_maxQualityloss.txt>
 	* Note that this info string identifies the mechanism as PMAT mechanism instead of MAXADVRISK. This
 	* is useful to avoid invoking the optimization again when the mechanism is reconstructed from the 
 	* info string, for example by the estimators when they process the Empdata objects (which store the 
 	* mechanisms infos together with the data).  
 	*/
	std::string get_info() const override; 
	std::string getTypeAndPriv(const int p_precision) const override;
	
private:
 	
	// we store these information to allow constructing the info string in get_info()
    double cell_side;          // this is important to interpret the loss values
	double maxUtilityLoss;     // maximum allowed utility loss value

	/* returns the mechanism matrix */
	void construct(const Prior p_prior, const Loss* p_advLoss, 
		const Loss* p_loss);
};


#endif /* SRC_MECHANISMS_PLAN_MAXRISKGIVENMAXLOSSMECHANISM_H_ */
