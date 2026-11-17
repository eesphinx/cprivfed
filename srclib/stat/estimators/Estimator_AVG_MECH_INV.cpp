#include <iostream>
#include "DoubleRow.h"
#include "DoubleColumn.h"
#include "GridOutputLocalMechanism.h"
#include "FullGridSpace.h"
#include "EmpData.h"
#include "Distributions.h"
#include "Prior.h"
#include "Estimator_AVG_MECH_INV.h"




/**
 * Construct an estimator for the distribution on a given space of secrets,
 * and a set of empirical data objects resulting from obfuscation using various mechanisms. We construct
 * from the EmpData object the average (full) (noisy) empirical distribution which may have zero components
 * (for the cells that are not observed). P.S. This distribution is obtained by invoking a function in the
 * xspace (for every empdata object) to encounter all locations.
 * An average mechanism is obtained from the individual mechanisms of the given empirical data.
 */

/** constructs an average mechanism from the individual mechanisms in p_empdata and run INV on the average mechanism */
Estimator_AVG_MECH_INV::Estimator_AVG_MECH_INV (const FullGridSpace p_xspace, const std::vector<EmpData>& p_list_empdata) {
	
	// combine empdata objects that have the same mechanism 
	std::vector<EmpData> list_empdata = EmpData::compress(p_list_empdata);
	
	// if all mechanisms are kRR construct an estimator optimized for the kRR
	std::string krr_str = Mechanism::mechTypeToString(KRR);
    if (std::all_of(list_empdata.begin(), list_empdata.end(), 
        [krr_str](EmpData x){ return (x.getObfMechInfo().compare(0, krr_str.length(), krr_str) == 0);})) {
		all_mechanisms_are_krr = true;
		avg_mech_inv_krr = Estimator_AVG_MECH_INV_KRR(p_xspace,list_empdata);
	}
	else {
		DoubleMatrix pr_z_giv_x(p_xspace.getSize(), p_xspace.getSize());
		pr_z_giv_x.fillWith(0);
		DoubleRow emp_pr_z(p_xspace.getSize());
		emp_pr_z.fillWith(0);

		// obtain the total observations from the given data objects
		int total_observations = 0;
		for (size_t j=0; j< list_empdata.size(); j++)
		total_observations = total_observations + list_empdata[j].get_totalObservations();

		//todo: debug
		//std::cout << "pr_z_giv_x = \n" << pr_z_giv_x.toString() << "\n";
		//std::cout << "emp_pr_z = " << emp_pr_z.toString() << "\n\n";

		for (size_t j=0; j< list_empdata.size(); j++) {
			Mechanism* mech_ptr = Mechanism::createMechanism(list_empdata[j].getObfMechInfo());
			GridOutputLocalMechanism* local_mech_ptr = dynamic_cast<GridOutputLocalMechanism*>(mech_ptr);
			if (!local_mech_ptr) {  // the given mechanism is not local
				std::cerr << "unable to construct AVG_MECH_IBU: the mechanism of data object " << j << " is not local grid output mechanism";
				delete(mech_ptr); mech_ptr=nullptr;
				exit(1);
			}
			if (local_mech_ptr->hsize != p_xspace.hsize || local_mech_ptr->vsize != p_xspace.vsize) {
				std::cerr << "unable to construct an AVG_MECH_INV estimator: "
					<< "The grid of mechanism ("<< local_mech_ptr->getTypeStr()<<")"
					<< "that produced the empirical data object " << j << "is not compatible with the given space";
				delete(mech_ptr); mech_ptr=nullptr;
				exit(1);
			}
			// Now everything should be fine

			//todo: debug
			//std::cout << "processing noisy data: " << j<< "\n";
			double ratio = (double) list_empdata[j].get_totalObservations()/total_observations;
			//std::cout << "ratio = " << std::to_string(ratio) <<"\n";
			//std::cout << "mechanism matrix (" <<j<<")= " << grid_mech_ptr-> get_mech_matrix().selectRow(0).toString() << "\n";
			pr_z_giv_x.addi( (local_mech_ptr-> get_mech_matrix()).mul(ratio)  );

			//std::cout << "after seeing data ("<<j<<"): pr_z_giv_x = \n" << pr_z_giv_x.selectRow(0).toString() << "\n";

			//todo: debug
			//std::cout << "emp_pr_z("<<j<<") = " << p_xspace.get_distribution_from_emp_data(p_empdata[j]).toString() << "\n";
			//Prior(p_xspace.get_distribution_from_emp_data(p_empdata[j]), p_xspace.hsize, p_xspace.vsize)
			//.write_prior_as_matrix("noisy_distribution"+ std::to_string(j)+".txt", false);

			emp_pr_z.addi( p_xspace.get_distribution_from_emp_data(list_empdata[j]).mul(ratio)  );

			//std::cout << "after seeing data ("<<j<<"): emp_pr_z = " << emp_pr_z.toString() << "\n";
			//std::cout <<"sums of elements = " << emp_pr_z.sum() << "\n\n";
			delete(mech_ptr); mech_ptr=nullptr;
			//P.S. the full empirical distribution: contains an probability for every element in xspace including those that do not appear in p_emp_data.
			// construct an estimator from the average metrix
			inv = Estimator_INV(pr_z_giv_x, emp_pr_z);  
		}
		//todo: debug
		//exit(0);
		//Prior(emp_pr_z, p_xspace.hsize, p_xspace.vsize)
		//.write_prior_as_matrix("noisy_distribution_avg.txt", false);
	}
}


/**
 * Returns an MLE on the elements of the Full grid secrets space (xspace) using matrix inversion method.
 * using the full empirical data distribution emp_pr_z (computed from xspace), and the mechanism matrix
 * (pr_z_giv_x). The -ve components in the resulting vector are corrected to a valid distribution
 * using normalization 'N' or projection 'P', Or 'O' for no-post processing. 
 */
DoubleRow Estimator_AVG_MECH_INV::getMaximizer (char p_post_processing) const {
	if (all_mechanisms_are_krr)
		return avg_mech_inv_krr.getMaximizer(p_post_processing);
	else
		return inv.getMaximizer(p_post_processing);
}


