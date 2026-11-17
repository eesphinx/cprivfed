/** constructs an average mechanism from the individual mechanisms in 
 * p_empdata and run IBU on the average mechanism. This estimator works 
 * only for GridOutputMechanisms so that the average mechanism can be constructed */
#include <iostream>
#include "GridOutputLocalMechanism.h"
#include "EmpData.h"
#include "Estimator_IBU.h"
#include "Estimator_AVG_MECH_IBU.h"

/** constructs an average mechanism from the individual mechanisms in p_empdata and run IBU on the average mechanism */
Estimator_AVG_MECH_IBU::Estimator_AVG_MECH_IBU (const FullGridSpace p_xspace, const std::vector<EmpData>& p_list_empdata) {

	// combine empdata objects that have the same mechanism 
	std::vector<EmpData> list_empdata = EmpData::compress(p_list_empdata);

	DoubleMatrix avg_mech_matrix(p_xspace.hsize*p_xspace.vsize, p_xspace.hsize*p_xspace.vsize);
    avg_mech_matrix.fillWith(0);
    DoubleRow  avg_q(p_xspace.hsize*p_xspace.vsize);
    avg_q.fillWith(0);

	// compute total observations
	int total_observations =0;
	for (size_t j=0; j<list_empdata.size(); j++)
		total_observations = total_observations + list_empdata[j].get_totalObservations();

	for (size_t j=0; j<list_empdata.size(); j++) {
		if (list_empdata[j].getSize() == 0) {
			std::cerr << "unable to construct AVG_MECH_IBU: data object " << j << " is empty";
			exit(1);
		}
		Mechanism* mech_ptr = Mechanism::createMechanism(list_empdata[j].getObfMechInfo());
		GridOutputLocalMechanism* local_mech_ptr = dynamic_cast<GridOutputLocalMechanism*>(mech_ptr);
		if (!local_mech_ptr) {  // the given mechanism is not local
			std::cerr << "unable to construct AVG_MECH_IBU: the mechanism of data object " << j << " is not local grid output mechanism";
			delete(mech_ptr); mech_ptr=nullptr;
			exit(1);
		}
		if (local_mech_ptr->hsize != p_xspace.hsize || local_mech_ptr->vsize != p_xspace.vsize) {
			std::cerr << "unable to construct an AVG_MECH_IBU estimator: "
					<< "The grid of mechanism ("<< local_mech_ptr->getTypeStr()<<")"
					<< "that produced the empirical data object " << j << "is not compatible with the given space";
			delete(mech_ptr); mech_ptr=nullptr;
			exit(1);
		}
		// Now everything should be fine
		double m = (double) list_empdata[j].get_totalObservations()/total_observations;
		//std::cout << "ratio (m) for data object " << j << " = " << m << "\n";
		avg_mech_matrix.addi(local_mech_ptr->get_mech_matrix().muli(m));
		avg_q.addi(p_xspace.get_distribution_from_emp_data(list_empdata[j]).muli(m));
		delete(mech_ptr); mech_ptr=nullptr;
	}

	ibu = Estimator_IBU(avg_mech_matrix, avg_q);
}

/**
 * Starts with p_init_pr_x and iterates until the total variation between two successive
 * iterations is within p_min_diff or p_max_iterations is reached.
 */
DoubleRow Estimator_AVG_MECH_IBU::getMaximizer(const DoubleRow& p_init_pr_x, double p_min_diff, int p_max_iterations) const {
	return ibu.getMaximizer(p_init_pr_x,p_min_diff,p_max_iterations);
}


/**
 * Starts with the uniform distribution on the secrets space and iterates until the
 * total variation between successive estimates is within p_min_diff or p_max_iterations is reached.
 */
DoubleRow Estimator_AVG_MECH_IBU::getMaximizer(double p_min_diff, int p_max_iterations) const {
	return ibu.getMaximizer(p_min_diff, p_max_iterations);
}


/**
 * Starts with the uniform distribution on the secrets and iterates for the given number : p_iterations
 */
DoubleRow Estimator_AVG_MECH_IBU::getMaximizer(int p_iterations) const {
	return ibu.getMaximizer(p_iterations);   // call previous version with p_min_diff=0
}

