/**
 * This class runs IBU on every given empirical data and returns the weighted average of resulting distributions
 */

#include <iostream>
#include "EmpData.h"
#include "Estimator_AVG_RESULT_IBU.h"
#include "GridObservation.h"

/**
 * IBU is applied to each emp data object and the resulting distributions are averaged.
 */
Estimator_AVG_RESULT_IBU::Estimator_AVG_RESULT_IBU (const FullGridSpace p_xspace, const std::vector<EmpData>& p_list_empdata) {
	
	// combine empdata objects that have the same mechanism 
	std::vector<EmpData> list_empdata = EmpData::compress(p_list_empdata);
	for (size_t j=0; j<list_empdata.size(); j++) {
		ibus.push_back(Estimator_IBU(p_xspace, list_empdata[j]));
		empdata_total_observations.push_back(list_empdata[j].get_totalObservations());
	}


// 	// start
// 	// modified version to handle nonhomogeneous sampling. 
// 	// We assume we have only two mechanisms/ibus: ibu[0] processes the first half of the xspace (x < length/2), 
// 	// and ibu[1] processes the upper half (x>= length/2).
// 	// Here we remove for each mechanism, the irrelevant observations
// 	std::vector<EmpData> list_empdata = EmpData::compress(p_list_empdata);
// 	// std::cout << "size of data[0] before filtering: " << list_empdata[0].get_totalObservations() << "\n"; 
// 	// std::cout << "size of data[1] before filtering: " << list_empdata[1].get_totalObservations() << "\n"; 
// 	// list_empdata[0].write_data_to_file("data0_before.txt");
// 	// list_empdata[1].write_data_to_file("data1_before.txt");
// 	if (list_empdata.size()!=2) {std::cerr << "error: number of empdata objects is not 2"; exit(1);}
// 	for (size_t l=0; l<list_empdata[0].data.size(); l++) 
// 		if (GridObservation* grid_obs_ptr = dynamic_cast<GridObservation*>(list_empdata[0].data[l].obs_ptr)) 
// 			if (grid_obs_ptr->location>=p_xspace.getSize()/2 ){
// 				list_empdata[0].data[l].frequency=0;
// 			}
// 	for (size_t l=0; l<list_empdata[1].data.size(); l++)
// 		if (GridObservation* grid_obs_ptr = dynamic_cast<GridObservation*>(list_empdata[1].data[l].obs_ptr)) 
// 			if (grid_obs_ptr->location<p_xspace.getSize()/2 ){
// 				list_empdata[1].data[l].frequency=0;
// 			}
// 	for (size_t j=0; j<list_empdata.size(); j++) {
// 		ibus.push_back(Estimator_IBU(p_xspace, list_empdata[j]));
// 		empdata_total_observations.push_back(list_empdata[j].get_totalObservations());
// 	}
// 	// std::cout << "size of data[0] after filtering: " << list_empdata[0].get_totalObservations() << "\n"; 
// 	// std::cout << "size of data[1] after filtering: " << list_empdata[1].get_totalObservations() << "\n"; 
// 	// list_empdata[0].write_data_to_file("data0_after.txt");
// 	// list_empdata[1].write_data_to_file("data1_after.txt");
// 	// exit(0);
// //end


}

/**
 * Starts with the uniform distribution on the secrets space and iterates until the
 * total variation between successive estimates is within p_min_diff or p_max_iterations is reached.
 */
DoubleRow Estimator_AVG_RESULT_IBU::getMaximizer(double p_min_diff, int p_max_iterations) const {

	DoubleRow pr_x(ibus[0].get_number_of_secrets());
	pr_x.fillWith(0);

	for (size_t i= 0; i< ibus.size(); i++)
		pr_x.addi( (ibus[i].getMaximizer(p_min_diff, p_max_iterations)).muli( empdata_total_observations[i] ) );

	return pr_x.div(pr_x.sum());

// // start
// 	// modified version to handle nonhomogeneous sampling. 
// 	// We assume we have only two mechanisms/ibus: ibu[0] processes the first half of the xspace (x < length/2), 
// 	// and ibu[1] processes the upper half (x>= length/2).
// 	if (ibus.size()!=2) {std::cerr << "error: number of mechanisms/ibus is not 2"; exit(1);}
// 	std::vector<DoubleRow> list_init_pr_x;
// 	DoubleRow t(ibus[0].get_number_of_secrets()); t.fillWith(1.0);
// 	for (int j=t.getLength()/2; j<t.getLength(); j++)  t.put(j,0);
// 	list_init_pr_x.push_back(t.div(t.sum()));
// 	t.subi(1.0).muli(-1.0);
// 	list_init_pr_x.push_back(t.div(t.sum()));
// 	DoubleRow pr_x(ibus[0].get_number_of_secrets());
// 	pr_x.fillWith(0);
// 	for (size_t i= 0; i< ibus.size(); i++)
// 	 	pr_x.addi( (ibus[i].getMaximizer(list_init_pr_x[i], p_min_diff, p_max_iterations)).muli(empdata_total_observations[i])  );
// 	return pr_x.div(pr_x.sum());
// //end


}

