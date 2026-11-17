#include <iostream>
// #include "Mechanism.h"
// #include "PRRespMechanism.h"
#include "Distributions.h"
#include "EmpData.h"
#include "Estimator_AVG_RESULT_INV.h"

/**
 * INV is applied to each emp data object and the resulting distributions are averaged.
 */
Estimator_AVG_RESULT_INV::Estimator_AVG_RESULT_INV (const FullGridSpace p_xspace, const std::vector<EmpData>& p_list_empdata) {
	// combine empdata objects that have the same mechanism 
	std::vector<EmpData> list_empdata = EmpData::compress(p_list_empdata);

	invs.resize(list_empdata.size());
	for (size_t j=0; j<list_empdata.size(); j++) {
		empdata_obs_count.push_back(list_empdata[j].get_totalObservations());
		invs[j] = Estimator_INV(p_xspace, list_empdata[j] );
		mech_infos.push_back(list_empdata[j].getObfMechInfo());  // save the mechanisms
	}
}

/**
 * using normalization 'N' or projection 'P'.
 */
DoubleRow Estimator_AVG_RESULT_INV::getMaximizer (char p_post_processing) const {

	DoubleRow pr_x(invs[0].get_number_of_secrets());
	pr_x.fillWith(0);

	for (size_t i= 0; i< invs.size(); i++)
		pr_x.addi( (invs[i].getMaximizer(p_post_processing)).muli(empdata_obs_count[i])  ); 
	pr_x.divi(pr_x.sum());


	// bool all_krr = false; // assume given mechanisms are kRR until we find a different one
	// std::vector<double> krr_leps;  // record local epsilon for krr mechanisms
	// for (size_t i=0; i<mech_infos.size() && all_krr; i++) {
	// 	Mechanism* mech_ptr = Mechanism::createMechanism(mech_infos[i]);
	// 	if (PRRespMechanism* krr_mech_ptr = dynamic_cast<PRRespMechanism*>(mech_ptr))
	// 		krr_leps.push_back(krr_mech_ptr->get_eps());
	// 	else
	// 		all_krr = false;
	// 	delete mech_ptr; mech_ptr=nullptr;
	// }

	// I abserved that invoking the sub inv estimators with no post-processing incurs bad utility 
	// Therefore I pass the post-processing argument to the sub estimators. 
	
	// if (all_krr) {
	// 	std::cout << "all mechanisms are kRR \n"; 
	// 	int k = pr_x.getLength(); // domain size
	// 	for (size_t i= 0; i< invs.size(); i++) {
	// 		// evaluate an approximation for the kRR inv variance given in 
	// 		// Kairouz paper (proposition 3). The variance v_d in the proposition 
	// 		// is for the estimated distribution (h/n). That is Var(h/n)= vd.
	// 		// Therefore the variance for the histogram h is Var(h)= n^2 * vd 
	// 		double t = exp(krr_leps[i]) - 1; 
	// 		double var = empdata_obs_count[i]*(1+(k-1)*(k + 2*t)/t*t); 
	// 		// In the following we multiply by empdata_total_observation[i] to get h[i]
	// 		// and then multiply by var
	// 		pr_x.addi( (invs[i].getMaximizer('O')).muli(empdata_obs_count[i]/var)  ); 
	// 	}
	// 	pr_x.divi(pr_x.sum());
	// }
	// else {
	// 	std::cout << "all mechanisms are not kRR \n"; 
	// 	for (size_t i= 0; i< invs.size(); i++)
	// 		pr_x.addi( (invs[i].getMaximizer(p_post_processing)).muli(empdata_obs_count[i])  ); 
	// 	pr_x.divi(pr_x.sum());
	// }
	
	return pr_x;
}
