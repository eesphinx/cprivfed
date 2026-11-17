#include <iostream>
#include "EmpData.h"
#include "Estimator_AVG_RESULT_RAPPOR.h"

/**
 * Estimator_RAPPOR  is applied to each emp data object and the resulting distributions are averaged.
 */
Estimator_AVG_RESULT_RAPPOR::Estimator_AVG_RESULT_RAPPOR (const FullGridSpace p_xspace, const std::vector<EmpData>& p_list_empdata) {
	
	// combine empdata objects that have the same mechanism 
	std::vector<EmpData> list_empdata = EmpData::compress(p_list_empdata);

	raps.resize(list_empdata.size());
	weights.resize(list_empdata.size());
	std::vector<int> empdata_total_observation(list_empdata.size());   // temp vector to hold the obs number of each empdata

	int total_obs=0;
	for (size_t j=0; j<list_empdata.size(); j++) {
		empdata_total_observation[j] = list_empdata[j].get_totalObservations();
		total_obs = total_obs + empdata_total_observation[j];
	}

	for (size_t j=0; j<list_empdata.size(); j++) {
		raps[j] = Estimator_RAPPOR(p_xspace, list_empdata[j] );
		weights[j] = (double) empdata_total_observation[j]/total_obs;
	}
}

/**
 * using normalization 'N' or projection 'P'.
 */
DoubleRow Estimator_AVG_RESULT_RAPPOR::getMaximizer (char p_post_processing) const {

	DoubleRow pr_x(raps[0].get_number_of_secrets());
	pr_x.fillWith(0);

	for (size_t i= 0; i< raps.size(); i++)
		pr_x.addi( (raps[i].getEstimator(p_post_processing)).muli(weights[i])  );

	return pr_x;
}
