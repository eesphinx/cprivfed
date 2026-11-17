/**
 * This class evaluates the maximum likelihood distribution on a given space of secrets (xspace),
 * using multiple empirical data objects. Each object describes a collection of locations and the (nonzero)
 * frequency of every location. Each empirical data object contains the mechanism that produced the empirical data.
 * A GIBU object consists of multiple IBU objects, each one is constructed from one empdata object.
 */
#include <iostream>
#include <cmath>
#include <fstream>
#include "AbstractOutputLocalMechanism.h"
#include "DoubleColumn.h"
#include "EmpData.h"
#include "Estimator_GIBU.h"

/**
 * Construct a GIBU estimator for computing the max-likelihood distribution
 * on a given space of secrets, with respect to empirical grid data objects
 * (with nonzero frequencies) resulting by obfuscation by a given mechanism
 * Every Empdata contains info about the mechanism that produced the data.
 */
Estimator_GIBU::Estimator_GIBU (const FullGridSpace p_xspace, const std::vector<EmpData>& p_list_empdata) {
	// combine empdata objects that have the same mechanism to optimize the performance of GIBU
	std::vector<EmpData> list_empdata = EmpData::compress(p_list_empdata);

  	// we set original distribution with number of columns=0 to check if set or not
  	original_distribution_on_secrets = DoubleRow(0); 
	
	// now work on the new list
	for (size_t j=0; j< list_empdata.size(); j++) {
		//std::cout << "mechanism info: " << p_empdata[j].getObfMechInfo() << "\n";
		//std::cout << "different observations: " << p_empdata[j].data.size() << "\n";
		//std::cout << "total observations: " << p_empdata[j].get_totalObservations() << "\n\n";
		ibus.push_back(Estimator_IBU(p_xspace, list_empdata[j]));
		empdata_total_observations.push_back(list_empdata[j].get_totalObservations());
	}
}

/**
* Construct a GIBU estimator for given <p_empdata> objects, where
* the mechanisms are selected probabilstically according to the conditional 
* probabilities matrix <p_omega>, which maps the elements of xspace (rows) to 
* the mechanisms in the <p_empdata> (columns).
*/
Estimator_GIBU::Estimator_GIBU (const FullGridSpace p_xspace, const std::vector<EmpData>& p_list_empdata, const DoubleMatrix& p_omega) {
	
	// validate p_omega
	if ((p_omega.getRows() != p_xspace.getSize()) || ((size_t) p_omega.getColumns() != p_list_empdata.size())) {
		std::cerr << "error: given p_omega matrix is incompatible with the space or the list of empdata objects";
		exit(1);
	}
	
	// combine empdata objects that have the same mechanism to optimize the performance of GIBU
	std::vector<EmpData> list_empdata = EmpData::compress(p_list_empdata);
	
	// Accordingly combine the columns of p_omega corresponding to clients whose mechanism is the same into one column, 
	// yielding a compressed matrix omega.
	DoubleMatrix omega(p_omega.getRows(),0);
	for (size_t i=0; i<list_empdata.size(); i++) {
		DoubleColumn col(p_omega.getRows());
		for (size_t j=0; j<p_list_empdata.size(); j++) {
			if (p_list_empdata[j].getObfMechInfo() == list_empdata[i].getObfMechInfo() ) {
				col.addi(p_omega.selectColumn(j));
			}
		}
		omega = omega.concatHorizontally(col);
	}

	// debug:
	// omega.write_to_file("omega1.txt");
	// exit(0);

	// we set original distribution with length 0 to check if set or not
  	original_distribution_on_secrets = DoubleRow(0); 
	
	// now work on the new list_empdata and omega
	for (size_t j=0; j< list_empdata.size(); j++) {
		DoubleRow q = list_empdata[j].getDistributionOnObservables(); // q[i] corresponds to the frequency in empdata.data[i]
		Mechanism* mech_ptr = Mechanism::createMechanism(list_empdata[j].getObfMechInfo());
		if (auto local_mech_ptr = dynamic_cast<AbstractOutputLocalMechanism*>(mech_ptr)) {
			DoubleMatrix A(p_xspace.getSize(),0);
			for (size_t z=0; z<list_empdata[j].data.size(); z++) {
				DoubleColumn col = local_mech_ptr->get_mech_column(list_empdata[j].data[z].obs_ptr);
				A = A.concatHorizontally(col);
			}
			A = A.mulColumn(omega.selectColumn(j));
			ibus.push_back(Estimator_IBU(A , q));
			empdata_total_observations.push_back(list_empdata[j].get_totalObservations());
		}
		else {
			std::cerr << "GIBU construction error: Unable to recover the local mechanism";
			exit(1);
		}
	}

}


/** one iteration of GIBU refines p_pr_x */
DoubleRow Estimator_GIBU::update_estimate(const DoubleRow& p_pr_x) const {
	// This transformation is a linear composition of the update procedures of the underlying IBU objects.

	DoubleRow new_pr_x(p_pr_x.getLength());
	new_pr_x.fillWith(0);
	for (size_t j=0; j<ibus.size(); j++) {
		DoubleRow temp = ibus[j].update_estimate(p_pr_x);
		temp.muli(empdata_total_observations[j]);
		new_pr_x.addi(temp);
	}
	return new_pr_x.div(new_pr_x.sum());

}

/** returns the log_likelihood of a given distribution with respect to the mechanism matrix and per one observation
 * P.S. The log_likelihood with respect to all observations is obtained by multiplying the result by number of observations */
double Estimator_GIBU::get_log_likelihood_per_obs(const DoubleRow& p_pr_x) const {
	double ret=0;
	int tot_observations = 0;
	for (size_t j=0; j<ibus.size(); j++) {
		ret = ret + empdata_total_observations[j] * ibus[j].get_log_likelihood_per_obs(p_pr_x);
		tot_observations = tot_observations + empdata_total_observations[j];
	}
	return ret/tot_observations;
}

/**
 * computes the maximum-likelihood estimate on the elements of the secrets space (xspace), using
 * the available empirical distribution on observations emp_pr_z, and the reduced G matrix
 * (pr_z_giv_x). This computation starts with a uniform distribution on xspace. This method is used
 * for debugging as to see the behavior of the EM algorithm on a large number of iterations.
 * @param p_iterations: number of iterations to run
 * @param p_cycle: a number of iterations after which, the distance and likelihood are recorded.
 * @param p_distType: type of the distance to be plotted
 * @param p_xspace_ptr: a pointer to a generic grid space specifying hsize,vsize,cell_side, required if p_distType=EMD
 * @param p_logLikelihoodFileName : file to store the progress of the likelihood with the iterations.
 * @param p_DistanceFileName : to store the progress of the progress of the distance to original distribution.
 * @param p_DifferenceFileName: to store differences (TV) between successive distributions in iterations (slope)
 * @return a distribution on the secrets space, as a row vector
 */
DoubleRow Estimator_GIBU::getMaximizer(int p_iterations, int p_cycle,
		StatDistanceType p_distType, AbstractGridSpace* p_xspace_ptr,
		std::string p_logLikelihoodFileName, std::string p_DistanceFileName, 
		std::string p_DifferenceFileName) const {

	std::ofstream logLikelihoodPlot_fs(p_logLikelihoodFileName);  // to save likelihood plot data
	if (!logLikelihoodPlot_fs.is_open()) {
		std::cerr << "unable to open likelihood plot file: " << p_logLikelihoodFileName << "\n";
		exit(0);
	}

	std::ofstream diffPlot_fs(p_DifferenceFileName);  // to save the TV distance between successive cycles
	if (!diffPlot_fs.is_open()) {
		std::cerr << "unable to open difference plot file: " << p_DifferenceFileName << "\n";
		exit(0);
	}

	// The following is to save the plot data of distance from original plot data.
	// Note that ofstream is only constructed but does not open the file (opened later if needed)
	std::ofstream distancePlot_fs;
	if (original_distribution_on_secrets.getColumns() > 0) {   // i.e. the original distribution has been set
		distancePlot_fs.open(p_DistanceFileName);
		if (!distancePlot_fs.is_open()) {
			std::cerr << "unable to open distance-from-original plot file: " << p_DistanceFileName << "\n";
			exit(0);
		}
	}

	// now the files are successfully open. Proceed to evaluating the MLE

	// Initialize pr_x to uniform distribution on secrets
	int num_secrets = ibus[0].get_number_of_secrets();
	DoubleRow pr_x(num_secrets);
	pr_x.fillWith(1.0/num_secrets);

	// write initial log data to log files
	logLikelihoodPlot_fs << "0 " << get_log_likelihood_per_obs(pr_x) << "\n";

	// if original distribution is available, the method logs the difference to the current distribution.
	if (original_distribution_on_secrets.getColumns() > 0)
		distancePlot_fs << "0 " << StatDistance::getDist(original_distribution_on_secrets , pr_x, p_distType, p_xspace_ptr) << "\n";

	// Iterate the following steps for a given number of iterations or until likelihood difference
	// is small (evaluating likelihoods is not too expensive)
	DoubleRow old_pr_x = pr_x;
	for (int k=1;k<=p_iterations;k++)
	{
		// update pr_x
		pr_x = update_estimate(pr_x);

		// write the iteration results to log files
		logLikelihoodPlot_fs << k << " "<< get_log_likelihood_per_obs(pr_x) << "\n";
		if (k % p_cycle == 0 && original_distribution_on_secrets.getColumns() > 0) {
			distancePlot_fs << k << " " << StatDistance::getDist(original_distribution_on_secrets , pr_x, p_distType, p_xspace_ptr) << "\n";
			std::cout << "Processing iteration " << k << " ... \n";
			double diff =  StatDistance::getTotalVariationDist(pr_x, old_pr_x);
			diffPlot_fs << k <<" "<< diff << "\n";
			old_pr_x = pr_x;
			// write the iteration results to log files
		}
	}
	logLikelihoodPlot_fs.close();
	diffPlot_fs.close();
	distancePlot_fs.close();

	return pr_x;
}


/**
 * Starts with p_init_pr_x and iterates until the total variation between two successive
 * iterations is within p_min_tv or p_max_iterations is reached.
 */
DoubleRow Estimator_GIBU::getMaximizer(const DoubleRow& p_init_pr_x, double p_min_tv, int p_max_iterations) const {

	//todo: debug
	// std::cout << "init distribution = " << p_init_pr_x.toString() << "\n";

	DoubleRow pr_x = p_init_pr_x;  // a row vector
	DoubleRow old_pr_x = pr_x;
	//double likelihood;
	//double old_likelihood = get_log_likelihood_per_obs(old_pr_x);
	double err;
	int i=0;
	do
	{
		i++;  // count the iterations
		//todo: debug
		// std::cout << "Processing iteration  "<< i << "\n";

		// update pr_x
		pr_x = update_estimate(pr_x);
		//likelihood = get_log_likelihood_per_obs(pr_x);


         //err = std::abs(likelihood - old_likelihood);
		err = StatDistance::getTotalVariationDist(pr_x, old_pr_x);

		//todo:debug:
		// std::cout << "error = " << err << "\n";
		// std::cout << "log-likelihood = " << get_log_likelihood_per_obs(pr_x) << "\n";

		old_pr_x = pr_x;
		//old_likelihood = likelihood;

		//todo: debug
		// if (i>200) exit(0);

	} while (err > p_min_tv && i< p_max_iterations);
	std::cout << "GIBU terminated at iteration " << i <<" with error ="<< err
			<< " loglikelihood per observation = "<< get_log_likelihood_per_obs(pr_x) << "\n";
	
	// exit(0);
	return pr_x;
}


/**
 * Starts with the uniform distribution on the secrets space and iterates until the
 * total variation between successive estimates is within p_min_diff or p_max_iterations is reached.
 */
DoubleRow Estimator_GIBU::getMaximizer(double p_min_diff, int p_max_iterations) const {
	// get a uniform distribution on the set of locations described by xspace.
	int num_secrets = ibus[0].get_number_of_secrets();
	DoubleRow init_pr_x(num_secrets);
	init_pr_x.fillWith(1.0/num_secrets);
	return getMaximizer(init_pr_x, p_min_diff, p_max_iterations);
}


/**
 * Starts with the uniform distribution on the secrets and iterates for the given number : p_iterations
 */
DoubleRow Estimator_GIBU::getMaximizer(int p_iterations) const {
	return getMaximizer(0, p_iterations);   // call previous version with p_min_diff=0
}

/**
 * set the original distribution if needed to track the distance to the distribution at every iteration.
 */
void Estimator_GIBU::setOriginalDistributionOnSecrets(const DoubleRow& p_orig_dist) {
	original_distribution_on_secrets = p_orig_dist;
}



