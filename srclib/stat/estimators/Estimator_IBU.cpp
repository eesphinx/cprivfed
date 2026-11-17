/**
 * This class evaluates the maximum likelihood distribution on a given space of secrets (xspace),
 * using empirical data which describes a collection of generic observations and the (nonzero)
 * frequency of every one, and the mechanism that produced this empirical data.
 */
#include <iostream>
#include <fstream>
#include "DoubleColumn.h"
#include "AbstractOutputLocalMechanism.h"
#include "EmpData.h"
#include "Estimator_IBU.h"

/** Constructs an IBU object from a matrix A (of probabilities) and empirical distribution q */
Estimator_IBU::Estimator_IBU(const DoubleMatrix& A, const DoubleRow& q) {
	if (A.getColumns() != q.getLength()) {
		std::cerr << "error: given A and q must have the same columns";
		exit(0);
	}
	/* Execluding elements z with q(z)=0 from both q and A, to improve the performance */
	std::vector<int> nz_indexes;
	for (int i=0; i<q.getLength(); i++)
		if (q.get(i)>0) nz_indexes.push_back(i);
	emp_pr_z = q.selectColumns(nz_indexes);
	pr_z_giv_x = A.selectColumns(nz_indexes);
}


/**
 * Construct an IBU estimator for computing the max-likelihood distribution
 * on a given space of secrets, with respect to empirical grid data
 * (with nonzero frequencies) resulting by obfuscation by the mechanism (specified in the empdata object)
 * The mechanism information is included in the empdata object
 */
Estimator_IBU::Estimator_IBU (const FullGridSpace p_xspace, const EmpData& p_empdata)
	:reference_distribution_on_secrets(0) { // we set reference distribution with columns=0 to check if set or not

	if (p_empdata.getSize() == 0) {
		std::cerr << "unable to construct IBU estimator: the data object is empty";
		exit(0);
	}
	// to debug
	// std::cout << "Inside IBU \n";
	// std::cout << "p_empdata: totalObs = " << p_empdata.get_totalObservations() << ", mech = " << p_empdata.getObfMechInfo() << "\n";
	// std::cout << "empirical distribution: " << p_empdata.getDistributionOnObservables() << "\n";

	Mechanism* mech_ptr = Mechanism::createMechanism(p_empdata.getObfMechInfo());
	AbstractOutputLocalMechanism* local_mech_ptr = dynamic_cast<AbstractOutputLocalMechanism*>(mech_ptr);
	if (!local_mech_ptr) {  // the given mechanism is not local
		std::cerr << "unable to construct IBU estimator: the mechanism of data object is not local";
		delete(mech_ptr); mech_ptr=nullptr;
		exit(1);
	}
	if (local_mech_ptr->hsize != p_xspace.hsize || local_mech_ptr->vsize != p_xspace.vsize) {
		std::cerr << "unable to construct IBU estimator: the grid of mechanism ("<< local_mech_ptr->getTypeStr()
				<<") that produced the empirical data is incompatible with the given space";
		delete(mech_ptr); mech_ptr=nullptr;
		exit(1);
	}
	// Now everything should be fine

	// get the mechanism matrix columns p_empdata, i.e. a slim version of the matrix G
	// called 'slim' because (unlike normal G) no column is repeated.
	// Every column corresponds to one datum (with nonzero frequency).
	pr_z_giv_x = local_mech_ptr->get_mech_column(p_empdata.data[0].obs_ptr);
	for (size_t i=1; i< p_empdata.getSize();i++)
		pr_z_giv_x = pr_z_giv_x.concatHorizontally(local_mech_ptr->get_mech_column(p_empdata.data[i].obs_ptr));

	// This allows secret spaces to be subset of the full grid.
	// This might impact the performance.
	std::vector<int> secretsRows = p_xspace.getSecretsIndexes();
	pr_z_giv_x.selectiRows(secretsRows);
	delete(mech_ptr); mech_ptr = nullptr;

	// The following is faster but specific to grid-outputs mechanisms.
	//std::vector<int> spaceRows = p_xspace_ptr->getSecretsIndexes();
	//std::vector<int> obsCols = p_empdata.getObservables();
	//pr_z_giv_x = mechPtr->get_mech_matrix().selectiRows(spaceRows).selectiColumns(obsCols);


	// the empirical distribution: every element of emp_pr_z has a nonzero value.
	// if an observable does not show up in p_empdata, then emp_pr_z does not include an entry.
	// This is important because the components of emp_pr_z must correspond to the columns of pr_z_giv_x given above.
	emp_pr_z = p_empdata.getDistributionOnObservables();
}

/*
 * one iteration of IBU: transforms p_pr_x in the light of the nonzero empirical distribution p_pr_z_giv_x
 * and the probability matrix p_pr_z_giv_x (restricted to elements in p_emp_pr_z)
 * P.S. If for some x, p_pr_x[x] = 0, then the update will have alsp pr[x]=0. 
 * */
DoubleRow Estimator_IBU::update_estimate(const DoubleRow& p_pr_x) const {
	// Compute a row vector containing the probability of
	// every z in the empirical data (only observed ones, not every z in Z)
	// with respect to the current distribution, i.e. pr_z = pr_x.G
	// We exclude non-observed elements of Z to avoid division by zero. 

	//todo: debug
	//std::cout << "p_pr_x = "<< p_pr_x.toString() << "\n";
	DoubleRow pr_z = p_pr_x.mmul(pr_z_giv_x);
	//std::cout << "pr_z = " << pr_z.toString() << "\n";

	// Construct a joint probability matrix where each entry is p(x,z)
	// multiply (point-wise) every column of pr_z_giv_x by pr_x in a column form, i.e. transposed.
	DoubleMatrix pr_x_and_z = pr_z_giv_x.mulColumn(p_pr_x.transpose());
	//std::cout << "pr_x_and_z = \n" << pr_x_and_z.toString()<< "\n";

	// Construct a conditional matrix of the form P(x|z), with z indexing rows and x indexing columns
	DoubleMatrix pr_x_giv_z = (pr_x_and_z.divRow(pr_z)).transpose();
	//std::cout << "pr_x_giv_z = \n" << pr_x_giv_z.toString() << "\n";

	// update the current distribution
	return emp_pr_z.mmul(pr_x_giv_z);
}

// /** returns the log_likelihood of a given distribution with respect to the mechanism matrix and per one observation
//  * P.S. The log_likelihood with respect to all observations, the result has to be multiplied by number of observations
//  * It is made public because it is accessed by GIBU to compute its likelihood, and also from experiments to inspect the
//  * likelihood of estimated distributions */
// double Estimator_IBU::get_log_likelihood_per_obs(const DoubleRow& p_pr_x) const {
// 	 return Likelihood::get_log_likelihood(p_pr_x, pr_z_giv_x, emp_pr_z, 1);
// }

/** returns the number of secrets on which a distribution is estimated. Used to initialize distributions */
int Estimator_IBU::get_number_of_secrets() const {
	return pr_z_giv_x.getRows();
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
 * @param p_DistanceFileName : to store the progress of the progress of the distance to the reference distribution.
 * @param p_DifferenceFileName: to store differences (TV) between successive distributions in iterations (slope)
 * @return a distribution on the secrets space, as a row vector
 */
DoubleRow Estimator_IBU::getMaximizer(int p_iterations, int p_cycle,
		StatDistanceType p_distType, AbstractGridSpace* p_xspace_ptr,
		std::string p_logLikelihoodFileName, std::string p_DistanceFileName, 
		std::string p_DifferenceFileName) const {

	std::ofstream logLikelihoodPlot_fs(p_logLikelihoodFileName);  // to save likelihood plot data
	if (!p_logLikelihoodFileName.empty()) {
		if (!logLikelihoodPlot_fs.is_open()) {
			std::cerr << "unable to open likelihood plot file: " << p_logLikelihoodFileName << "\n";
			exit(0);
		}
	}

	std::ofstream diffPlot_fs(p_DifferenceFileName);  // to save the TV distance between successive cycles
	if (!p_logLikelihoodFileName.empty()) {
		if (!diffPlot_fs.is_open()) {
			std::cerr << "unable to open difference plot file: " << p_DifferenceFileName << "\n";
			exit(0);
		}
	}

	// The following is to save the distance from the reference distribution 
	// into the distance file to allow plotting this distance.
	// Note that ofstream is only constructed but does not open the file (opened later if needed)
	std::ofstream distancePlot_fs;
	if (!p_DistanceFileName.empty()) {
		if (reference_distribution_on_secrets.getColumns() > 0) {   // i.e. the reference distribution has been set
			distancePlot_fs.open(p_DistanceFileName);
			if (!distancePlot_fs.is_open()) {
				std::cerr << "unable to open distance-from-reference plot file: " << p_DistanceFileName << "\n";
				exit(0);
			}
		}
	}
	// now the files are successfully open. Proceed to evaluating the MLE


	// Initialize pr_x to uniform distribution on secrets
	int num_secrets = get_number_of_secrets();
	DoubleRow pr_x(num_secrets);
	pr_x.fillWith(1.0/num_secrets);

	// write initial log data to log files
	if (!p_logLikelihoodFileName.empty()) {
		logLikelihoodPlot_fs << "0 " << get_log_likelihood_per_obs(pr_x) << "\n";
	}

	// if reference distribution is available, the method logs the difference to the current distribution.
	if (!p_DistanceFileName.empty() && 
	    reference_distribution_on_secrets.getColumns() > 0) {
		// to debug
		//std::cout << "step 0:\n";
		//std::cout << "reference distribution:\n";
		//std::cout << reference_distribution_on_secrets.toString()<< "\n";
		//std::cout << "ibu distribution at setp 0:\n";
		//std::cout << pr_x.toString()<< "\n";
		//std::cout << "distance = "
		//		  << StatDistance::getDist(reference_distribution_on_secrets , pr_x, p_distType, p_xspace_ptr);
		//std::cout << "\n\n";

		distancePlot_fs << "0 " << StatDistance::getDist(reference_distribution_on_secrets , pr_x, p_distType, p_xspace_ptr) << "\n";
	}

	// Iterate the following steps for a given number of iterations or until likelihood difference
	// is small (evaluating likelihoods is not too expensive)
	DoubleRow old_pr_x = pr_x;
	for (int k=1;k<=p_iterations;k++)
	{
		// update pr_x
		pr_x = update_estimate(pr_x);
		// write the iteration results to log files
		
		// save likelihood data if required
		if (!p_logLikelihoodFileName.empty())
			logLikelihoodPlot_fs << k << " "<< 
			get_log_likelihood_per_obs(pr_x)<<"\n";
		
		// save distances to the reference and previous distributions
		// at the beginnig of every cycle 
		if (k % p_cycle == 0) {
			//std::cout << "Processing iteration " << k << " ... \n";

		    // save the distance to the reference distribution if required
			if (!p_DistanceFileName.empty() && 
	    			reference_distribution_on_secrets.getColumns() > 0) {
				distancePlot_fs << k << " " << StatDistance::getDist(reference_distribution_on_secrets , pr_x, p_distType, p_xspace_ptr) << "\n";
				// to debug
				//std::cout << "reference distribution:\n";
				//std::cout << reference_distribution_on_secrets.toString()<< "\n";
				//std::cout << "ibu distribution:\n";
				//std::cout << pr_x.toString()<< "\n";
				//std::cout << "distance = "
				//  << StatDistance::getDist(reference_distribution_on_secrets , pr_x, p_distType, p_xspace_ptr);
				//std::cout << "\n\n";
			}
			
			// save the difference to previous distribution if required
			if (!p_DifferenceFileName.empty()) {
				double diff =  StatDistance::getTotalVariationDist(pr_x, old_pr_x);
				diffPlot_fs << k <<" "<< diff << "\n";
				old_pr_x = pr_x;
				// write the iteration results to log files
			}
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
DoubleRow Estimator_IBU::getMaximizer(const DoubleRow& p_init_pr_x, double p_min_tv, int p_max_iterations) const {

	//todo: debug
	//std::cout << "emp_pr_z = " << emp_pr_z.toString() << "\n";
	//std::cout << "pr_z_giv_x = \n" << pr_z_giv_x.toString() << "\n";
	//std::cout << "init distribution = " << p_init_pr_x.toString() << "\n";

	DoubleRow pr_x = p_init_pr_x;  // a row vector
	DoubleRow old_pr_x = pr_x;
	double tv_diff;
	double likelihood_per_obs_diff;
	int i=0;
	do
	{
		i++;  // count the iterations
		//todo: debug
		// std::cout << "Processing iteration  "<< i << "\n";

		// update pr_x
		pr_x = update_estimate(pr_x);
		tv_diff = StatDistance::getTotalVariationDist(pr_x, old_pr_x);
		likelihood_per_obs_diff = get_log_likelihood_per_obs(pr_x)
						- get_log_likelihood_per_obs(old_pr_x);
		//todo:debug:
		// std::cout << "TV error = " << tv_diff << "\n";
		old_pr_x = pr_x;
		//todo: debug
		//if (i>4) exit(0);

		//FullGridSpace s(p_init_pr_x.getLength(),1,1);
		//std::cout << "EMD distance from original = " << StatDistance::getGeoEMDDist(pr_x,
		//		             original_distribution_on_secrets, &s ) << "\n";
		//   //StatDistance::getTotalVariationDist(pr_x, original_distribution_on_secrets) << "\n";
		// if (i%1000==0) {
		// std::cout << "iteration: " << i << "\n";
		// std::cout << "TV distance from previous estimate = " << tv_diff << "\n";
		// // std::cout << "Likelihood/obs = " << get_log_likelihood_per_obs(pr_x) << "\n\n";
		// std::cout << "Likelihood/obs Diff = " << likelihood_per_obs_diff << "\n\n";
		// }
		//

	} while (tv_diff > p_min_tv && i< p_max_iterations);
	std::cout << "IBU terminated at iteration " << i <<" with total variation diff ="<< tv_diff
			<< " Log Likelihood difference = "<< likelihood_per_obs_diff << "\n"
			<< " loglikelihood per observation = "<< get_log_likelihood_per_obs(pr_x) << "\n";
	return pr_x;
}

/**
 * Starts with p_init_pr_x and iterates until the increase of log_likelihood between two successive
 * iterations is within p_max_logLikelihood_dff or p_max_iterations is reached.
 */
DoubleRow Estimator_IBU::getMaximizer_by_logLikelihood_per_obs(const DoubleRow& p_init_pr_x, 
	double p_max_logLikelihood_per_obs_dff, int p_max_iterations) const {

	//todo: debug
	//std::cout << "emp_pr_z = " << emp_pr_z.toString() << "\n";
	//std::cout << "pr_z_giv_x = \n" << pr_z_giv_x.toString() << "\n";
	//std::cout << "init distribution = " << p_init_pr_x.toString() << "\n";

	DoubleRow pr_x = p_init_pr_x;  // a row vector
	DoubleRow old_pr_x = pr_x;
	double tv_diff;
	double logLikelihood_per_obs_diff;
	int i=0;
	do
	{
		i++;  // count the iterations
		//todo: debug
		//std::cout << "Processing iteration  "<< i << "\n";

		// update pr_x
		pr_x = update_estimate(pr_x);
		tv_diff = StatDistance::getTotalVariationDist(pr_x, old_pr_x);

		// debug:
		// std::cout << "run: " << i << ", TV = " << tv_diff << "\n";

		// std::cout<< "old likelihood = " << get_log_likelihood_per_obs(old_pr_x) << "\n";
		// std::cout<< "new likelihood = " << get_log_likelihood_per_obs(pr_x) << "\n";
		logLikelihood_per_obs_diff = get_log_likelihood_per_obs(pr_x)
						- get_log_likelihood_per_obs(old_pr_x);
		//todo:debug:
		//std::cout << "em error = " << err << "\n";
		old_pr_x = pr_x;
		//todo: debug
		//if (i>4) exit(0);

		// to debug: 17/11/2024
		//FullGridSpace s(p_init_pr_x.getLength(),1,1);
		//std::cout << "EMD distance from original = " << StatDistance::getGeoEMDDist(pr_x,
		//		             original_distribution_on_secrets, &s ) << "\n";
		//   //StatDistance::getTotalVariationDist(pr_x, original_distribution_on_secrets) << "\n";
		// if (i%100==0) {
		// std::cout << "iteration: " << i << "\n";
		// std::cout << "TV distance from previous estimate = " << tv_diff << "\n";
		// std::cout << "Likelihood/obs = " << get_log_likelihood_per_obs(pr_x) << "\n\n";
		// std::cout << "Likelihood/obs Diff = " << logLikelihood_per_obs_diff << "\n\n";
		// }
		//

	} while (logLikelihood_per_obs_diff > p_max_logLikelihood_per_obs_dff && i< p_max_iterations);
	std::cout << "IBU terminated at iteration " << i <<" with log_likelikood_per_obs_diff ="<< logLikelihood_per_obs_diff
			<< " loglikelihood per observation = "<< get_log_likelihood_per_obs(pr_x) << "\n"
			<< " Total Variation Diff = "<< tv_diff << "\n";
	
	return pr_x;
}



/**
 * Starts with the uniform distribution on the secrets space and iterates until the
 * total variation between successive estimates is within p_max_tv or p_max_iterations is reached.
 */
DoubleRow Estimator_IBU::getMaximizer(double p_max_tv, int p_max_iterations) const {
	// Initialize pr_x to uniform distribution on secrets
	int num_secrets = get_number_of_secrets();
	DoubleRow init_pr_x(num_secrets);
	init_pr_x.fillWith(1.0/num_secrets);
	return getMaximizer(init_pr_x, p_max_tv, p_max_iterations);
}


/**
 * Starts with the uniform distribution on the secrets and iterates for the given number : p_iterations
 */
DoubleRow Estimator_IBU::getMaximizer(int p_iterations) const {
	return getMaximizer(0, p_iterations);   // call previous version with p_min_diff=0
}

/** returns the log_likelihood of a given distribution with respect to the probability 
 * matrix pr_z_giv_x and per one observation. P.S. To obtain the log_likelihood with 
 * respect to all observations, the result has to be multiplied by number of observations */
double Estimator_IBU::get_log_likelihood_per_obs(const DoubleRow& p_pr_x) const {
	double ret=0;
	for (int j=0;j<emp_pr_z.getLength();j++) {
		DoubleColumn col = pr_z_giv_x.selectColumn(j);
		double d = p_pr_x.dot( col.transpose() );  // the two vectors/matrices have to have the same dimensions.
		ret = ret + emp_pr_z.get(j)* log(d );
	}
	return ret;
}

/**
 * set a reference distribution if needed to track the distance to this distribution at every iteration.
 */
void Estimator_IBU::setReferencelDistributionOnSecrets(const DoubleRow& p_ref_distribution) {
	reference_distribution_on_secrets = p_ref_distribution;
}


DoubleMatrix Estimator_IBU::getAsympFisher(const DoubleRow& p_pr_x, int observationCount) const {
	// compute the expected Fisher matrix I(p_pr_x)
	std::vector<int> rows;
//	int k = p_pr_x.getLength()/2;
//	for (int x=0  ; x < k ; x++)  rows.push_back(x);
//	for (int x=k+1; x < p_pr_x.getLength(); x++)  rows.push_back(x);
	int k = p_pr_x.getLength()-1;
	for (int x=0  ; x < k ; x++)  rows.push_back(x);

	//for (int x=0; x<p_pr_x.getLength()-1; x++)  rows.push_back(x);

	DoubleMatrix fisherMatrix(p_pr_x.getLength()-1, p_pr_x.getLength()-1);
	fisherMatrix.fillWith(0);
	// for every observation <obs> construct corresponding FisherMatrix
	for (int z=0; z<pr_z_giv_x.getColumns(); z++) {

//		std::cout << " z = " << z << "\n";
        //  way 1: using the given distribution p_pr_x to compute p(z). This allows using the original distribution
		DoubleColumn score = pr_z_giv_x.selectColumn(z); // col
//		std::cout << "col(z) = " << score.toString() << "\n";
		double pz = p_pr_x.transpose().dot(score);      // sum_x theta(x) col(x)
//		std::cout << "pz = " << pz << "\n";
		double a_kz = score.get(k);
        // formulate score column
		score.selectiRows(rows).subi(a_kz).divi(pz);

		DoubleMatrix m = score.mul(pz).mmul(score.transpose());


//		std::cout << "fisher(z) = " << m.toString() << "\n";
		fisherMatrix.addi(m);


		////////////

		//  way 2: using the empirical p(z) directly as approximation of the exact p(z).
		//  This allows estimating variance with no knowledge about the real distribution
//		DoubleColumn col = pr_z_giv_x.selectColumn(z); // col
//		double d = p_pr_x.transpose().dot(col);        // sum_x theta(x) col(x)
//		col.divi(d);                                   // score column
//		DoubleMatrix m = col.mul(emp_pr_z.get(z)).mmul(col.transpose());
//		fisherMatrix.addi(m);
		////////////

	}
	return fisherMatrix;

}

DoubleMatrix Estimator_IBU::getAsympFisherInv(const DoubleRow& p_pr_x, int observationCount) const {
	return getAsympFisher(p_pr_x,observationCount).getInverse();
}



/** returns an approximate variance, also called the expected mean squared error (MSE), also
 * called the expected norm-2 error of the MLE for the given mechanism and data based on the
 * asymptotic normality of the MLE. The argument p_pr_x may be the true distribution or an MLE
 * obtained from IBU  */
double Estimator_IBU::getAsympRootVariance(const DoubleRow& p_pr_x, int observationCount) const {
	DoubleMatrix cov = getAsympFisherInv(p_pr_x, observationCount);
	//  std::cout << "variance of sum_x theta_x = " << cov.sum() << "\n";
	return sqrt((cov.getDiag().sum()+ cov.sum() )/observationCount);
}

/** returns an approximate variance, also called the expected mean squared error (MSE), also
 * called the expected norm-2 error of the MLE for the given mechanism and data based on the
 * asymptotic normality of the MLE. The argument p_pr_x may be the true distribution or an MLE
 * obtained from IBU  */
double Estimator_IBU::getApproxVariance1(const DoubleRow& p_pr_x, int observationCount) const {
	// compute the expected Fisher matrix I(p_pr_x)
	DoubleMatrix fisherMatrix(p_pr_x.getLength(),p_pr_x.getLength());
	fisherMatrix.fillWith(0);
	// for every observation <obs> construct corresponding FisherMatrix
	for (int z=0; z<pr_z_giv_x.getColumns(); z++) {

//		std::cout << " z = " << z << "\n";
        //  way 1: using the given distribution p_pr_x to compute p(z). This allows using the original distribution
		DoubleColumn col = pr_z_giv_x.selectColumn(z); // col

//		std::cout << "col(z) = " << col.toString() << "\n";
		double d = p_pr_x.transpose().dot(col);        // sum_x theta(x) col(x)
//		std::cout << "d = " << d << "\n";
		col.divi(d);                                   // score column
		DoubleMatrix m = col.mul(d).mmul(col.transpose());

//		std::cout << "fisher(z) = " << m.toString() << "\n";
		fisherMatrix.addi(m);


		////////////

		//  way 2: using the empirical p(z) directly as approximation of the exact p(z).
		//  This allows estimating variance with no knowledge about the real distribution
//		DoubleColumn col = pr_z_giv_x.selectColumn(z); // col
//		double d = p_pr_x.transpose().dot(col);        // sum_x theta(x) col(x)
//		col.divi(d);                                   // score column
//		DoubleMatrix m = col.mul(emp_pr_z.get(z)).mmul(col.transpose());
//		fisherMatrix.addi(m);
		////////////

	}

//	std::cout << "fisher = " << fisherMatrix.toString() << "\n";
	// get the inverse matrix
	DoubleMatrix cov = fisherMatrix.getInverse();
//	std::cout << "fisher inverse = " << cov.toString() << "\n";
//	std::cout << "fisher.fisher^-1 = " << fisherMatrix.mmul(cov).toString() << "\n";
//	std::cout << "variance of sum_x theta_x = " << cov.sum() << "\n";
	return sqrt((cov.getDiag().sum())/observationCount);
}


