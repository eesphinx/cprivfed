/**
* This receives the clients configuration, represented by a list of records.
 * Each one has the form <mechanism:dssize:freq> describing a central mechanism, 
 * a dataset size, and the number of clients having this pair. This may be 
 * constructed by associating one record (mechanism,dssize,1) to each client, 
 * and then grouping these records on (mechanism,dssize).  
 * 
 * The data points of each client are sampled from a master_distribution. 
 * The overall histogram of data is the sum of histograms of individual datasets, 
 * and the original distribution is obtained by normalizing this histogram.
 * Each client obfuscates its points using its local mechanism to produce a noisy histogram, 
 * then the original distribution may be estimated from the noisy histograms using:
 *     a. AVG_MECH_INV, AVG_MECH_IBU
 *     b. AVG_MECH_INV_KRR (with projection or normalization when all the mechanisms are kRR),
 *     c. GIBU
 *     d. AVG_RESULT_INV, AVG_RESULT_IBU
 *     e. COMBINED_RAPPOR
 * The sampling-obfuscation-estimations procedure is repeated <runs> times. 
 * The statistical distance, and optinally the execution time, for each run 
 * is recorded in a file corresponding to each estimator.
 
 * This program receives a set of central mechanisms (CLaplace, CGauss, or CGeom), 
 * a number of datasets (numds), and the average size of each dataset (dssize). 
 * Each dataset 
 * is sampled homogeneosely (h) or non-homogeneousely (nh) from the same master distribution. 
 * Every dataset is obfuscated by the associated central mechanism into a noisy histogram. 
 * Finally the EMP (or weigthed-EMP) estimator is applied to obtain the overall 
 * estimate by summing the noisy histigrams and normalizing the sum into a 
 * distribution. Results are stored in a file with 
 * 
 */

 /**
 * Original data is always sampled from a master distribution that depends on the space.
 * The datasets are populated homogeneousely (h) or non-homogeneousely depending on the 
 * parameter <sampling_mode>.
 * For instance in the "linear_binomial" space, the master distribution is the binomial 
 * on the linear 100 points with success probaility 0.5. In the planar dataset, the 
 * master distribution is the one computed from the real points in a specific region. 
 * Afterwards, the master distribution is used to construct the original datasets. 
 * */ 

#include <iostream>
#include "FullGridSpace.h"
#include "CentralMechanism.h"
#include "Estimator_EMP.h"
#include "StatDistance.h"
#include "StatUtils.h"

//results directory: updated later depending on the choice of space, mechanisms, and base_size
std::string results_folder; // e.g. "./results/many_databases";

// choose one of these space
// std::string space = "uci_adults_ages";
// std::string space = "linear_binomial"
// std::string space = "manhattan_10km_7km";
// std::string space = "sf_12km_8km";
std::string space_str;

// space variables: set in the main
FullGridSpace xspace(0,0,0);
DoubleRow masterDistribution;

std::string sampling_mode; // takes values "h" for homogeneous or "nh" for non-homogeneous
StatDistanceType statDistanceType; //e.g. EMD, TV, SQD 

// mechanisms
std::vector<CentralMechanism*> mechanism_ptr_list;   // only central: CLAP, CGAUSS

// Client struct is defined in StatUtils.h
std::vector<Client> client_list; // list of the above structure obtained from the argument

std::vector<std::string> estimators;
std::vector<std::ofstream> distance_file_streams;

DoubleMatrix omega; // omega_xs is the probability of assigning data point x to server s

// runs
int runs = 20; // 20 

// forward declarations
void process();

// This program assumes we have multiple datasets, the size of each is 
// dssize. Takes the parameters: dssize, number of datasets (numds), 
// and a list of info_strings for arbitrary (central) mechanisms. 
// The estimator EMP_CLAP is evaluated if the mechanisms are CLaplace
// and the estimator EMP_CGAUSS is evaluated if the mechanisms are CGauss.
int main(int argc, char* argv[]) {
if (argc < 7) {
        std::cerr << "number of arguments = " << argc << "\n";
        std::cerr << "Usage: " << argv[0] 
        << " <space> <sampling_mode> <statDist> "
        << " ""<mechinfo1:dssize1:freq1> ; <mechinfo2:dssize2:freq2> ; ... "" "  // client items
		<< " <estimators> <distanceFiles> " << "\n";

		for (int j=0; j<argc; j++)
			std::cout << "ARG["<<j<<"]= " << argv[j] <<"\n";
		return 1;
    }
	// //To debug:
	// for (int i = 0; i< argc ; i++)
	// 	std::cout << "arg[" << i << "] = " << argv[i] << "\n";
	// exit(0);

	space_str = argv[1];
	sampling_mode = argv[2];
	statDistanceType = StatDistance::stringToDistType(argv[3]);

	// parse arg[4] to obtain the client_item of the form <mechinfo1:dssize1:freq>
	// These tokens are separated by ";"
	// elements of client_item are separated by":" (because an info string may contain spaces)
	// A client item has the information for one client (mech and dssize) and the frequency of this among clients 
	std::istringstream client_list_iss(argv[4]);  
    std::string client_item;
    while (std::getline(client_list_iss, client_item, ';')) {
		// now parse the client_item to extract its elements which are separated by ":"
		std::istringstream client_item_iss(client_item);
		std::string mech_info_str, dssize_str,freq_str;
		if (std::getline(client_item_iss, mech_info_str, ':') &&
            std::getline(client_item_iss, dssize_str, ':') &&
            std::getline(client_item_iss, freq_str, ':')) {
				Client client;
				// get the mechanism index in the mechanism_ptr_list
				size_t i=0; while (i< mechanism_ptr_list.size() && mechanism_ptr_list[i]->get_info() != mech_info_str) i++;
				if (i < mechanism_ptr_list.size())   // a mechanism was already created before
					client.mech_idx = i;
				else {
					// construct every mechanism from given info strings, 
					// e.g. "CLAP 1 100 2.0" and ensure they are central
					Mechanism* mech_ptr = Mechanism::createMechanism(mech_info_str);
					auto central_mech_ptr = dynamic_cast<CentralMechanism*>(mech_ptr);
					if (central_mech_ptr) 
        				mechanism_ptr_list.push_back(central_mech_ptr);
					else {
						std::cerr << " error: [" << mech_info_str << "] is an invalid central mechanism info string";
						delete(mech_ptr); mech_ptr = nullptr;
						exit(1);
					}
					client.mech_idx = mechanism_ptr_list.size()-1;  // the index of the last item in the list
				}
				client.dssize = std::stoi(dssize_str);
				client.freq = std::stoi(freq_str);
            	client_list.push_back(client);
			}
		else {
			std::cerr << "Unable to parse client item: " << client_item;
			exit(1);
		}
	}

	// parse arg[5] to obtain the required estimators separated by spaces
    std::istringstream est_iss(argv[5]);
    std::string est_token;
    while (est_iss >> est_token) {
		estimators.push_back(est_token);	
	}

	// get the distance file for each estimator from argv[6] . 
    // These files are used at every run, in which a new noisy distribution is produced 
	std::istringstream distFiles_iss(argv[6]);
    std::string token;
    while (distFiles_iss >> token) {
		// std::cout << "opening file: " << distFileNames[j] <<"\n";
		distance_file_streams.emplace_back(token, std::ios::app);  // create directly in the vector because ofstream is not copyable 
		if (!distance_file_streams.back().is_open()) {
			std::cerr << "unable to open distance file " << token << "\n";
			exit(1);
		}
	}
	if (distance_file_streams.size() != estimators.size()) {
		std::cerr << "number of given distance files must be equal to the number of estimators\n";
		exit(1);
	}

    // The original data obtained from some dataset. 
	StatUtils::construct_space(space_str, xspace, masterDistribution);  // from utils
	process();
	
	// close files
	for (size_t j=0; j< distance_file_streams.size(); j++) distance_file_streams[j].close();
	// delete mechanisms pointers
	for (size_t j=0; j<mechanism_ptr_list.size(); j++) {
		 delete mechanism_ptr_list[j]; mechanism_ptr_list[j] = nullptr;
	}
}

/** construct datasets by sampling from the master distribution, assign mechanisms from <mechanisms>. 
 * Then obfuscate and estimate for <runs> times*/
void process() {

	// total number of clients and data points
	int numClients = 0; 
	int numPoints = 0;
	for (size_t k=0; k < client_list.size(); k++) {
		numClients = numClients + client_list[k].freq;
		numPoints = numPoints + client_list[k].freq * client_list[k].dssize;
	}

	for (int k=0; k<runs;k++) {
		std::cout << "alphabet size = " << xspace.getSize()
		<< ", dataset size of client0 = " << client_list[0].dssize
		<< ", number of clients = " << numClients  
		<< ", total number of points = " << numPoints  
        << ", run: " << k << "\n\n";

		// original datasets construction
		// The "construct_omega1" uses some threshold in the range of input domain
		// to distribute the points on the clients. So, the marginal probabilities 
		// of the two halfs of clients are not generally equal. 
		// However, the "construct_omega2" uses "masterDistribution" to ensure that 
		// the two halfs of clients have equal marginal probabilities to receive points. 
		//
		// The "sampleOrigData_nh_halfs" sampling function uses the constructed omega 
		// matrix to distribute the points sampled from (masterDistribution) over clients, 
		// but additionally ensures that the two halfs of clients have the same number 
		// of points.
		// 
		std::vector<GridEmpData> origEmpData;
		std::cout << "constructing original datasets (" << sampling_mode << ") \n";
		if (sampling_mode == "h") 
			origEmpData = StatUtils::sampleOriginalData_h(xspace, masterDistribution, client_list);
		else if (sampling_mode == "nh1") {
			omega = StatUtils::construct_omega1(xspace,numClients);
			origEmpData = StatUtils::sampleOriginalData_nh_halfs(xspace, masterDistribution,numPoints,omega);
		}
		else if (sampling_mode == "nh2") { // special sampling that depends on the master distribution
			omega = StatUtils::construct_omega2(xspace,masterDistribution,numClients);
			origEmpData = StatUtils::sampleOriginalData_nh_halfs(xspace, masterDistribution,numPoints,omega);
		}
		else {
			std::cerr << "unrecognized sampling mode";
			exit(1);
		}

        // original distribution (initially filled with zeros)
		DoubleRow original_distribution(xspace.getSize());
		// for (size_t i=0; i< origEmpData.size(); i++)
		// 	original_distribution.addi(xspace.get_histogram_from_emp_data(origEmpData[i]));
		// original_distribution.divi(original_distribution.sum());

		// set the original distribution to the master 
		original_distribution = masterDistribution;
		
		// // to debug: plot the original distribution
		// original_distribution.transpose().write_to_file("original_distribution.txt");
		
		
		// construct a list of noisy histograms and a list of noise_inv_variance
		// each central mechanism obfuscates a client dataset where clients follow
		// the same order of original datasets creation to keep the mapping 
		// between the original dataset and its noisy version.
		int client_idx=0;
		std::vector<DoubleRow> list_noisy_histograms;
		std::vector<double> list_variances; // used for weighted estimates and analytic error evaluation
		std::vector<int> list_subcounts;  // the number of points ate each client
		for (size_t k=0; k<client_list.size(); k++)
			for (int j=0; j< client_list[k].freq; j++) {
				list_noisy_histograms.push_back(mechanism_ptr_list[client_list[k].mech_idx]->obfuscateDataset(origEmpData[client_idx])); 
				list_variances.push_back(mechanism_ptr_list[client_list[k].mech_idx]->getVariance());
				list_subcounts.push_back(client_list[k].dssize);
				client_idx++;
			}


		// estimate:
		// compute the estimated distribution by summing (aggregating) the noisy histograms, 
		// and then normalizing with no -ve correction (O), normalizing with -ve values correction (N),
		// or by normalizing and then projecting to the simplex (P) 
		for (size_t j=0; j<estimators.size(); j++) {
			// plain emp estimate obtainded by aggregating noisy histograms and normalizing
			if (estimators[j]==std::string("EMP")) {
				Estimator_EMP emp(list_noisy_histograms);
				std::cout << "Estimating using EMP method....\n";
				DoubleRow pr_x_emp = emp.getEstimate('O');
        		double emp_distance = StatDistance::getDist(original_distribution, pr_x_emp, statDistanceType, &xspace);
				distance_file_streams[j] << std::to_string(emp_distance) << "\n";
			}
			// normalized emp estimate obtained by aggregating noisy histograms, truncating negative bins, and finally normalizing
			else if (estimators[j]==std::string("EMPN")) {
				Estimator_EMP emp(list_noisy_histograms);
				std::cout << "Estimating using EMPN method....\n";
				DoubleRow pr_x_empn = emp.getEstimate('N');
        		double empn_distance = StatDistance::getDist(original_distribution, pr_x_empn, statDistanceType, &xspace);
				distance_file_streams[j] << std::to_string(empn_distance) << "\n";
			}
			// projected emp estimate obtained by aggregating noisy histograms, normalizing, and finally projecting onto the simplix
			else if (estimators[j]==std::string("EMPP")) {
				Estimator_EMP emp(list_noisy_histograms);
				std::cout << "Estimating using EMPP method....\n";
				DoubleRow pr_x_empp = emp.getEstimate('P');
        		double empp_distance = StatDistance::getDist(original_distribution, pr_x_empp, statDistanceType, &xspace);
				distance_file_streams[j] << std::to_string(empp_distance) << "\n";
			}
			// weighted emp estimate obtained by weighted aggregation of noisy histograms and then normalizing
			else if (estimators[j]==std::string("EMPW")) {
				Estimator_EMP emp(list_noisy_histograms);
				std::cout << "Estimating using EMPW method....\n";
				DoubleRow pr_x_empw = emp.getEstimate(list_variances,'O');
        		double empw_distance = StatDistance::getDist(original_distribution, pr_x_empw, statDistanceType, &xspace);
				distance_file_streams[j] << std::to_string(empw_distance) << "\n";
			}
			// normalized weighted emp estimate obtained by constructing a weighted sum of noisy histograms, truncating negative bins, and finally normalizing
			else if (estimators[j]==std::string("EMPWN")) {
				Estimator_EMP emp(list_noisy_histograms);
				std::cout << "Estimating using EMPWN method....\n";
				DoubleRow pr_x_empwn = emp.getEstimate(list_variances,'N');
        		double empwn_distance = StatDistance::getDist(original_distribution, pr_x_empwn, statDistanceType, &xspace);
				distance_file_streams[j] << std::to_string(empwn_distance) << "\n";
			}
			// projected emp estimate obtained by weighted aggregation of noisy histograms, normalizing, and finally projecting onto the simplix
			else if (estimators[j]==std::string("EMPWP")) {
				Estimator_EMP emp(list_noisy_histograms);
				std::cout << "Estimating using EMPWP method....\n";
				DoubleRow pr_x_empwp = emp.getEstimate(list_variances,'P');
        		double empwp_distance = StatDistance::getDist(original_distribution, pr_x_empwp, statDistanceType, &xspace);
				distance_file_streams[j] << std::to_string(empwp_distance) << "\n";
			}
			// pseudo emp estimates obtained by aggregating the noisy histograms and dividing by the original number of data points.
			// Note it is not a distribution.
			else if (estimators[j]==std::string("PSEMP")) {
				Estimator_EMP emp(list_noisy_histograms);
				std::cout << "Estimating using PseudoEMP method....\n";
				DoubleRow pr_x_psemp = emp.getPseudoEstimate(numPoints);
        		double psemp_distance = StatDistance::getDist(original_distribution, pr_x_psemp, statDistanceType, &xspace);
				distance_file_streams[j] << std::to_string(psemp_distance) << "\n";
			}
			// pseudo emp estimates obtained by weighted aggregation of the noisy histograms and dividing by the weighted sum of subcounts.
			// Not necessarily a distribution
			else if (estimators[j]==std::string("PSEMPW")) {
				Estimator_EMP emp(list_noisy_histograms);
				std::cout << "Estimating using weighted Pseudo-EMP method....\n";
				DoubleRow pr_x_psempw = emp.getPseudoEstimate(list_variances,list_subcounts);
        		double psempw_distance = StatDistance::getDist(original_distribution, pr_x_psempw, statDistanceType, &xspace);
				distance_file_streams[j] << std::to_string(psempw_distance) << "\n";
			}

			// Analytic errors (only for SQD)

			// Analytic error of the pseudo emp estimates evaluated using the theoritical error expression.
			else if (estimators[j]==std::string("APSEMP")) {
				if (statDistanceType != SQD) {
					std::cerr << "Unable to evaluate the analytic error: statistical distance must be SQD";
					exit(1);
				}
				std::cout << "Evaluating analytic error of the PseudoEMP method....\n";
        		double apsemp_distance = Estimator_EMP::getPseudoEstimateAnalyticMeanSQD(original_distribution, list_variances, numPoints);
				distance_file_streams[j] << std::to_string(apsemp_distance) << "\n";
			}
			// pseudo emp estimates obtained by weighted aggregation of the noisy histograms and dividing by the weighted sum of subcounts.
			// Not necessarily a distribution
			else if (estimators[j]==std::string("APSEMPW")) {
				if (statDistanceType != SQD) {
					std::cerr << "Unable to evaluate the analytic error: statistical distance must be SQD";
					exit(1);
				}
				std::cout << "Evaluating the analytic error of the weighted Pseudo-EMP method....\n";
        		// for (int j=0; j< list_variances.size(); j++) std::cout << "p_list_variances["<<j<<"] = " << list_variances[j] << "\n";
        		double apsempw_distance = Estimator_EMP::getPseudoEstimateAnalyticMeanSQD(original_distribution, list_variances, list_subcounts);
				distance_file_streams[j] << std::to_string(apsempw_distance) << "\n";
			}
			else {
				std::cerr << "Estimator " << estimators[j] << " is not recognized";
				exit(1); 
			}
		}
	}
}
