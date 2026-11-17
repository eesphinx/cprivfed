/**
 * runs one time and saves the master distribution, original histograms, noisy_histogram, and estimated distribution using GIBU 
 * */ 
 

#include <iostream>
#include "AbstractOutputLocalMechanism.h"
#include "GridEmpData.h"
#include "Estimator_COMBINED_RAPPOR.h"
#include "Estimator_AVG_MECH_INV.h"
#include "Estimator_AVG_MECH_IBU.h"
#include "Estimator_AVG_RESULT_RAPPOR.h"
#include "Estimator_AVG_RESULT_INV.h"
#include "Estimator_AVG_RESULT_IBU.h"
#include "Estimator_GIBU.h"
#include "StatUtils.h"


// Choose the space, i.e. the dataset:
std::string space_str;
// may take one of the following values:
// std::string space = "uci_adults_ages";
// std::string space = "linear_binomial"
// std::string space = "manhattan_10km_7km";
// std::string space = "sf_12km_8km";

// space variables: set in the main
FullGridSpace xspace(0,0,0);
DoubleRow masterDistribution;

std::string sampling_mode; // takes values "h" for homogeneous or "nh" for non-homogeneous
StatDistanceType statDistanceType; //e.g. EMD, TV, SQD

// mechanisms
std::vector<AbstractOutputLocalMechanism*> mechanism_ptr_list;   // only local

// Client struct is defined in StatUtils.h
std::vector<Client> client_list; // list of the above structure obtained from the argument

std::vector<std::string> estimators;
std::vector<std::ofstream> distance_file_streams;  // must has size equal to estimators
std::vector<std::ofstream> time_file_streams; // must has size equal to estimators


DoubleMatrix omega; // omega_xs is the probability of assigning data point x to server s
// ibu/gibu parameters
int max_iterations = 1e6; // max iterations (should be large enough to reach the termination error)
double error = 1e-5;    //default: 1e-5;  // the TV error at which the em terminates (use smaller e.g. 1e-5 for faster execution)

// runs
int runs = 1; // make it smaller for faster results

bool recordTimes=false;  // set if required

// forward declarations
void process();

// This program assumes we have <numds> datasets, each populated by <dssize> points. 
int main(int argc, char* argv[]) {
	if (argc < 7) {
        std::cerr << "number of arguments = " << argc << "\n";
        std::cerr << "Usage: " << argv[0] 
        << " <space> <sampling_mode> <statDist> "
        << " ""<mechinfo1:dssize1:freq1> ; <mechinfo2:dssize2:freq2> ; ... "" "  // client items
		<< " <estimators> <distanceFiles> [optinal <timeFiles>] " << "\n";

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
					// e.g. "KRR 1 100 10.0" and ensure they are local
					Mechanism* mech_ptr = Mechanism::createMechanism(mech_info_str);
					AbstractOutputLocalMechanism* local_mech_ptr = dynamic_cast<AbstractOutputLocalMechanism*>(mech_ptr);
					if (local_mech_ptr) 
        				mechanism_ptr_list.push_back(local_mech_ptr);
					else {
						std::cerr << " error: [" << mech_info_str << "] is not a local mechanism info string";
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

	// //debug:
	// for (int i=0; i<distFileNames.size(); i++) {
	// 	std::cout << "item "<< i << " is: " << distFileNames[i] <<"\n";
	// }
	// exit(0);

    // optional execution time files
	if (argc==8) {
		std::istringstream timeFiles_iss(argv[7]);
   		std::string token;
    	while (timeFiles_iss >> token) {
			// std::cout << "opening time file: " << token <<"\n";
			time_file_streams.emplace_back(token, std::ios::app);  // create directly in the vector because ofstream is not copyable 
			if (!time_file_streams.back().is_open()) {
				std::cerr << "unable to open time file " << token << "\n";
				exit(1);
			}
		}
		if (time_file_streams.size() != estimators.size()) {
			std::cerr << "number of given time files must be equal to the number of estimators\n";
			exit(1);
		}
		recordTimes=true;
	}
	
    // The original data obtained from some dataset. 
	StatUtils::construct_space(space_str, xspace, masterDistribution);  // from utils
	process();
	
	// close files
	for (size_t j=0; j< distance_file_streams.size(); j++) distance_file_streams[j].close();
	for (size_t j=0; j< time_file_streams.size(); j++) time_file_streams[j].close();
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
			std::cerr << "unrecognized sampling mode : " << sampling_mode;
			exit(1);
		}

		
		// original distribution (initially filled with zeros)
		DoubleRow original_distribution(xspace.getSize());
		// for (size_t i=0; i< origEmpData.size(); i++)
		// 	original_distribution.addi(xspace.get_histogram_from_emp_data(origEmpData[i]));
		// std::cout << "overall distribution constructed from " << original_distribution.sum() << " data points\n";
		// original_distribution.divi(original_distribution.sum());

		// set the original distribution to the master 
		original_distribution = masterDistribution;
		
		// to debug: write the original distribution
		original_distribution.transpose().write_to_file("original_distribution.txt");

		// // Writing original distribution and histograms: 
		// write the overall original histograms, 
		// and assuming we have two groups, write the original histogram/distribution for the 1st group of servers
		DoubleRow original_histogram(xspace.getSize());
		for (size_t i=0; i< origEmpData.size(); i++)
			original_histogram.addi(xspace.get_histogram_from_emp_data(origEmpData[i]));
		original_histogram.transpose().write_to_file("original_histogram.txt");

		DoubleRow original_histogram0(xspace.getSize());
		for (size_t i=0; i< origEmpData.size()/2; i++)
			original_histogram0.addi(xspace.get_histogram_from_emp_data(origEmpData[i]));
		// std::cout << "overall group1 histogram constructed from " << original_histogram0.sum() << " data points\n";
		original_histogram0.transpose().write_to_file("original_histogram0.txt");

		DoubleRow original_histogram1(xspace.getSize());
		for (size_t i=origEmpData.size()/2; i< origEmpData.size(); i++)
			original_histogram1.addi(xspace.get_histogram_from_emp_data(origEmpData[i]));
		// std::cout << "overall group2 histogram constructed from " << original_histogram1.sum() << " data points\n";
		original_histogram1.transpose().write_to_file("original_histogram1.txt");

		
		// construct noisy empdata for each client: 
		// we need to follow the same order of constructing origEmpData from from client_list 
		// This is important to keep the mapping between original dataset and its noisy version
		int client_idx=0;
		std::vector<EmpData> list_noisy_data;
		for (size_t k=0; k<client_list.size(); k++)
			for (int j=0; j< client_list[k].freq; j++) {
				list_noisy_data.push_back(EmpData(origEmpData[client_idx], mechanism_ptr_list[client_list[k].mech_idx])); 
				client_idx++;
			}

		// Write the noisy distribution
		DoubleRow noisy_distribution(xspace.getSize());
		for (size_t i=0; i< list_noisy_data.size(); i++)
			noisy_distribution.addi(xspace.get_histogram_from_emp_data(list_noisy_data[i]));
		noisy_distribution.divi(noisy_distribution.sum());
		noisy_distribution.transpose().write_to_file("noisy_distribution.txt");


		// estimate
		for (size_t j=0; j<estimators.size(); j++) {
			if (estimators[j]==std::string("COMBINED_RAPPORN")) {
            	auto start = std::chrono::high_resolution_clock::now();		
				Estimator_COMBINED_RAPPOR combined_rappor(xspace, list_noisy_data);
				std::cout << "Estimating using COMBINED_RAPPORN method....\n";
				DoubleRow pr_x_combined_rapporn = combined_rappor.getMaximizer('N');
            	auto end = std::chrono::high_resolution_clock::now();
            	std::chrono::duration<double> elapsed = end - start;
        		double combined_rapporn_distance = StatDistance::getDist(original_distribution, pr_x_combined_rapporn, statDistanceType, &xspace);
				distance_file_streams[j] << std::to_string(combined_rapporn_distance) << "\n";
            	if (recordTimes)
					time_file_streams[j] << elapsed.count() << "\n"; 
			}
			else if (estimators[j]==std::string("COMBINED_RAPPORP")) {
            	auto start = std::chrono::high_resolution_clock::now();		
				Estimator_COMBINED_RAPPOR combined_rappor(xspace, list_noisy_data);
				std::cout << "Estimating using COMBINED_RAPPORP method....\n";
				DoubleRow pr_x_combined_rapporp = combined_rappor.getMaximizer('P');
            	auto end = std::chrono::high_resolution_clock::now();
            	std::chrono::duration<double> elapsed = end - start;
        		double combined_rapporp_distance = StatDistance::getDist(original_distribution, pr_x_combined_rapporp, statDistanceType, &xspace);
				distance_file_streams[j] << std::to_string(combined_rapporp_distance) << "\n";
            	if (recordTimes)
					time_file_streams[j] << elapsed.count() << "\n"; 
			}
		
			else if (estimators[j]==std::string("AVG_MECH_INV")) {
            	auto start = std::chrono::high_resolution_clock::now();		
				Estimator_AVG_MECH_INV avg_mech_inv(xspace, list_noisy_data);
				std::cout << "Estimating using AVG_MECH_INV method....\n";
				DoubleRow pr_x_avg_mech_inv = avg_mech_inv.getMaximizer('O');
            	auto end = std::chrono::high_resolution_clock::now();
            	std::chrono::duration<double> elapsed = end - start;
        		double avg_mech_inv_distance = StatDistance::getDist(original_distribution, pr_x_avg_mech_inv, statDistanceType, &xspace);
				distance_file_streams[j] << std::to_string(avg_mech_inv_distance) << "\n";
            	if (recordTimes)
					time_file_streams[j] << elapsed.count() << "\n"; 
			}
			else if (estimators[j]==std::string("AVG_MECH_INVN")) {
            	auto start = std::chrono::high_resolution_clock::now();		
				Estimator_AVG_MECH_INV avg_mech_inv(xspace, list_noisy_data);
				std::cout << "Estimating using AVG_MECH_INVN method....\n";
				DoubleRow pr_x_avg_mech_invn = avg_mech_inv.getMaximizer('N');
            	auto end = std::chrono::high_resolution_clock::now();
            	std::chrono::duration<double> elapsed = end - start;
        		double avg_mech_invn_distance = StatDistance::getDist(original_distribution, pr_x_avg_mech_invn, statDistanceType, &xspace);
				distance_file_streams[j] << std::to_string(avg_mech_invn_distance) << "\n";
            	if (recordTimes)
					time_file_streams[j] << elapsed.count() << "\n"; 
			}
			else if (estimators[j]==std::string("AVG_MECH_INVP")) {
            	auto start = std::chrono::high_resolution_clock::now();		
				Estimator_AVG_MECH_INV avg_mech_inv(xspace, list_noisy_data);
				std::cout << "Estimating using AVG_MECH_INVP method....\n";
				DoubleRow pr_x_avg_mech_invp = avg_mech_inv.getMaximizer('P');
            	auto end = std::chrono::high_resolution_clock::now();
            	std::chrono::duration<double> elapsed = end - start;
        		double avg_mech_invp_distance = StatDistance::getDist(original_distribution, pr_x_avg_mech_invp, statDistanceType, &xspace);
				distance_file_streams[j] << std::to_string(avg_mech_invp_distance) << "\n";
            	if (recordTimes)
					time_file_streams[j] << elapsed.count() << "\n"; 
			}
		
			else if (estimators[j]==std::string("AVG_MECH_IBU")) {
            	auto start = std::chrono::high_resolution_clock::now();		
				Estimator_AVG_MECH_IBU avg_mech_ibu(xspace, list_noisy_data);
				std::cout << "Estimating using AVG_MECH_IBU method....\n";
				DoubleRow pr_x_avg_mech_ibu = avg_mech_ibu.getMaximizer(error, max_iterations);
            	auto end = std::chrono::high_resolution_clock::now();
            	std::chrono::duration<double> elapsed = end - start;
				double avg_mech_ibu_distance = StatDistance::getDist(original_distribution, pr_x_avg_mech_ibu, statDistanceType, &xspace);
				distance_file_streams[j] << std::to_string(avg_mech_ibu_distance) << "\n";
            	if (recordTimes)
					time_file_streams[j] << elapsed.count() << "\n"; 
			}

			else if (estimators[j]==std::string("AVG_RESULT_RAPPORN")) {
            	auto start = std::chrono::high_resolution_clock::now();		
				Estimator_AVG_RESULT_RAPPOR avg_result_rappor(xspace, list_noisy_data);
				std::cout << "Estimating using AVG_RESULT_RAPPORN method....\n";
				DoubleRow pr_x_avg_result_rapporn = avg_result_rappor.getMaximizer('N');
            	auto end = std::chrono::high_resolution_clock::now();
            	std::chrono::duration<double> elapsed = end - start;
				double avg_result_rapporn_distance = StatDistance::getDist(original_distribution, pr_x_avg_result_rapporn, statDistanceType, &xspace);
				distance_file_streams[j] << std::to_string(avg_result_rapporn_distance) << "\n";
            	if (recordTimes)
					time_file_streams[j] << elapsed.count() << "\n"; 
			}
			else if (estimators[j]==std::string("AVG_RESULT_RAPPORP")) {
            	auto start = std::chrono::high_resolution_clock::now();		
				Estimator_AVG_RESULT_RAPPOR avg_result_rappor(xspace, list_noisy_data);
				std::cout << "Estimating using AVG_RESULT_RAPPORP method....\n";
				DoubleRow pr_x_avg_result_rapporp = avg_result_rappor.getMaximizer('P');
            	auto end = std::chrono::high_resolution_clock::now();
            	std::chrono::duration<double> elapsed = end - start;
				double avg_result_rapporp_distance = StatDistance::getDist(original_distribution, pr_x_avg_result_rapporp, statDistanceType, &xspace);
				distance_file_streams[j] << std::to_string(avg_result_rapporp_distance) << "\n";
            	if (recordTimes)
					time_file_streams[j] << elapsed.count() << "\n"; 
			}
			
			else if (estimators[j]==std::string("AVG_RESULT_INVN")) {
            	auto start = std::chrono::high_resolution_clock::now();		
				Estimator_AVG_RESULT_INV avg_result_inv(xspace, list_noisy_data);
				std::cout << "Estimating using AVG_RESULT_INVN method....\n";
				DoubleRow pr_x_avg_result_invn = avg_result_inv.getMaximizer('N');
            	auto end = std::chrono::high_resolution_clock::now();
            	std::chrono::duration<double> elapsed = end - start;
				double avg_result_invn_distance = StatDistance::getDist(original_distribution, pr_x_avg_result_invn, statDistanceType, &xspace);
				distance_file_streams[j] << std::to_string(avg_result_invn_distance) << "\n";
            	if (recordTimes)
					time_file_streams[j] << elapsed.count() << "\n"; 
			}
			else if (estimators[j]==std::string("AVG_RESULT_INVP")) {
            	auto start = std::chrono::high_resolution_clock::now();		
				Estimator_AVG_RESULT_INV avg_result_inv(xspace, list_noisy_data);
				std::cout << "Estimating using AVG_RESULT_INVP method....\n";
				DoubleRow pr_x_avg_result_invp = avg_result_inv.getMaximizer('P');
            	auto end = std::chrono::high_resolution_clock::now();
            	std::chrono::duration<double> elapsed = end - start;
				double avg_result_invp_distance = StatDistance::getDist(original_distribution, pr_x_avg_result_invp, statDistanceType, &xspace);
				distance_file_streams[j] << std::to_string(avg_result_invp_distance) << "\n";
            	if (recordTimes)
					time_file_streams[j] << elapsed.count() << "\n"; 
			}

			else if (estimators[j]==std::string("AVG_RESULT_IBU")) {
            	auto start = std::chrono::high_resolution_clock::now();		
				Estimator_AVG_RESULT_IBU avg_result_ibu(xspace, list_noisy_data);
				std::cout << "Estimating using AVG_RESULT_IBU method....\n";
				DoubleRow pr_x_avg_result_ibu = avg_result_ibu.getMaximizer(error, max_iterations);
            	auto end = std::chrono::high_resolution_clock::now();
            	std::chrono::duration<double> elapsed = end - start;
				double avg_result_ibu_distance = StatDistance::getDist(original_distribution, pr_x_avg_result_ibu, statDistanceType, &xspace);
				distance_file_streams[j] << std::to_string(avg_result_ibu_distance) << "\n";
            	if (recordTimes)
					time_file_streams[j] << elapsed.count() << "\n"; 
			}
	
			else if (estimators[j]==std::string("GIBU")) {
				DoubleRow pr_x_gibu;
				if (omega.getLength()==0) {    // if no omega is specified
					Estimator_GIBU gibu(xspace, list_noisy_data);  
					std::cout << "Estimating using GIBU method with no omega matrix....\n";
					auto start = std::chrono::high_resolution_clock::now();		
					pr_x_gibu = gibu.getMaximizer(error, max_iterations);
            		auto end = std::chrono::high_resolution_clock::now();
            		std::chrono::duration<double> elapsed = end - start;
					double distance_gibu = StatDistance::getDist(original_distribution, pr_x_gibu, statDistanceType, &xspace);
					distance_file_streams[j] << std::to_string(distance_gibu) << "\n";
            		if (recordTimes)
					time_file_streams[j] << elapsed.count() << "\n"; 
				} // The only difference is to use omega
     			else {
					// std::cout << "omega rows = " << omega.getRows() << ", columns = " << omega.getColumns() << "\n";
					// for (int i=0; i< omega.getRows(); i++) {
					// 	std::cout << "sum of row " << i << " = " << omega.selectRow(i).sum() << "\n";
					// }
					// omega.write_to_file("omega.txt");
					Estimator_GIBU gibu(xspace, list_noisy_data, omega);  
					std::cout << "Estimating using GIBU method with omega matrix....\n";
					auto start = std::chrono::high_resolution_clock::now();		
					pr_x_gibu = gibu.getMaximizer(error, max_iterations);
            		auto end = std::chrono::high_resolution_clock::now();
            		std::chrono::duration<double> elapsed = end - start;
					double distance_gibu = StatDistance::getDist(original_distribution, pr_x_gibu, statDistanceType, &xspace);
					distance_file_streams[j] << std::to_string(distance_gibu) << "\n";
            		if (recordTimes)
					time_file_streams[j] << elapsed.count() << "\n"; 
				}
			// Write the GIBU estimate: 
			pr_x_gibu.transpose().write_to_file("pr_x_gibu.txt");
			}
			else {
				std::cerr << "Estimator " << estimators[j] << " is not recognized";
				exit(1); 
			}
		}

	}

}



// some old settings for kRR to produce plots of the email: 

    // fixed eps = 1.0:
    //	eps[0] = 1.0;  	eps[1] = 1.0;   eps[2] = 1.0;  	eps[3] = 1.0;     eps[4] = 1.0;
    //	eps[5] = 1.0;  	eps[6] = 1.0;   eps[7] = 1.0;  	eps[8] = 1.0;     eps[9] = 1.0;


    // fixed eps = 2.0:
    //	eps[0] = 2.0;  	eps[1] = 2.0;   eps[2] = 2.0;  	eps[3] = 2.0;     eps[4] = 2.0;
    //	eps[5] = 2.0;  	eps[6] = 2.0;   eps[7] = 2.0;  	eps[8] = 2.0;     eps[9] = 2.0;

    // fixed eps = 1.0:
    //  eps[0] = 10.0;  	eps[1] = 10.0;   eps[2] = 10.0;  	eps[3] = 10.0;     eps[4] = 10.0;
    //	eps[5] = 10.0;  	eps[6] = 10.0;   eps[7] = 10.0;  	eps[8] = 10.0;     eps[9] = 10.0;


    // variable eps between 2.0 and 20.0
    //	eps[0] = 2.0;  	    eps[1] = 4.0;    eps[2] = 6.0;  	eps[3] = 8.0;      eps[4] = 10.0;
    //	eps[5] = 12.0;  	eps[6] = 14.0;   eps[7] = 16.0;  	eps[8] = 18.0;     eps[9] = 20.0;

    // variable eps between 1.0 and 10.0 (default configuration)
    //	eps[0] = 1.0;  	eps[1] = 2.0;   eps[2] = 3.0;  	eps[3] = 4.0;     eps[4] = 5.0;
    //	eps[5] = 6.0;  	eps[6] = 7.0;   eps[7] = 8.0;  	eps[8] = 9.0;     eps[9] = 10.0;

    // variable eps between 1.0 and 10.0 (configuration 2) which shows a bigger difference between methods
    //	eps[0] = 1.0;  	eps[1] = 1.0;   eps[2] = 1.0;  	eps[3] = 1.0;     eps[4] = 1.0;
    //	eps[5] = 5.0;  	eps[6] = 5.0;   eps[7] = 10.0;  eps[8] = 10.0;    eps[9] = 10.0;



    // variable eps between 0.2 and 10.0 (configuration 3) which all are amplified by 400-size datasets.
    //	eps[0] = 0.2;  	eps[1] = 0.2;   eps[2] = 0.4;  	eps[3] = 0.4;     eps[4] = 0.6;
    //	eps[5] = 0.6;  	eps[6] = 10.0;   eps[7] = 10.0;   eps[8] = 10.0;     eps[9] = 10.0;

