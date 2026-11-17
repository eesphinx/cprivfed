#include <iostream>
#include <sstream>
#include <fstream>
#include <regex>
#include <cmath>
#include <algorithm>
#include "DoubleRow.h"
#include "GDataset.h"
#include "GridEmpData.h"
#include "UCIAgeEmpData.h"
#include "UCILabTestsEmpData.h"
#include "FullGridSpace.h"
#include "Distributions.h"
#include "StatUtils.h"

/** 
 * 
 * This source includes utility functions such as constructing the space and sampling data
 * 
 **/ 


// construction of the xspace and master_distribution using the given space_str
void StatUtils::construct_space(const std::string space_str, FullGridSpace &p_xspace, DoubleRow &p_masterDistribution) {
    std::smatch match; // used for matching patterns

	// special spaces with variable domain size, used for non-homogeneous sampling, with generic size
	// case [1] : inverted discrete Gauss with large variance: Catuscia's choice 
	if (std::regex_match(space_str, match, std::regex(R"(linear_InvDiscGaussA(\d+))") )) {  
		// the regex defines two matching groups. match[1] corresponds to (\d+), meaning 1 or more digits  
        // match[1] contains the number part
		p_xspace.cell_side = 1.0;
		p_xspace.hsize=1;   // to force the grid to be one dimensional
		p_xspace.vsize = std::stoi(match[1].str());
        p_masterDistribution = Distributions::getDiscGaussDistribution(p_xspace.vsize-1, (double) p_xspace.vsize/2, (double) p_xspace.vsize/3);
		p_masterDistribution.subi(1.5*p_masterDistribution.max()).mul(-1);
		p_masterDistribution.divi(p_masterDistribution.sum()); 
	}

	// case [2] : discrete Beta 
	if (std::regex_match(space_str, match, std::regex(R"(linear_DiscBeta(\d+))") )) {  
		p_xspace.cell_side = 1.0;
		p_xspace.hsize=1;   // to force the grid to be one dimensional
		p_xspace.vsize = std::stoi(match[1].str());
		p_masterDistribution = Distributions::getDiscBetaDistribution(p_xspace.vsize-1, 5.0, 5.0);
	}

	// case [3] : inverted discrete gauss with small variance 
	if (std::regex_match(space_str, match, std::regex(R"(linear_InvDiscGaussB(\d+))") )) {  
		p_xspace.cell_side = 1.0;
		p_xspace.hsize=1;   // to force the grid to be one dimensional
		p_xspace.vsize = std::stoi(match[1].str());
		p_masterDistribution = Distributions::getDiscGaussDistribution(p_xspace.vsize-1, (double) p_xspace.vsize/2, (double) p_xspace.vsize/6);
		p_masterDistribution.subi(1.1*p_masterDistribution.max()).mul(-1);
		p_masterDistribution.divi(p_masterDistribution.sum()); 
	}

	// case [4] : bimodal distribution (combination of two binomials)
	if (std::regex_match(space_str, match, std::regex(R"(linear_Bimodal(\d+))") )) {  
		p_xspace.cell_side = 1.0;
		p_xspace.hsize=1;   // to force the grid to be one dimensional
		p_xspace.vsize = std::stoi(match[1].str());
		p_masterDistribution = Distributions::getBinomialDistribution(p_xspace.vsize-1, 0.3).muli(0.5);
		p_masterDistribution.addi(Distributions::getBinomialDistribution(p_xspace.vsize-1, 0.7).muli(0.5));
	}

	// generic size binomial
	if (std::regex_match(space_str, match, std::regex(R"(linear_binomial(\d+))") )) {
		p_xspace.cell_side = 1.0;
		p_xspace.hsize=1;   // to force the grid to be one dimensional
		p_xspace.vsize = std::stoi(match[1].str());
		p_masterDistribution = Distributions::getBinomialDistribution(p_xspace.vsize-1, 0.5);
    } 
	
	// generic size discrete gauss
	if (std::regex_match(space_str, match, std::regex(R"(DiscGauss(\d+))") )) {
		p_xspace.cell_side = 1.0;
		p_xspace.hsize=1;   // to force the grid to be one dimensional
		p_xspace.vsize = std::stoi(match[1].str());
		p_masterDistribution = Distributions::getDiscGaussDistribution(p_xspace.vsize-1, (double) p_xspace.vsize/2, (double) p_xspace.vsize/16);
    } 

	// special space used for Rappor in IBU paper
	if (space_str == "linear_uniform10") {
		p_xspace.cell_side =1.0;
		p_xspace.hsize=1;   // to force the grid to be one dimensional
		p_xspace.vsize = 10;
		GridEmpData masterEmpData = sampleUniformData(p_xspace, 100000, 0.3*p_xspace.vsize, 0.6*p_xspace.vsize);
		p_masterDistribution = p_xspace.get_distribution_from_emp_data(masterEmpData);
	}
	
	// // Other possible distributions

	// // inverted discrete beta
	// masterDistribution = Distributions::getDiscBetaDistribution(xspace.vsize-1, 5.0, 5.0);
	// masterDistribution.subi(masterDistribution.max()).mul(-1);
	// masterDistribution.divi(masterDistribution.sum()); 

	// // inverted binomial
	// masterDistribution = Distributions::getBinomialDistribution(xspace.vsize-1, 0.5).sub(1).mul(-1);
	// masterDistribution.divi(masterDistribution.sum()); 

	// // Discrete Gauss
	// masterDistribution = Distributions::getDiscGaussDistribution(xspace.vsize-1, (double) xspace.vsize/2, (double) xspace.vsize/7);

	// to debug:
	// masterDistribution.transpose().write_to_file("master.txt");

	// real datasets
	if (space_str == "uci_adults_ages") {
		// Ages space specification
		int min_age = 15; int max_age = 95; int age_bin_size = 1;  // settings for k-RR, gemetric, and Laplace mechanisms
		//int min_age = 0; int max_age = 99; int age_bin_size = 10;  // a setting for rappor
		// grid space specification
		p_xspace.cell_side = age_bin_size;
		p_xspace.hsize=1;   // to force the grid to be one dimensional
		// In the following (max_age-min_age+1) is the number of different ages between min,max
		// Casted to Double, and divided by age_bin_size to get the number of bins (with fractions).
		// We take the "ceil" because if the number of bins has fractions, this fraction should be considered a final bin.
		p_xspace.vsize = (int) ceil( (double)(max_age-min_age+1)/age_bin_size );
		UCIAgeEmpData age_empdata("../../datasets/adult.data",min_age,max_age,age_bin_size);
		p_masterDistribution = p_xspace.get_distribution_from_emp_data(age_empdata);
	}

	if (space_str == "uci_diabetes_labTests") {
		// LabTests space specification
		int min_labTests = 1; int max_labTests = 140; int labTests_bin_size = 1;  // settings for k-RR, gemetric, and Laplace mechanisms
		// grid space specification
		p_xspace.cell_side = labTests_bin_size;
		p_xspace.hsize=1;   // to force the grid to be one dimensional
		// In the following (max_-min_labTests+1) is the number of different labTests between min,max
		// Casted to Double, and divided by labTests_bin_size to get the number of bins (with fractions).
		// We take the "ceil" because if the number of bins has fractions, this fraction should be considered a final bin.
		p_xspace.vsize = (int) ceil( (double)(max_labTests-min_labTests+1)/labTests_bin_size );
		UCILabTestsEmpData labTests_empdata("./datasets/diabetes_data.csv",min_labTests,max_labTests,labTests_bin_size);
		p_masterDistribution = p_xspace.get_distribution_from_emp_data(labTests_empdata);
		// debug:
		// masterDistribution.transpose().write_to_file("original.txt");
		// exit(0);
	}

	if (space_str== "sf_12km_8km") {
		// San Francisco (used in EuroS&P)
		// containing 123273 checkins for 5674 users
		// east-to-west: 12.0 km, north-to-south: 8.0 km, 
		// Grid: hsize = 24 , vsize = 16, cell side :0.5 km
		double south = 37.7228;
		double north = 37.7946;
		double west = -122.5153;
		double east = -122.3789;
		// grid space specification
		p_xspace.cell_side = 0.5;  // 0.5 for KRR,GEOM  
		std::cout << "\nConstructing dataset...\n";
		std::string datasetFile("./datasets/Gowalla_totalCheckins.txt");
    	GDataset ds(datasetFile, south, north, west, east, p_xspace.cell_side, std::regex(".*"));
		// ds.write_data_to_file("sf_12km_8km.txt");
		//std::cout << ds.get_info() << "\n";
		//std::cout << "The data contains " << ds.get_number_of_distinct_checkins() <<" distinct checkins.\n";
    	p_xspace.hsize = ds.get_hsize(p_xspace.cell_side);
    	p_xspace.vsize = ds.get_vsize(p_xspace.cell_side);
		GridEmpData masterEmpData = ds.get_grid_emp_data(p_xspace.cell_side);   //c++ handles this by default assignment operator which copies the members
    	p_masterDistribution = p_xspace.get_distribution_from_emp_data(masterEmpData);
	}

	// This is the same space of sf_12km_8km except that the cell_side = 1.0 instead of 0.5
	// to allow experiments on Shokri's mechanism
	if (space_str== "sf1_12km_8km") {
		// San Francisco (used in EuroS&P)
		// containing 123273 checkins for 5674 users
		// east-to-west: 12.0 km, north-to-south: 8.0 km, 
		// Grid: hsize = 12 , vsize = 8, cell side :1.0 km
		double south = 37.7228;
		double north = 37.7946;
		double west = -122.5153;
		double east = -122.3789;
		// grid space specification
		p_xspace.cell_side = 1.0;  // 1.0 for Shokri's mechanisms
		std::cout << "\nConstructing dataset...\n";
		std::string datasetFile("./datasets/Gowalla_totalCheckins.txt");
    	GDataset ds(datasetFile, south, north, west, east, p_xspace.cell_side, std::regex(".*"));
		//std::cout << ds.get_info() << "\n";
		//std::cout << "The data contains " << ds.get_number_of_distinct_checkins() <<" distinct checkins.\n";
    	p_xspace.hsize = ds.get_hsize(p_xspace.cell_side);
    	p_xspace.vsize = ds.get_vsize(p_xspace.cell_side);
		GridEmpData masterEmpData = ds.get_grid_emp_data(p_xspace.cell_side);   //c++ handles this by default assignment operator which copies the members
    	p_masterDistribution = p_xspace.get_distribution_from_emp_data(masterEmpData);
	}

	if (space_str=="manhattan_10km_7km") {
		// planar space specification
		// Manhattan (used for journal paper)
		// south-to-north: 10km, west-to-est:7km 
		// # Number of checkins = 98060
		// Total number of users in the region = 5737
		double south= 40.7044;
		double north= 40.7942;
		double west= -74.0205;
		double east= -73.9375;
		// grid space specification
		p_xspace.cell_side = 0.5;
		std::cout << "\nConstructing dataset...\n";
		std::string datasetFile("./datasets/Gowalla_totalCheckins.txt");
    	GDataset ds(datasetFile, south, north, west, east, p_xspace.cell_side, std::regex(".*"));
		// std::cout << ds.get_info() << "\n";
		// std::cout << "The data contains " << ds.get_number_of_distinct_checkins() <<" distinct checkins.\n";
    	p_xspace.hsize = ds.get_hsize(p_xspace.cell_side);
    	p_xspace.vsize = ds.get_vsize(p_xspace.cell_side);
		GridEmpData masterEmpData = ds.get_grid_emp_data(p_xspace.cell_side);   //c++ handles this by default assignment operator which copies the members
    	p_masterDistribution = p_xspace.get_distribution_from_emp_data(masterEmpData);
	}

	if (space_str=="paris_10km_12km") {
		// planar space specification
		// Paris (used for journal paper)
		// south-to-north: 10km, west-to-est:12km 
		// Grid: hsize = 20 , vsize = 24
		// cell side :0.5 km
		// checkins: 13979	    
		double south = 48.8176;
		double north = 48.9074;
		double west = 2.2509;
		double east = 2.4149;

		//The side length in km of one cell in the grid. Used to partition the region in square cells. 
		p_xspace.cell_side = 0.5;   //in km 	   

		std::cout << "\nConstructing dataset...\n";
		std::string datasetFile("./datasets/Gowalla_totalCheckins.txt");
    	GDataset ds(datasetFile, south, north, west, east, p_xspace.cell_side, std::regex(".*"));

		// to debug:
		// ds.write_map_script("paris.html","MARKERS",true);
		// exit(0);

		// std::cout << ds.get_info() << "\n";
		// std::cout << "The data contains " << ds.get_number_of_distinct_checkins() <<" distinct checkins.\n";
    	p_xspace.hsize = ds.get_hsize(p_xspace.cell_side);
    	p_xspace.vsize = ds.get_vsize(p_xspace.cell_side);
		GridEmpData masterEmpData = ds.get_grid_emp_data(p_xspace.cell_side);   //c++ handles this by default assignment operator which copies the members
    	p_masterDistribution = p_xspace.get_distribution_from_emp_data(masterEmpData);
	}

	// generic size linear_binomial with flexible parameter used in the program: ibu_with_multiple_settings.cpp
    std::smatch match_flex_binomial;
	if (std::regex_match(space_str, match, std::regex(R"(linear_space(\d+)_BINOMIAL(0\.\d|1\.0))") )) {
		//the regex defines two matching groups:
        // match[1], corresponding to (\d+), contains the first number (size), and match[2] contains the second (binomial parameter)
		p_xspace.cell_side = 1.0;
		p_xspace.hsize=1;   // to force the grid to be one dimensional
		p_xspace.vsize = std::stoi(match[1].str());
		double  binomial_param  = std::stod(match[2].str());
		p_masterDistribution = Distributions::getBinomialDistribution(p_xspace.vsize-1, binomial_param);
	} 
	
	
	// if failed to match the space name at all, raise an error
	if (p_xspace.vsize==0) {
		std::cerr << "The space and true distribution are not recognized ";
		exit(1);
	}

}

 
// construction of the omega matrix which specifies the probabilites Pr[client j | x],
// describing the probability that an original data point x is assigned to client j.
// This matrix is used in the non-homogeneous sampling of clients' data points.


/**
 * construct omega1 matrix based on a threshold in the input alphabet.
 * Points below this threshold are assigned to first half of the clients with probability p. 
 * Points above or equal this threshold are mapped to the second half of clients with probability p. 
 **/
DoubleMatrix StatUtils::construct_omega1(const FullGridSpace &p_xspace, const int p_num_clients) {
	int threshold = p_xspace.vsize/2;
	double p = 0.75;
	DoubleMatrix omega(p_xspace.getSize(), p_num_clients);
	for (int x=0; x< threshold; x++) {
		for (int j=0; j<p_num_clients/2; j++) omega.put(x,j, (double) p/(p_num_clients/2) );
		for (int j=p_num_clients/2; j<p_num_clients; j++) omega.put(x,j,(double) (1-p)/(p_num_clients/2));
	}
	for (int x=threshold; x< p_xspace.getSize(); x++) {
		for (int j=0; j<p_num_clients/2; j++) omega.put(x,j, (double) (1-p)/(p_num_clients/2) );
		for (int j=p_num_clients/2; j<p_num_clients; j++) omega.put(x,j, (double) p/(p_num_clients/2) );
	}
	return omega;
}

/**
 * construct omega2 matrix based on a tendency to specific point x1 in the input alphabet:
 * Clients are partitioned in two halfs s1,s2. Original data points that are close to x1 are 
 * more likely to be mapped to half 1 of clients than other points. That is 
 * P(s1 | x) = c exp(- (x - x1)*(x - x1)/b ), P(s2 | x) = 1 - P(s1 | x). 
 * We set b to the smallest value that satisfies P(s1) = P(s2) = 1/2. Note that, 
 * this construction of omega depends on the distribution on X, i.e. masterDistribution.  
 **/
DoubleMatrix StatUtils::construct_omega2(const FullGridSpace &p_xspace, const DoubleRow &p_masterDistribution, const int p_num_clients) {
    if (p_xspace.hsize != 1) {
		std::cerr << "Unable to construct the matrix omega: The given space is not linear";
		exit(1); 
	}
	// special settings chosen by Catuscia: x1 = 5;
    int x1 = 5;
	double c = 2.0; // any initial value > 1
	double b = 0;   // we need to tune this to ensure that  c <= 1
	while (c > 1.0) { 
		b++;
		double sum = 0;
		for (int x = 0; x<p_xspace.getSize(); x++)
			sum = sum + exp(- (x - x1)*(x - x1)/b ) * p_masterDistribution.get(x);
		c = 0.5/sum; // ensures that sum_x P(x) P(s1 |x)  = 1/2, where s1 is the 1st half of servers
	}
	std::cout << "c = " << c << ", b = " << b << "\n";

	DoubleMatrix mat(p_xspace.getSize(), p_num_clients);
	for (int x=0; x<p_xspace.getSize(); x++) {
		for (int j=0; j<p_num_clients/2; j++)
			mat.put(x,j, c * exp(- (x - x1)*(x - x1)/b )/(p_num_clients/2) );
		for (int j=p_num_clients/2; j<p_num_clients; j++)
			mat.put(x,j, (1 - c * exp(- (x - x1)*(x - x1)/b ))/(p_num_clients/2) );
	}
	return mat; 
}


// sampling original points: two types: 
// homogeneous, where data the data is sampled for every client independently 
// non-homogeneous, where all data points are sampled from masterDistribution and mapped to clients according to omega.  

// homogeneous sampling:
// construct a list of empdata objects for all clients, represented by the <client_list>
// Each client in <client_list> is speciefied with dssize (number of data points)
// each empdata object contains original data points sampled from <masterDistribution>
std::vector<GridEmpData> StatUtils::sampleOriginalData_h (const FullGridSpace &p_xspace, const DoubleRow &p_masterDistribution, const std::vector<Client> &p_client_list) {
	std::vector<GridEmpData> ret;
	for (size_t i=0; i<p_client_list.size(); i++)
		for (int j=0; j< p_client_list[i].freq; j++) {
			ret.push_back(sampleData(p_xspace,p_masterDistribution,p_client_list[i].dssize));
		}
	return ret;
}

// /** OLD:  construct a list of numds empdata objects, each has size dssize (this is replaced to use Clients list), 
//  * sampled from the same master distribution (i.e. homogeneousely) 
//  */
// std::vector<GridEmpData> sampleOriginalData_h (FullGridSpace &xspace, DoubleRow &masterDistribution, int &numds, int &dssize) {
// 	std::vector<GridEmpData> ret;
// 	for (int i=0; i < numds; i++) {
// 		std::vector<int> v = Distributions::sampleDiscrete(masterDistribution,dssize);
// 	    std::vector<int> freq_array(xspace.getSize(),0);
// 		for (int i=0; i<dssize;i++) freq_array[ v[i] ]++;
// 		ret.push_back(GridEmpData(freq_array, xspace.hsize, xspace.vsize));
// 	}	
// 	return ret;
// }


/** 
 * construct a list of empdata objects for the clients from given master distribution 
 * and number of points. The original points are sampled from the given master distribution, 
 * and assigned to the clients according to the given omega.
 */ 
// std::vector<GridEmpData> StatUtils::sampleOriginalData_nh (const FullGridSpace &p_xspace, const DoubleRow &p_masterDistribution, const int p_numpoints, const DoubleMatrix &p_omega) {
// }





/** 
 * A special sampling that considers the first and second halfs of clients, each consisting of num_clients/2 clients. 
 * It samples p_numpoints from p_masterDistribution. It directs the sampled points according to p_omega until one of the 
 * halfs reaches exactly p_numpoints/2, then it continues sampling until the other half aquires the same number of points.  
 * The purpose is to apply shuffling to the poinst of each half (if they have same epsilon)
 */
std::vector<GridEmpData> StatUtils::sampleOriginalData_nh_halfs(const FullGridSpace &p_xspace, const DoubleRow &p_masterDistribution, const int p_numpoints, const DoubleMatrix &p_omega) {
	// sampling points and assigning to clients using omega
	std::vector<std::vector<int>> freq_array;
	//initialization
	int numClients = p_omega.getColumns();   // number of clients is the same asthe number of columns in p_omega.
	for (int j=0; j<numClients; j++) {
		freq_array.push_back(std::vector<int>(p_xspace.getSize(),0)); // each client has a vector that contains the frequencies of sampled values
	}

	std::vector<int> v = Distributions::sampleDiscrete(p_masterDistribution, p_numpoints); // sample all points
	int half0_count=0; int half1_count=0;  
	int i = 0;  // index of the sampled data point.

	// sample until one of the two halfs reachs p_numpoints/2
	while ((half0_count < p_numpoints/2) && (half1_count < p_numpoints/2) ) {
		DoubleRow omega_row = p_omega.selectRow(v[i]); // get the distribution, over clients, corresponding to the value of point i.
		std::vector<int> s = Distributions::sampleDiscrete(omega_row, 1); // sample a client s[0] from the distribution on them
		freq_array[ s[0] ][ v[i] ]++;
		if (s[0] < numClients/2) {
			half0_count++;
		}
		else {
			half1_count++;
		}
		i++;
	}
	// check the number of points at each half
	std::cout << "half0 = " << half0_count << "\n";
	std::cout << "half1 = " << half1_count << "\n";
	
	// continue sampling for the smaller half to ensure exact balance.
	if (half0_count < p_numpoints/2) {
		std::cout << "completing the data of half0 ...\n";
		while (half0_count < p_numpoints/2) {
			std::vector<int> v = Distributions::sampleDiscrete(p_masterDistribution,1); // sample one data point v[0]
			DoubleRow omega_row = p_omega.selectRow(v[0]); // distribution over clients depending on the value of the data point
			std::vector<int> s = Distributions::sampleDiscrete(omega_row, 1);// select the server
			if (s[0]< numClients/2) {
				half0_count++;
				freq_array[ s[0] ][ v[0] ]++;
			}
		}
	}
	if (half1_count < p_numpoints/2) {
		std::cout << "completing the data of half1 ...\n";
		while (half1_count < p_numpoints/2) {
			std::vector<int> v = Distributions::sampleDiscrete(p_masterDistribution,1); // sample one data point v[0]
			DoubleRow omega_row = p_omega.selectRow(v[0]); // distribution over clients depending on the value of the data point
			std::vector<int> s = Distributions::sampleDiscrete(omega_row, 1);// v is a vector of singe entry v[0]
			if (s[0] >= numClients/2) {
				half1_count++;
				freq_array[ s[0] ][ v[0] ]++;
			}
		}
	}

	std::cout << "finished: half0_count = " << half0_count << "\n"; 
	std::cout << "half1_count = " << half1_count << "\n"; 
	// omega.write_to_file("omega.txt");
	
	std::vector<GridEmpData> ret;
	for (int j=0; j<numClients; j++) 
		ret.push_back(GridEmpData(freq_array[j], p_xspace.hsize, p_xspace.vsize));
		
	return ret;
}


/** sample p_numpoints data points from the master distribution on p_xspace and return GridEmpData */
GridEmpData StatUtils::sampleData (const FullGridSpace &p_xspace, const DoubleRow &p_masterDistribution, const int p_numpoints) {
		std::vector<int> v = Distributions::sampleDiscrete(p_masterDistribution,p_numpoints);
	    std::vector<int> freq_array(p_xspace.getSize(),0); // initialize to 0
		for (int i=0; i<p_numpoints;i++) freq_array[ v[i] ]++;
		return GridEmpData(freq_array, p_xspace.hsize,p_xspace.vsize);
}

/** sample p_numpoints original data points from a uniform distribution on the range p_min to p_max and return GridEmpData object*/
GridEmpData StatUtils::sampleUniformData (const FullGridSpace &p_xspace, const int p_numpoints, const int p_min, const int p_max) {
	    std::vector<int> freq_array(p_xspace.vsize,0); // initialized to 0
	    	for (int i=0;i<p_numpoints;i++) {
	    		int sampled = Distributions::sampleUniformInt(p_min, p_max);
	    		//geometric: use 0.2..0.4. Rappor: 0.3 .. 0.6
	    		freq_array[sampled]++;
			}
		return GridEmpData(freq_array, p_xspace.hsize,p_xspace.vsize);
}

