/*
 * StatUtils.h
 *
 *  Created on: Sep 14, 2026
 *      
 */

#ifndef SRCLIB_STAT_TOOLS_STATUTILS_H_
#define SRCLIB_STAT_TOOLS_STATUTILS_H_

/*
Utility functions used in the experiments.
*/


// client info (mechanism and dataset size) that may be repeated freq (number of clients) times
struct Client {
    int mech_idx;  // the index of the mechanism in the separate mechansim list mechanism_ptr_list
	int dssize;
	int freq;
};



class StatUtils {
	public:
	/** construction of the xspace and master_distribution using the given space_str */
	static void construct_space(const std::string space_str, FullGridSpace &p_xspace, DoubleRow &p_masterDistribution);
	
	/**
	 * construct omega1 matrix based on a threshold in the input alphabet.
	 * Points below this threshold are assigned to first half of the clients with probability p. 
	 * Points above or equal this threshold are mapped to the second half of clients with probability p. 
	 */
	static DoubleMatrix construct_omega1(const FullGridSpace &p_xspace, const int p_num_clients);

	/**
	 * construct omega2 matrix based on a tendency to specific point x1 in the input alphabet:
	 * Clients are partitioned in two halfs s1,s2. Original data points that are close to x1 are 
	 * more likely to be mapped to half 1 of clients than other points. That is 
	 * P(s1 | x) = c exp(- (x - x1)*(x - x1)/b ), P(s2 | x) = 1 - P(s1 | x). 
	 * We set b to the smallest value that satisfies P(s1) = P(s2) = 1/2. Note that, 
	 * this construction of omega depends on the distribution on X, i.e. masterDistribution. 
	 */
 	static DoubleMatrix construct_omega2(const FullGridSpace &p_xspace, const DoubleRow &p_masterDistribution, const int p_num_clients);

	/** homogeneous sampling:
	 * construct a list of empdata objects for all clients, represented by the <client_list> 
	 * Each client in <client_list> is speciefied with dssize (number of data points) 
	 * each empdata object contains original data points sampled from <masterDistribution>
	 */
	static std::vector<GridEmpData> sampleOriginalData_h(const FullGridSpace &p_xspace, const DoubleRow &p_masterDistribution, const std::vector<Client> &p_client_list);


	/** 
	 * non-homogeneous sampling:
	 * construct a list of empdata objects for the clients from given master distribution 
	 * and number of points. The original points are sampled from the given master distribution, 
	 * and assigned to the clients according to the given omega. Number od clients are exactly the 
	 * number of columns in the matrix <p_omega>.
	 */ 
	// static std::vector<GridEmpData> sampleOriginalData_nh(const FullGridSpace &p_xspace, const DoubleRow &p_masterDistribution,  const int p_numpoints, const DoubleMatrix &p_omega);

	/** 
 	 * A special sampling that considers the first and second halfs of clients, each consisting of num_clients/2 clients. 
 	 * It samples p_numpoints from p_masterDistribution. It directs the sampled points according to p_omega until one of the 
   	 * halfs reaches exactly p_numpoints/2, then it continues sampling until the other half aquires the same number of points.  
 	 * The purpose is to apply shuffling to the poinst of each half (if they have same epsilon)
 	 */
	static std::vector<GridEmpData> sampleOriginalData_nh_halfs(const FullGridSpace &xspace, const DoubleRow &p_masterDistribution, const int p_numpoints, const DoubleMatrix &p_omega);

	/** sample p_numpoints data points from the master distribution on p_xspace and return GridEmpData */
	static GridEmpData sampleData (const FullGridSpace &p_xspace, const DoubleRow &p_masterDistribution, const int p_numpoints);

	/** sample p_numpoints original data points from a uniform distribution on the range p_min to p_max and return GridEmpData object*/
	static GridEmpData sampleUniformData (const FullGridSpace &p_xspace, const int p_numpoints, const int p_min, const int p_max);
};

#endif /* SRCLIB_STAT_TOOLS_STATUTILS_H_ */


