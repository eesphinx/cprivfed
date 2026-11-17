/*
 * Distributions.h
 *
 *  Created on: Jun 27, 2019
 *      
 */

#ifndef STAT_DISTRIBUTIONS_DISTRIBUTIONS_H_
#define STAT_DISTRIBUTIONS_DISTRIBUTIONS_H_

#include <vector>


enum DistributionType { UNIFORM, BINOMIAL };


class DoubleRow;

class Distributions {
public:
	/**
	 * returns a vector of q integers sampled from a probability distribution p_probabilities
	 * @param p_probabilities
	 * @param q
	 * @return
	 */
	static std::vector<int> sampleDiscrete(const DoubleRow& p_probabilities, int q);


	/**
	 * samples an integer between 0 and n using binomial distributions with parameters n,p
	 */
	static int sampleBinomial(int p_maxnum, double p);

	// /** get the probability of element k in [0 .. p_max]  according to the binomial distribution
	//  * on the range [0 .. p_max] with parameter p.
	//  */
	// static double getBinomialProbability(int p_k, int p_max, double p);

	static DoubleRow getBinomialDistribution(int p_max, double p);

	// Discretized Beta distribution over {0, 1, ..., p_max}
    static DoubleRow getDiscBetaDistribution(int p_max, double alpha_param, double beta_param);

	// Discretized Gauss distribution over {0, 1, ..., p_max}
    static DoubleRow getDiscGaussDistribution(int p_max, double mu_param, double sigma_param);
	/**
 	* return a vector of q integers between 1 and p_n sampled from a zipf distributions with parameters p_n>=1 and p_s>=0
 	*/
	static std::vector<int> sampleZipf(double p_alpha, int p_n, int q);

    /**
	 * samples an integer between p_min and p_max using the uniform distribution.
	 */
	static int sampleUniformInt(int p_min, int p_max);

	/** samples a double from uniform distribution on [p_min and p_max] */
	static double sampleUniformDouble(double p_min, double p_max);

	/** samples a double from normal (Gaussian) distribution with (p_mu, p_sigma) */
	static double sampleNormalDouble(double p_mu, double p_sigma);

	/** sample a double from laplace distribution with (p_mu, p_beta)*/
	static double sampleLaplaceDouble(double p_mu, double p_beta);

	/** sample an integer in {...,-2,-1,0,1,2,...} from a double sided 
	 * (symmetric) geometric distribution, where P[k] = (1-a)/(1+a) a^|k|
	 * which is also the geometric mechanism. */
	static int sampleSymGeometricInt(double p);
	
	/** returns a string representation for the given statistical distance type */
	static std::string distributionTypeToString(DistributionType p_type);

	/** normalize a vector to be a valid distribution */
	static void normalize(DoubleRow& p_vec);

	/** project the given distribution onto the simplex using the method in:
	 *  "Projection onto the probability simplex: An efficient algorithm with a simple proof, and an application" by Weiran Wang et al. **/
	static void project_onto_simplex(DoubleRow& p_vec);


};



#endif /* STAT_DISTRIBUTIONS_DISTRIBUTIONS_H_ */
