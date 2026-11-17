#include <iostream>
#include <cmath>
#include <random>
#include <vector>
#include "Distributions.h"
#include "DoubleRow.h"


/** get a random number generator gen which is initialized by a random seed obtained using the random device rd */
static std::random_device rd;  // random device engine, usually based on /dev/random on UNIX-like systems
static std::mt19937 gen(rd()); // initialize Mersennes' twister using rd to generate the seed

/**
 * returns a vector of q integers sampled from {0,1,...,n} according to the probability distribution p_probabilities
 * n = len(p_probabilities)-1
 */
std::vector<int> Distributions::sampleDiscrete(const DoubleRow& p_probabilities, int q) {
	std::vector<int> samples;
	std::vector<double> probabilities = p_probabilities.toVector();  // convert the DoubleRow to Vector<double>
	std::discrete_distribution<int> distribution (probabilities.begin(), probabilities.end());
	for (int i = 0; i < q; i++) {
		int s = distribution(gen);
		samples.push_back(s);
	}
	return samples;
}


/**
 * samples an integer between 0 and p_max using binomial distributions with parameters p_max and p
 */
int Distributions::sampleBinomial(int p_max, double p) {
	//std::random_device rd;  // random device engine, usually based on /dev/random on UNIX-like systems
	//std::mt19937 gen(rd()); // initialize Mersennes' twister using rd to generate the seed

	std::binomial_distribution<int> distribution(p_max,p);
	return distribution(gen);
}

// Compute n choose k using the formula n/1 * (n-1)/2 * ... * (n-k+1)/k
double binom_coeff(int n, int k) {
    double res = 1.0;
    for (int i = 1; i <= k; i++) {
        res *= (n - i + 1) / (double)i;
    }
    return res;
}
// This constructs the binomial distribution on 0..p_max
DoubleRow Distributions::getBinomialDistribution(int p_max, double p) {
	DoubleRow dist(p_max + 1);
    for (int j = 0; j <= p_max; j++) {
        double coeff = binom_coeff(p_max, j);
        dist.put(j, coeff * std::pow(p, j) * std::pow(1 - p, p_max - j) );
    }
    return dist;
}

// Discretized Gauss distribution over {0, 1, ..., p_max}
DoubleRow Distributions::getDiscGaussDistribution(int p_max, double mu_param, double sigma_param) {
    const double PI = 3.14159265358979323846;
    double coeff = 1.0 / (sigma_param * std::sqrt(2.0 * PI));

	DoubleRow dist(p_max + 1);
    // Compute unnormalized probabilities
    for (int x = 0; x <= p_max; ++x)
    {
        // Map i to [0,1] and use the Gaussian PDF
    	double exponent = -((x - mu_param) * (x - mu_param)) / (2.0 * sigma_param * sigma_param);
        dist.put(x, coeff * std::exp(exponent));
    }

	// Normalize and return
	return dist.div(dist.sum());
}



// Discretized Beta distribution over {0, 1, ..., p_max}
DoubleRow Distributions::getDiscBetaDistribution(int p_max, double alpha_param, double beta_param)
{
    DoubleRow dist(p_max + 1);
    double B = std::tgamma(alpha_param) * std::tgamma(beta_param) / std::tgamma(alpha_param + beta_param);

    // Compute unnormalized probabilities
    for (int j = 0; j <= p_max; ++j)
    {
        double x = (double)(j) / p_max;
        // Handle boundaries explicitly
        if ((x == 0.0 && alpha_param < 1.0) ||
            (x == 1.0 && beta_param < 1.0))
        {
            dist.put(j,0.0);
        }
        else
        {
            dist.put(j, std::pow(x, alpha_param - 1.0)
                 * std::pow(1.0 - x, beta_param - 1.0)
                 / B);
        }
    }
	if (dist.sum() <= 0 ) {
		std::cerr << "Unable to construct Beta distribution";
		exit(1);
	}

    // Normalize and return
	return dist.div(dist.sum());
}


/**
 * return a vector of q integers between 1 and p_n sampled from a zipf distributions with parameters p_n>=1 and p_s>=0
 */
std::vector<int> Distributions::sampleZipf(double p_alpha, int p_n, int q)
{
	double c = 0;          // The normalization constant
  	// Compute normalization constant on first call only
  	for (int i=1; i<=p_n; i++)
		c = c + (1.0 / std::pow((double) i, p_alpha));
	c = 1.0 / c;
	
	// compute the cdf
	std::vector<double> cdf(p_n + 1, 0.0);  // the value of the cummulative density function at 1,2,...,p_n. Note cdf[0] = 0.  
	cdf[1] = c; 
	for (int i=2; i<= p_n; i++)
		cdf[i] = cdf[i-1] + c / std::pow((double) i, p_alpha);

	// generate a uniform random number (0 < z < 1)
  	std::uniform_real_distribution<> dis(0.0, 1.0);  // define a uniform distribution between 0 and 1

	// sample q numbers
	std::vector<int> zipf_values;
	zipf_values.reserve(q);
	for (int i = 0; i < q; i++) {
		double z = dis(gen);
		// Map z to the value using binary search
		auto it = std::lower_bound(cdf.begin() + 1, cdf.end(), z); // find the smallest k such that z <= cdf(k) 
		int k = std::distance(cdf.begin(), it);
		zipf_values.push_back(k);
		
		// manual search 
  		//int low = 1, high = p_n, mid;
  		//while (low <= high) {
    	//	mid = (low + high) / 2;
    	//	if (cdf[mid] >= z && cdf[mid - 1] < z) {
        //		zipf_values.push_back(mid);
        //		break;
    	//	} else if (cdf[mid] >= z) {
        //		high = mid - 1;
    	//	} else {
        //		low = mid + 1;
    	//	}
		//}
	}
  return(zipf_values);
}



/**
 * samples an integer between p_min and p_max using the uniform distribution
 */
int Distributions::sampleUniformInt(int p_min, int p_max) {
	// construct a trivial random generator engine from a time-based seed. Note the seed must change in every call
	//unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	//std::default_random_engine gen (seed);

	std::uniform_int_distribution<int> distribution(p_min,p_max);
	return distribution(gen);
}

/** samples a double from uniform distribution on [p_min and p_max] */
double Distributions::sampleUniformDouble(double p_min, double p_max) {
	std::uniform_real_distribution<double> distribution(p_min, p_max);
	return distribution(gen);
}

/** samples a double from normal (Gaussian) distribution with (mu, sigma) */
double Distributions::sampleNormalDouble(double p_mu, double p_sigma) {
	std::normal_distribution<double> distribution(p_mu,p_sigma);
	return distribution(gen);
}

/** samples a double from Laplace distribution with (p_mu, p_beta) */
double Distributions::sampleLaplaceDouble(double p_mu, double p_beta) {
	std::exponential_distribution<double> exp_distribution(1/p_beta);
	double r1 = exp_distribution(gen);
	double r2 = exp_distribution(gen);
	return p_mu + r1 - r2;
}

/** sample an integer in {...,-2,-1,0,1,2,...} from a double sided 
 * (symmetric) geometric distribution, where P[k] = (1-a)/(1+a) a^|k|
 * which is also the geometric mechanism. */
int Distributions::sampleSymGeometricInt(double a) {
	std::uniform_real_distribution<double> unif(0.0, 1.0);
    // The following implements the distribution is P[k]=(1-a) a^k.
	std::geometric_distribution<int> geom(1.0-a); // support is 0,1,2,...
	if (unif(gen) < (1.0 - a)/(1.0 + a))
        return 0;
	// Otherwise...
	int k = geom(gen)+1;
	// observe that P[0]= (1.0 - a)/(1.0 + a), and for k=1,2,...
	// P[k]=( 1- (1.0 - a)/(1.0 + a) ) * (1-a) a^(k-1) = 2 (1-a)/(1+a) a^k
	// then we choose k and -k with equal probability. 
	if (unif(gen) < 0.5)
        return k;
	else
		return -k;
}

/** returns a string representation for the given statistical distance type */
std::string Distributions::distributionTypeToString(DistributionType p_type) {
	switch(p_type) {
	case UNIFORM : return "UNIFORM";
	case BINOMIAL : return "BINOMIAL";
	default :
		std::cerr << "invalid distribution type \n";
		exit(0);
	}
	return "";
}

/** normalize a vector to be a valid distribution */
void Distributions::normalize(DoubleRow& p_vec) {
	p_vec.muli(p_vec.ge(0));    // set -ve entries to zeros
	double sum = p_vec.sum();
	p_vec.divi(sum);
	//System.out.println("euclidean distance to the -ve distribution = "+ pr_x.distance2(pr_x1));
}


/** project the given distribution onto the simplex using the method in:
 *  "Projection onto the probability simplex: An efficient algorithm with a simple proof, and an application" by Weiran Wang et al. **/
void Distributions::project_onto_simplex(DoubleRow& p_vec) {
	// projection
	std::vector<double> y = p_vec.toVector();
	std::sort(y.begin(), y.begin() + y.size());  // sorts the vector in ascending order
	double lambda =0;
	double sum_u =0;
	for (size_t j=1; j<=y.size(); j++) {
		double u = y[y.size()-j];
		sum_u = sum_u + u;
		if (u + (1- sum_u)/j >0)
			lambda = (1- sum_u)/j; // this turns to be lambda in the end
	}
	p_vec.addi(lambda);         // add lambda to every element
	p_vec.muli(p_vec.ge(0));    // truncate -ve elements to zeros

	/* debug:
	System.out.println("After projection:\n"+pr_x);
	System.out.println("Total probability = "+pr_x. sum());
	System.out.println("Size =" + pr_x.length);
	for (int i=0;i<pr_x.length;i++)
		if (pr_x.get(i)<0)
			System.out.println("-ve element at "+i +":"+pr_x.get(i));
	 */
	//System.out.println("euclidean distance to the -ve distribution = "+ pr_x.distance2(pr_x1));
}

