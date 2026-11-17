/*
 *  Created on: Aug 28, 2025
 *      Author: ehab
 * 
 * This program receives local eps0, delta, and returns the amplified epsilon
 * using the expression for generic LDP mechanism in our GIBU paper Theorem 3
 *
 */
//include <cmath>
#include <iostream>

// use as krr_clones_amp <eps0> <delta> <n> 
int main(int argc, char* argv[]) {
	if (argc != 4) {
        std::cerr << "number of arguments = " << argc;
        std::cerr << "Usage: " << argv[0] 
        << " <eps0> <n> <delta> \n"; 
        return 1;
    }
	// To debug:
	// for (int i = 0; i< argc ; i++)
	// 	std::cout << "arg[" << i << "] = " << argv[i] << "\n";
	// exit(0);

	double eps0 = std::atof(argv[1]);
	double deltap = std::atof(argv[2]);
    int n = std::atoi(argv[3]);
	
	double p = 1.0/(exp(eps0)+1); 
    double q = exp(eps0)/(exp(eps0)+1);
    double delta = deltap/(2*q-1);        // adjusting delta to account for its amplification 

    // Our method (uses p,q,n,delta)
	std::cerr << "our generic clones2022: eps0 threshold is " <<
	 			log( (n-1)/(6.0*log(2.0/delta) - 1) ) << "\n";
	if (2*p < (12.0/(n-1))*log(2.0/delta)) {
	   std::cerr << " eps is not amplified \n ";
       std::cout << eps0 << std::endl;
	}
    else {
	    double c_min = 2*p*(n-1) - sqrt(6*p*(n-1)*log(2.0/delta));
	    double c_max = 2*p*(n-1) + sqrt(6*p*(n-1)*log(2.0/delta));
	    double eps1 = log(c_max/2 + sqrt(c_max/2 * log(2.0/delta) +1) ) 
                        - log(c_min/2 - sqrt(c_max/2 * log(2.0/delta)));
        double eps = log(1 + (2 * q - 1) * ( exp(eps1) -1 ) ); 
		std::cerr << " eps is amplified: " << eps0 << " -> " << eps << " \n"; 
    	std::cout << eps << std::endl;
    }
}

