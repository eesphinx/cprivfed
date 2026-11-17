/*
 *  Created on: Aug 28, 2025
 *      Author: ehab
 * 
 * This program receives local eps0, delta, and returns the amplified epsilon
 * using the expression for generic LDP mechanism in the clones paper 
 * Clones Feldman2020 (arxiv version dated 7 sep 2021) - Theorem 3.1
 *
 */
//include <cmath>
#include <iostream>

// use as krr_clones_amp <eps0> <delta> <n> 
int main(int argc, char* argv[]) {
	if (argc != 4) {
        std::cout << "number of arguments = " << argc;
        std::cerr << "Usage: " << argv[0] 
        << " <eps0> <n> <delta> \n"; 
        return 1;
    }
	// To debug:
	// for (int i = 0; i< argc ; i++)
	// 	std::cout << "arg[" << i << "] = " << argv[i] << "\n";
	// exit(0);

	double eps0 = std::atof(argv[1]);
	double delta = std::atof(argv[2]);
    int n = std::atoi(argv[3]);

  	std::cerr << "generic clones2020: eps0 threshold is " 
	          << log( n/(16 * log(2/delta) )) << "\n";
	if ( eps0 > log( n/(16 * log(2/delta)) ) ) {
	   std::cerr << " eps is not amplified \n ";
       std::cout << eps0 << std::endl;
	}
  	else {
    	double term1 = (exp(eps0) - 1) / (exp(eps0) + 1);
    	double term2 = 8 * sqrt( exp(eps0) * log(4/delta) ) / sqrt(n);
    	double term3 = 8 * exp(eps0) / n;
		double eps = log( 1 + term1 * (term2 + term3)  );
		std::cerr << " eps is amplified: " << eps0 << " -> " << eps << " \n"; 
    	std::cout << eps << std::endl;
	}
}

