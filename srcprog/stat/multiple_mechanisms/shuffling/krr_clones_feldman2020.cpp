/*
 *  Created on: Aug 28, 2025
 *      Author: ehab
 * 
 * This program receives local eps0, delta, and returns the amplified epsilon
 * using the expression for the kRR mechanism in the clones paper 
 * Feldman et al 2020 (arxiv version dated 7 sep 2021) - Crollary 4.2
 *
 */
//include <cmath>
#include <iostream>

// use as krr_clones_feldman2020 <eps0> <delta> <n> <k> 
int main(int argc, char* argv[]) {
	if (argc != 5) {
        std::cout << "number of arguments = " << argc;
        std::cerr << "Usage: " << argv[0] 
        << " <eps0> <delta> <n> <k> \n"; 
        return 1;
    }
	// To debug:
	// for (int i = 0; i< argc ; i++)
	// 	std::cout << "arg[" << i << "] = " << argv[i] << "\n";
	// exit(0);

	double eps0 = std::atof(argv[1]);
	double delta = std::atof(argv[2]);
    int n = std::atoi(argv[3]);
    int k = std::atoi(argv[4]);
  	
    std::cerr << "krr specific clones2020: eps0 threshold is " 
	          << log( n/(16 * log(2/delta)) ) << "\n";
    if ( eps0 > log( n/(16 * log(2/delta)) ) ) {
	   std::cerr << " eps is not amplified \n ";
       std::cout << eps0 << std::endl;
    }
    else {
        double term1 = (4 * sqrt(2*(k+1) * log(4/delta))) / sqrt((exp(eps0) + k + 1) * k * n);
        double term2 = 4 * (k+1) / (k * n);
        double eps = log(1 + (exp(eps0) - 1) * (term1 + term2));
		std::cerr << " eps is amplified: " << eps0 << " -> " << eps << " \n"; 
    	std::cout << eps << std::endl;
    }
}

