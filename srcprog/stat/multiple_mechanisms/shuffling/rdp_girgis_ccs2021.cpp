#include <cmath>
#include <iostream>
#include <limits>

/*
 * Shuffle amplification using Girgis Theorem 2 (Upper Bound 2)
 * 
 * ε(λ) ≤ 1/(λ-1) log( e^{λ²(e^{ε₀} - 1)²/𝑛̄} + e^{ε₀ λ - (n-1)/(8e^{ε₀})} )
 * where 𝑛̄ = floor((n-1)/(2e^{ε₀})) + 1
 */
double girgis_shuffle_rdp(double eps0, double lambda, int n) {
    if (lambda <= 1.0) return std::numeric_limits<double>::infinity();
    
    double exp_eps = std::exp(eps0);
    
    // Compute 𝑛̄ = floor((n-1)/(2e^{ε₀})) + 1
    double n_bar = std::floor((n - 1.0) / (2.0 * exp_eps)) + 1.0;
    
    // First term: λ²(e^{ε₀} - 1)²/𝑛̄
    double term1 = (lambda * lambda) * (exp_eps - 1.0) * (exp_eps - 1.0) / n_bar;
    
    // Second term: ε₀ λ - (n-1)/(8e^{ε₀})
    double term2 = eps0 * lambda - (n - 1.0) / (8.0 * exp_eps);
    
    // Combine: log( e^{term1} + e^{term2} )
    // Use log-sum-exp trick for numerical stability
    double max_exp = std::max(term1, term2);
    double log_sum = std::log(std::exp(term1 - max_exp) + std::exp(term2 - max_exp)) + max_exp;
    
    // Final RDP: 1/(λ-1) * log_sum
    return log_sum / (lambda - 1.0);
}

/*
 * Convert RDP to (ε,δ)-DP using Girgis Lemma 1
 * 
 * ε = ε(λ) + (log(1/δ) + (λ-1)log(1 - 1/λ) - log(λ)) / (λ-1)
 */
double girgis_lemma1_conversion(double eps_rdp, double lambda, double delta) {
    if (lambda <= 1.0) return std::numeric_limits<double>::infinity();
    
    // Lemma 1 conversion
    double correction = (lambda - 1.0) * std::log(1.0 - 1.0 / lambda) - std::log(lambda);
    double lemma1_term = (std::log(1.0 / delta) + correction) / (lambda - 1.0);
    
    return eps_rdp + lemma1_term;
}

/*
 * Optimize over lambda to find minimal ε
 */
double find_minimal_epsilon(double eps0, int n, double delta) {
    double best_eps = std::numeric_limits<double>::infinity();
    
    // Search λ from 2 to 1000 with logarithmic spacing
    const int steps = 1000;
    const double lambda_min = 2.0;
    const double lambda_max = 1000.0;
    
    double log_min = std::log(lambda_min);
    double log_max = std::log(lambda_max);
    
    for (int i = 0; i <= steps; i++) {
        double lambda = std::exp(log_min + (log_max - log_min) * i / steps);
        
        // Compute RDP using Girgis Theorem 2
        double eps_rdp = girgis_shuffle_rdp(eps0, lambda, n);
        
        // Convert to (ε,δ)-DP using Lemma 1
        double eps = girgis_lemma1_conversion(eps_rdp, lambda, delta);
        
        // Update best if better
        if (eps < best_eps && !std::isinf(eps) && !std::isnan(eps)) {
            best_eps = eps;
        }
    }
    
    return best_eps;
}

/*
 * Compute amplified epsilon for given parameters
 */
double compute_amplified_epsilon(double eps0, double delta, int n) {
    // Edge cases
    if (eps0 <= 0.0) return 0.0;
    if (n <= 1) return eps0;  // No amplification with single user
    
    // // For very large ε₀, amplification is negligible
    // if (eps0 > 10.0) {
    //     return eps0 + std::sqrt(2.0 * std::log(1.0 / delta) / n);
    // }
    
    // Use Girgis Theorem 2 + Lemma 1
    return find_minimal_epsilon(eps0, n, delta);
}

int main(int argc, char* argv[]) {
	if (argc != 4) {
        std::cout << "number of arguments = " << argc;
        std::cerr << "Usage: " << argv[0] 
        << " <eps0> <delta> <n> \n"; 
        return 1;
    }
	// To debug:
	// for (int i = 0; i< argc ; i++)
	// 	std::cout << "arg[" << i << "] = " << argv[i] << "\n";
	// exit(0);

	double eps0 = std::atof(argv[1]);
	double delta = std::atof(argv[2]);
    int n = std::atoi(argv[3]);

    double eps = compute_amplified_epsilon(eps0, delta, n);
    std::cout << eps << std::endl;
}

