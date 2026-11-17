#include <iostream>
#include "DoubleRow.h"
#include "Distributions.h"
#include "Estimator_EMP.h"


	/**
	 * Construct an estimator for the distribution on a given space of secrets,
	 * and a set of histograms resulting from obfuscation using central mechanisms. 
	 */
	Estimator_EMP::Estimator_EMP (const std::vector<DoubleRow>& p_histograms) {
        histograms = p_histograms;
    }

    /**
	 * Returns an estimate by aggregating the underlying histograms, 
     * and post-processing using 
     * 'N': truncating negative components to 0 and then normalize, 
     * 'P': normalize then project onto the simplex, Or 
     * 'O': for only normalization without correcting -ve components.
	 */
	DoubleRow Estimator_EMP::getEstimate (char p_post_processing) const{
        DoubleRow overallHistogram(histograms[0].getLength()); // initialized to 0
        for (size_t i=0; i<histograms.size(); i++) {
            overallHistogram.addi(histograms[i]);
        }
        if (p_post_processing == 'O') 
            return overallHistogram.div(overallHistogram.sum());
	    else if (p_post_processing == 'N') {
		    Distributions::normalize(overallHistogram);
            return overallHistogram;
        }
	    else if (p_post_processing == 'P') {
            overallHistogram.divi(overallHistogram.sum()); 
		    Distributions::project_onto_simplex(overallHistogram);
            return overallHistogram;
        }
	    else {
		    std::cerr << "invalid post-processing parameter ("<< p_post_processing << ") for the EMP estimator";
		    exit(1);
	    }
    }

    /** apply given weights before post-processing */
	DoubleRow Estimator_EMP::getEstimate (const std::vector<double> p_list_variances, char p_post_processing) const {
        DoubleRow overallHistogram(histograms[0].getLength());
        for (size_t i=0; i<histograms.size(); i++) {
            overallHistogram.addi(histograms[i].div(p_list_variances[i]));
        }

        if (p_post_processing == 'O') 
            return overallHistogram.div(overallHistogram.sum());
	    else if (p_post_processing == 'N') {
		    Distributions::normalize(overallHistogram);
            return overallHistogram;
        }
	    else if (p_post_processing == 'P') {
            overallHistogram.divi(overallHistogram.sum()); 
		    Distributions::project_onto_simplex(overallHistogram);
            return overallHistogram;
        }
	    else {
		    std::cerr << "invalid post-processing parameter ("<< p_post_processing << ") for the EMP estimator";
		    exit(1);
	    }
    }

	/** returns the pseudo estimate defined by the aggregae divided by p_count */
	DoubleRow Estimator_EMP::getPseudoEstimate(int p_count) const {
        DoubleRow overallHistogram(histograms[0].getLength()); // initialized to 0
        for (size_t i=0; i<histograms.size(); i++) {
            overallHistogram.addi(histograms[i]);
        }
        return overallHistogram.div(p_count);
    }

	/** returns the pseudo estimate defined by the weighted aggregae divided by the weighted sum of subcounts */
	DoubleRow Estimator_EMP::getPseudoEstimate(const std::vector<double> p_list_variances, std::vector<int> p_subcounts) const {
        if (histograms.size() != p_list_variances.size() || histograms.size() != p_subcounts.size()) {
            std::cerr << "error: the vectors histograms, weights, and subcounts are of different sizes";
            exit(1);
        }
        DoubleRow overallHistogram(histograms[0].getLength());
        double overallCount = 0;
        for (size_t i=0; i<histograms.size(); i++) {
            overallHistogram.addi(histograms[i].div(p_list_variances[i]));
            overallCount = overallCount + (double)p_subcounts[i]/p_list_variances[i];
        }
        return overallHistogram.div(overallCount);
    }

	/** analytic Mean SQD without weighting */
    double Estimator_EMP::getPseudoEstimateAnalyticMeanSQD(const DoubleRow p_original_distribution, const std::vector<double> p_list_variances, int p_numPoints) {
        double sumVariances = 0;
        for (int j=0; j< p_list_variances.size(); j++) sumVariances = sumVariances +  p_list_variances[j];
        double ret = p_original_distribution.norm2();
        ret = (1 - ret*ret)/ (double) p_numPoints;
        ret = ret + sumVariances * p_original_distribution.getLength() / ((double) p_numPoints * (double) p_numPoints);
        return ret;
    }


	/** analytic Mean SQD with weighting */
    double Estimator_EMP::getPseudoEstimateAnalyticMeanSQD(const DoubleRow p_original_distribution, const std::vector<double> p_list_variances, std::vector<int> p_subcounts) {
        if (p_list_variances.size() != p_subcounts.size()) {
            std::cerr << "error: the list of variances , and list of subcounts are of different sizes";
            exit(1);
        }
        double sumInvVariances = 0; // sum_j 1/varj
        double wn = 0;              // sum_j nj/varj
        double r = 0;            // sum_j nj/varj^2
        for (int j=0; j< p_list_variances.size(); j++) {
            // std::cout << "p_list_variances["<<j<<"] = " << p_list_variances[j] << "\n";
            // std::cout << "p_subcounts["<<j<<"] = " << p_subcounts[j] << "\n";
            sumInvVariances = sumInvVariances +  1.0/p_list_variances[j];
            wn = wn + (double) p_subcounts[j] / p_list_variances[j];
            r = r + p_subcounts[j]/(p_list_variances[j] * p_list_variances[j]); 
            // std::cout << "sumInvVariances = " << sumInvVariances << "\n";
            // std::cout << "r = " << r << "\n";
            // std::cout << "wn = " << wn << "\n";
        }
        wn = wn / sumInvVariances;
        // std::cout << "sumInvVariances = " << sumInvVariances << ", r = " << r << ", wn = " << wn << "\n";
        double ret = p_original_distribution.norm2();
        // std::cout << "theta^2 = " << ret << "\n";  
        ret = (1 - ret*ret) * r / (wn * wn * sumInvVariances * sumInvVariances);
        // std::cout << "sampling error = " << ret << "\n";
        // std::cout << "k = " << p_original_distribution.getLength() << "\n";
        ret = ret + (1.0/sumInvVariances) * (double) p_original_distribution.getLength() / (wn * wn);
        // std::cout << "ret = " << ret; 
        return ret;
    }
