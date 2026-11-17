#include "iostream"
#include "PrivMetric.h"

PrivMetric::PrivMetric(double p_eps)
: eps(p_eps) {
}

/**
 * returns a string representation for the given statistical distance type
 */
std::string PrivMetric::privTypeToString(PrivacyType p_priv_type) {
	switch(p_priv_type) {
	case PRIV_GEO : return "PRIV_GEO";
	case PRIV_MAN : return "PRIV_MAN";
	case PRIV_MAX : return "PRIV_MAX";
	default :
		std::cerr << "error: invalid privacy type";
		exit(0);
	}
	return "";
}

