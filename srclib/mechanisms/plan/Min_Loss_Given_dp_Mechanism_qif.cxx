/*
 * Min_Loss_Given_dp_Mechanism.cpp
 *
 *  Created on: Feb 1, 2020
 *      Author: ehab
 */
#include <qif>
#include "Prior.h"
#include "PrivMetric.h"
#include "Loss.h"
#include "FullGridSpace.h"
#include "Min_Loss_Given_dp_Mechanism.h"

Min_Loss_Given_dp_Mechanism::Min_Loss_Given_dp_Mechanism(const Prior& p_prior, const FullGridSpace& p_xspace,
		const PrivMetric* p_privMetric, const Loss* p_lossMetric, const std::string p_mechFile)
:PMatrixMechanism(p_xspace.hsize, p_xspace.vsize) {

	// priv metrics in qif format
	qif::Metric<double,uint> priv = [=](uint a, uint b) {
		int a_x = a % hsize; int a_y = a/hsize;
		int b_x = b % hsize; int b_y = b/hsize;
		return p_privMetric->get_priv_distance(p_xspace.cell_side*a_x, p_xspace.cell_side*a_y, p_xspace.cell_side*b_x, p_xspace.cell_side*b_y);
	};

	qif::Metric<double,uint> loss = [=](uint a, uint b) {
		int a_x = a % hsize; int a_y = a/hsize;
		int b_x = b % hsize; int b_y = b/hsize;
		return p_lossMetric->get_loss(p_xspace.cell_side*a_x, p_xspace.cell_side*a_y, p_xspace.cell_side*b_x, p_xspace.cell_side*b_y);
	};

	uint n_secrets = p_xspace.hsize * p_xspace.vsize;
	arma::rowvec pi(p_prior.get_prior_vector().data, p_prior.get_prior_vector().getLength(), true, true);

	// qif function to compute the mechanism matrix
	arma::mat C = qif::mechanism::d_privacy::min_loss_given_d(pi,n_secrets,priv,loss);

	// convert C (in arma::mat format) to DoubleMatrix format - false because C is stored column by column
	mechanism_matrix = DoubleMatrix(C.n_rows, C.n_cols, C.memptr(), false);

	// Important: C.memptr() returns a pointer to the data block of C and the data of C are stored column by column (column major format).
	// However DoubleMatrix stores data row by row (row major format).

	mechanism_file = p_mechFile;
	write_mechanism_to_file(p_mechFile);

	std::cout << "Expected loss = " << get_exp_loss_for_prior(p_prior, p_lossMetric, p_xspace.cell_side) << "\n";

	//std::cout << "Expected loss = " << qif::utility::expected_distance(loss, pi, C) << "\n";
	// << "C:\n" << C
}



