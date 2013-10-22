#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>
#include "acor.h"
#include "exception.h"
#include "int2str.h"
#include "mean.h"
#include "model.h"
#include "quicksort.h"
#include "preprocessors.h"
#include "rng.h"
#include "var.h"

using namespace std;

double CovInterval(const Model & model, size_t i, size_t j);

void refine (Model & model) {
	
	model.level.LnPrims_r.resize(0);
	model.level.LnPrims_r.push_back(0);
	double Cx = 0;   // This number determines "how much" of the level building part is to be "trusted"
	for (size_t i = model.level.LnPrims_r.size(); i < model.level.num_level; ++i) {
		model.level.LnPrims_r.push_back( model.level.LnPrims_r[i-1] + log ((static_cast<double>(model.Above_Visits[i-1]) + Cx * exp(model.level.LnPrims[i] - model.level.LnPrims[i-1])) / (static_cast<double>(model.Level_Visits[i-1]) + Cx)) );
	}
	
	//for (size_t i = 0; i < model.level.num_level-1; ++i) {
		//cout << "refined   " << model.level.LnPrims_r[i] << endl;
		//cout << "unrefined " << model.level.LnPrims[i] << endl;
		//cout << model.Level_Visits[i] << endl;
		//cout << model.Above_Visits[i] << endl;
	//}
}


// This routine calculates the standard deviation of primor mass intervals.
void stdev_pm (Model & model, const double tau) {
	
	model.level.LnPrims_var.clear();
	model.level.LnPrims_var.push_back(-1.7e308);
	
	double var_xy, var_x, var_y, e_x, e_y;
	double p;
	for (size_t i = 1; i < model.level.num_level; ++i) {
		p = static_cast<double>(model.Above_Visits[i-1])/static_cast<double>(model.Level_Visits[i-1]);
		e_y   = p;
		e_x   = 1.0;               // scaled by exp(model.level.LnPrims_r[i-1])
		var_y = p * (1-p) / static_cast<double>(model.Level_Visits[i-1]) * tau;
		var_x = exp( model.level.LnPrims_var[i-1] - 2.0 * model.level.LnPrims_r[i-1] );
		var_xy = var_x * var_y + var_x * e_y * e_y + var_y * e_x * e_x;
		model.level.LnPrims_var.push_back( log(var_xy) + 2.0 * model.level.LnPrims_r[i-1] );
	}
	
	model.level.LnInterval_var.clear();
	for (size_t i = 0; i < model.level.num_level-1; ++i) {
		p = static_cast<double>(model.Above_Visits[i])/static_cast<double>(model.Level_Visits[i]);
		e_y   = 1. - p;
		e_x   = 1.0;                // scaled by exp(model.level.LnPrims_r[i]);
		var_y = p * (1-p) / static_cast<double>(model.Level_Visits[i]) * tau;
		var_x = exp( model.level.LnPrims_var[i] - 2.0 * model.level.LnPrims_r[i-1] );
		var_xy = var_x * var_y + var_x * e_y * e_y + var_y * e_x * e_x;
		model.level.LnInterval_var.push_back( log(var_xy) + 2.0 * model.level.LnPrims_r[i-1] );
	}
	model.level.LnInterval_var.push_back(model.level.LnPrims_var.back());
	
	for (size_t i = 0; i < model.level.num_level; ++i) {
		cout << model.level.LnPrims_var[i] << "    " << model.level.LnInterval_var[i] << endl;
	}
}
		
double evidence (Model & model, size_t ens_size) {
	
	model.level.LnThres.push_back(1.7e308);
	model.level.LnPrims.push_back(-1.7e308);
	vector< vector<double> > gaps(model.level.num_level);
	for (size_t i = 0; i < model.chain_LnLikelihood.size(); ++i) {
		for (size_t j = 0; j < model.level.num_level; ++j) {
			if ( model.chain_LnLikelihood[i] >= model.level.LnThres[j] && model.chain_LnLikelihood[i] < model.level.LnThres[j+1] ) {
				gaps[j].push_back(model.chain_LnLikelihood[i]);
				break;
			}
		}
	}
	
	double mean, sigma, tau;
	try {
		acor(mean, sigma, tau, model.Quant_Ens_Mean, model.Quant_Ens_Mean.size());
	}
	catch( Exception & e ){
		cout << e.ExceptionMessage() << endl;
		fstream acorout;
		acorout.open(("acor_err_"+model.model_name+"_"+model.time_label+"_"+int2str(model.Quant_Ens_Mean.size(),11)+".txt").c_str(), ios::out);
		for (size_t i = 0; i < model.Quant_Ens_Mean.size(); ++i) {
			acorout << model.Quant_Ens_Mean[i] << endl;
		}
		acorout.close();
		tau = 1;
	}
	refine(model);
	stdev_pm(model, tau);
	
	double integral = 0;
	double integral_r = 0;
	double integral_err1 = 0;   // caused by uncertainty of the integrand
	double integral_err2 = 0;   // caused by uncertainty of the prior masses
	double integral_err3 = 0;   // caused by uncertainty of the prior mass cross terms
	vector< double > lml(gaps.size());
	vector< double > lvl(gaps.size());
	vector< double > Lvar(gaps.size());
	for (size_t i = 0; i < lml.size(); ++i) {
		if (gaps[i].size() == 0) {
			lml[i] = -1e300;
			cout << i << " empty!" << endl;
			continue;
		}
		lml[i] = logMeanLog(gaps[i]);
		lvl[i] = logVarLog(gaps[i]);
		Lvar[i] = lvl[i] + log(tau) - log(static_cast<double>(gaps[i].size()));
	}
	double max_lml = *max_element( lml.begin(), lml.end() );
	for (size_t j = 0; j < model.level.num_level; ++j) {
		//cout << j << endl;
		if (j < model.level.num_level-1) {
			integral   += exp( lml[j]-max_lml ) * ( exp( model.level.LnPrims[j]-model.level.LnPrims[j+1] ) - 1 ) * exp( model.level.LnPrims[j+1] );
			integral_r += exp( lml[j]-max_lml ) * ( exp( model.level.LnPrims_r[j]-model.level.LnPrims_r[j+1] ) - 1 ) * exp( model.level.LnPrims_r[j+1] );
			integral_err1 += exp(Lvar[j] - 2.0 * max_lml) * (exp(model.level.LnPrims_r[j]-model.level.LnPrims_r[j+1]) - 1) * exp(model.level.LnPrims_r[j+1]) * (exp(model.level.LnPrims_r[j]-model.level.LnPrims_r[j+1]) - 1) * exp(model.level.LnPrims_r[j+1]);
			integral_err2 += exp(2.0 * (lml[j]-max_lml)) * exp(model.level.LnInterval_var[j]);
		}
		else {
			integral += exp(lml[j]-max_lml) * (exp(model.level.LnPrims[j]) - 0);
			integral_r += exp(lml[j]-max_lml) * (exp(model.level.LnPrims_r[j]) - 0);
			integral_err1 += exp(Lvar[j] - 2.0 * max_lml) * (exp(model.level.LnPrims_r[j]) - 0) * (exp(model.level.LnPrims_r[j]) - 0);
			integral_err2 += exp(2.0 * (lml[j]-max_lml)) * exp(model.level.LnInterval_var[j]);
		}
	}
	
	for (size_t i = 0; i < model.level.num_level - 1; ++i) {
		for (size_t j = i + 1; j < model.level.num_level; ++j) {
			integral_err3 += 2.0*exp( lml[i] + lml[j] - 2.0*max_lml ) * CovInterval(model, i, j) ;
			cout << "i = " << i << " j = " << j << " Cov = " << CovInterval(model, i, j) << endl;
		}
	}
	
	cout << "log unrefined evi = " << log(integral)+max_lml << endl;
	cout << "log   refined evi = " << log(integral_r)+max_lml << endl;
	cout << "err1  " << exp(log(integral_err1)+2.0*max_lml) << endl; 
	cout << "err2  " << exp(log(integral_err2)+2.0*max_lml) << endl; 
	cout << "err3  " << exp(log(integral_err3)+2.0*max_lml) << endl; 
	model.evidence = log(integral)+max_lml;
	model.evidence_r = log(integral_r)+max_lml;
	model.evidence_err1 = log(integral_err1)+2.0*max_lml;
	model.evidence_err2 = log(integral_err2)+2.0*max_lml;
	model.evidence_err3 = log(integral_err3)+2.0*max_lml;
	
	fstream out;
	//////////////////////////////////////////////////////
	out.open(("refines_"+model.model_name+"_"+model.time_label+".txt").c_str(), ios::out);
	//////////////////////////////////////////////////////
	out << setprecision(17) << "log unrefined evi = " << log(integral)+max_lml << endl;
	out << setprecision(17) << "log   refined evi = " << log(integral_r)+max_lml << endl;
	out << setprecision(17) << "log err1          = " << log(integral_err1)+2.*max_lml << endl; 
	out << setprecision(17) << "log err2          = " << log(integral_err2)+2.*max_lml << endl;
	out << setprecision(17) << "log err3          = " << log(integral_err3)+2.*max_lml << endl; 
	out << setprecision(17) << "autocorrelation t = " << tau << endl; 
	out << "Level, LnThres, LnPrims, LnPrims_r, Visits, Lvar, Above" << endl;
	for (size_t j = 0; j < model.level.num_level; ++j) {
		out << j << "   ";
		out << setprecision(17) << model.level.LnThres[j] << "   ";
		out << setprecision(17) << model.level.LnPrims[j] << "   ";
		out << setprecision(17) << model.level.LnPrims_r[j] << "   ";
		out << setprecision(17) << model.Level_Visits[j] << "   ";
		out << setprecision(17) << Lvar[j] << "   ";
		if (j < model.level.num_level - 1) {
			out << setprecision(17) << model.Above_Visits[j] << "   ";
		}
		out << endl;
	}
	model.level.LnThres.pop_back();
	model.level.LnPrims.pop_back();
	return model.evidence_r;
}

double CovInterval(const Model & model, size_t i, size_t j) {
	if ( i == j ) {
		return 0;
	}
	if ( i > j ) {
		double temp = j;
		j = i;
		i = temp;
	}
	double Cov = 0;
	if ( i == 0 ) {
		Cov = exp( model.level.LnPrims_var[i+1] - model.level.LnPrims_r[i+1] ) * ( - 1.0 + exp ( model.level.LnPrims_r[j+1] - model.level.LnPrims_r[j] ) ) * exp( model.level.LnPrims_r[j] );
		return Cov;
	}
	if ( j == model.level.num_level - 1 ) {
		Cov = ( exp( model.level.LnPrims_var[i] - model.level.LnPrims_r[i] ) - exp( model.level.LnPrims_var[i+1] - model.level.LnPrims_r[i+1] ) ) * exp( model.level.LnPrims_r[j] );
		return Cov;
	}
	Cov = ( exp(model.level.LnPrims_var[i] - model.level.LnPrims_r[i]) - exp(model.level.LnPrims_var[i+1] - model.level.LnPrims_r[i+1]) ) * ( 1.0 - exp ( model.level.LnPrims_r[j+1] - model.level.LnPrims_r[j] )) * exp( model.level.LnPrims_r[j] );
	return Cov;
}
