#include <cmath>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include "acor.h"
#include "data.h"
#include "evidence.h"
#include "exoplanetjd.h"
#include "int2str.h"
#include "model.h"
#include "preprocessors.h"
#include "quicksort.h"
#include "rng.h"
#include "sampling_DNest.h"
#include "uniformtest.h"


using namespace std;

void read_levels_from_file (Model & model, string level_file_name);
void write_levels_to_file (Model & model);
void write_levels_to_file_obo (Model & model);
void estimate_evidence (Model & model, double & upper, double & lower, double & med);

int main(void) {
	time_t begin = time(NULL);
	init_genrand(begin);
	time_t end;
	
	//string data_file_name("ploy_3_151_15_-12_9_-1"); size_t dim = 4;
	//string data_file_name("ploy_1_10_5_10"); size_t dim = 2;
	//string data_file_name("uni_2"); size_t dim = 2; size_t ens_size = 20; int num_level = 10; size_t step_size = 10;
	//string data_file_name("uni_10"); size_t dim = 10; size_t ens_size = 50; int num_level = 40; size_t step_size = 100;
	//string data_file_name("122"); size_t num_comp = 1; size_t ens_size = 150; int num_level = 100; size_t step_size = 10;
	//string data_file_name("122"); size_t num_comp = 2; size_t ens_size = 200; int num_level = 150; size_t step_size = 10;
	//string data_file_name("122"); size_t num_comp = 3; size_t ens_size = 200; int num_level = 150; size_t step_size = 10;
	//string data_file_name("122"); size_t num_comp = 5; size_t ens_size = 200; int num_level = 220; size_t step_size = 10;
	//string data_file_name("gliese_581"); size_t num_comp = 6; size_t ens_size = 200; int num_level = 0; size_t step_size = 10;
	//string data_file_name("gliese_581"); size_t num_comp = 5; size_t ens_size = 200; int num_level = 0; size_t step_size = 10;
	//string data_file_name("gliese_581"); size_t num_comp = 4; size_t ens_size = 200; int num_level = 0; size_t step_size = 10;
	//string data_file_name("gliese_581"); size_t num_comp = 3; size_t ens_size = 200; int num_level = 0; size_t step_size = 10;
	//string data_file_name("gliese_581"); size_t num_comp = 7; size_t ens_size = 200; int num_level = 150; size_t step_size = 10;
	
	//string data_file_name("x001"); size_t num_comp = 1; size_t ens_size = 50; int num_level = 150; size_t step_size = 10;
	string data_file_name("x001"); size_t num_comp = 2; size_t ens_size = 50; int num_level = 150; size_t step_size = 10;
	//string data_file_name("282"); size_t num_comp = 3; size_t ens_size = 50; int num_level = 150; size_t step_size = 10;
	//string data_file_name("282"); size_t num_comp = 4; size_t ens_size = 200; int num_level = 0; size_t step_size = 10;
	
	//string data_file_name("fake4"); size_t num_comp = 3; size_t ens_size = 200; int num_level = 150; size_t step_size = 10;
	Data data(data_file_name);
	string weight_type = "uniform";
	ExoplanetJD model(data, num_comp, weight_type);
	//UniformTest model(data, dim, weight_type); size_t num_comp = dim/2;
	matrix ensemble;
	model.init(ens_size, ensemble, 0.0000);
	cout << model.LnLikelihood(ensemble[0]) << endl; 
	cout << model.LnDensity(ensemble[0]) << endl;
	vector<size_t> levels(ens_size, 0);
	
	double ini = 0.00000001;
	size_t build_burn_in_steps = 5000;
	size_t build_chain_length  = 147513;
	size_t DNest_chain_length  = 100000;
	size_t num_loops = 200;
	bool succeed;
	
	fstream eout;
	//eout.open(("evi_"+model.model_name+"_"+int2str(time(NULL),0)+".txt").c_str(),ios::out);
	time_t b = time(NULL);
	int num_samples = 1;  // number of evidences to calculate
	bool redo = 0;
	bool redone = 0;
	for (int l = 0; l < num_samples; ++l) {
		cout << l << endl;
		model.level.clear_all();
		levels.clear();
		levels.resize(ens_size, 0);
		
		
		//read_levels_from_file(model, "output_levels_of_exop_gliese_581_mod_3_d_1366940332.txt");
		//read_levels_from_file(model, "output_levels_of_exop_gliese_581_mod_4_d_1366943322.txt");
		//read_levels_from_file(model, "output_levels_of_exop_gliese_581_mod_5_d_1366946223.txt");
		//read_levels_from_file(model, "output_levels_of_exop_gliese_581_mod_6_d_1366948780.txt");
		//read_levels_from_file(model, "output_levels_of_exop_282_mod_1_1373315182.txt");
		//read_levels_from_file(model, "output_levels_of_exop_282_mod_2_1366934165.txt");
		//read_levels_from_file(model, "output_levels_of_exop_282_mod_3_1366935956.txt");
		//read_levels_from_file(model, "output_levels_of_exop_282_mod_4_1366937633.txt");

		
		for (int i = model.level.num_level; i <= num_level; ++i) {
			model.init(ens_size, ensemble, ini);
			levels.clear();
			levels.resize(ens_size, model.level.num_level-1);
			
			
			succeed = new_level(model, ensemble, levels, 1.5, build_burn_in_steps, build_chain_length, build_chain_length*2, step_size);
			
			if (succeed == 0) {
				ini *= 0.01;
				--i;
				cout << "failed " << "new ini = " << ini << endl;
				continue;
			}
			cout << model.level.LnThres.size() << endl;
			cout << setprecision(15) << model.level.LnThres[i] << endl;
			cout << setprecision(15) << model.level.LnPrims[i] << endl;
			write_levels_to_file_obo(model);
			time_t used = time(NULL) - b;
			double upper, lower, med;
			upper = lower = med = 0;
			estimate_evidence(model, upper, lower, med);
			
			if (model.best_fit[model.dim] + model.level.LnPrims.back() < log(1.e-7) + lower) {
				cout << "OK to stop in " << l <<  endl;
				break;
			}
			size_t n = model.level.num_level;
			if (0 && redone == 0 && n>=4 && model.level.LnThres[n-1] - model.level.LnThres[n-2] > model.level.LnThres[n-2] - model.level.LnThres[n-3]) {
				model.level.pop_one_level();
				redo = 1;
				--i;
				continue;
			}
			//static fstream varout(("lkvar_" + model.model_name + "_" + model.time_label + ".txt").c_str(), ios::out);
			//varout << "between " << i-1 << " and " << i << "    " << setprecision(20) << model.level.LnLikelihood_var[i-1] << endl;
		}
		
		double upper, lower, med;
		upper = lower = med = 0;
		estimate_evidence(model, upper, lower, med);
		cout << "upper limit = " << upper << endl;
		cout << "lower limit = " << lower << endl;
		cout << "med   limit = " << med << endl;
		
		
		model.init(ens_size, ensemble, ini*10);
		vector<size_t> levels(ens_size, static_cast<size_t>(genrand_real2()*(model.level.num_level-1)));
		//write_levels_to_file(model);
		
		double evi;
		model.clear_all_visits();
		model.Level_Visits_T.resize(model.level.num_level);
		model.Level_Visits.resize(model.level.num_level);
		model.Above_Visits.resize(model.level.num_level-1);
		sampling_DNest(model, num_comp*40000+10000, step_size, ensemble, levels, 1.5);
		fstream acorout;
		acorout.open(("acor_burn_"+model.model_name+"_"+model.time_label+"_"+int2str(model.Quant_Ens_Mean.size(),11)+".txt").c_str(), ios::out);
		for (size_t i = 0; i < model.Quant_Ens_Mean.size(); ++i) {
			acorout << model.Quant_Ens_Mean[i] << endl;
		}
		acorout.close();
		
		model.clear_for_sampling();
		model.Level_Visits.resize(model.level.num_level);
		model.Above_Visits.resize(model.level.num_level-1);
		for (size_t i = 0; i < num_loops; ++i) {
			sampling_DNest(model, DNest_chain_length, step_size, ensemble, levels, 1.5);
			evi = evidence(model, ens_size); // refinement done inside the routine evidence
		}
		//cout << evi << endl;
		eout << setprecision(16) << model.evidence << "    " << model.evidence_r << "    " << model.evidence_err2 << endl;
		time_t e = time(NULL);
		cout << (double)(e - b)/(double)(l+1)*(double)(num_samples-l-1) << " seconds left!" << endl;
		
	
	}


	cout << "Time Cost: " << time(NULL) - b << " seconds!" << endl;
	
	
	
	return 0;
}

void read_levels_from_file (Model & model, string level_file_name) {
	model.level.clear_all();
	fstream input;
	input.open(level_file_name.c_str(), ios::in);
	if(input.fail() == 1) {
		cerr << level_file_name << " doesn't exist!" << endl;
		return;
	}
	double temp_L, temp_X;
	
	while(1) {
		input >> temp_L;
		input >> temp_X;
		if(input.eof()) break;
		model.level.add_level(temp_L, temp_X);
		//cout << temp_L << endl;
	}
	
	cout << "Read Levels from " << level_file_name << endl;
	for (size_t i = 0; i < model.level.num_level; ++i) {
		cout << "Level " << i << " " << model.level.LnThres[i] << " " << model.level.LnPrims[i] << endl;
	}
}

void write_levels_to_file_obo (Model & model) {
	
	string out_file_name("output_levels_of_" + model.model_name + "_" + model.time_label + ".txt");
	static fstream out(out_file_name.c_str(), ios::out);

	out << setprecision(15) << *(--model.level.LnThres.end());
	out << "   ";
	out << setprecision(15) << *(--model.level.LnPrims.end());
	out << endl;
}

void write_levels_to_file (Model & model) {
	time_t label = time(NULL);
	stringstream time_label;
	time_label << static_cast<unsigned long>(label);
		//cout << time_label.str() << endl;
	string out_file_name("obo_levels_of_" + model.model_name + "_" + model.time_label + ".txt");
	fstream out(out_file_name.c_str(), ios::out);
	// start from level 1 because level 0 is trivial
	for (size_t i = 1; i < model.level.num_level; ++i) {
		out << setprecision(15) << model.level.LnThres[i];
		out << "   ";
		out << setprecision(15) << model.level.LnPrims[i];
		out << endl;
	}
	out.close();
}

void estimate_evidence (Model & model, double & upper, double & lower, double & med) {
	size_t num_level = model.level.num_level;
	double max_thres = model.level.LnThres[num_level-1];
	model.level.LnPrims.push_back(-1.7e308);
	double right, left;
	for (size_t i = 1; i < num_level; ++i) {
		right = (exp(model.level.LnPrims[i-1]) + exp(model.level.LnPrims[i])) * 0.5;
		left =  (exp(model.level.LnPrims[i]) + exp(model.level.LnPrims[i+1])) * 0.5;
		med += exp(model.level.LnThres[i]-max_thres) * (right - left);
		upper += exp(model.level.LnThres[i]-max_thres) * (exp(model.level.LnPrims[i-1]) - exp(model.level.LnPrims[i]));
		//cout << i << "    " << exp(model.level.LnThres[i]-max_thres) * (exp(model.level.LnPrims[i-1]) - exp(model.level.LnPrims[i])) << endl;
		lower += exp(model.level.LnThres[i]-max_thres) * (exp(model.level.LnPrims[i]) - exp(model.level.LnPrims[i+1]));
		//cout << i << "    " << exp(model.level.LnThres[i]-max_thres) * (right - left) << endl;
	}
	upper = log(upper) + max_thres;
	lower = log(lower) + max_thres;
	med = log(med) + max_thres;
	model.level.LnPrims.pop_back();
}
