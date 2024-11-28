#ifndef ASIAN_OPT_HPP
#define ASIAN_OPT_HPP

class AsianCallOption {
public:
	AsianCallOption(double S0, double K, double r, double q, double sigma, double T, double m); // constructor

	double price(int N, double& std_err);

private:
	double S0;
	double K;
	double r;
	double q;
	double sigma;
	double T;
	double m;

	double sample_path_generator();
	double gaussian_box_muller(); // random number generator
};


#endif ASIAN_OPT_HPP