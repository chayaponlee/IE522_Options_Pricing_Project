#include "asian_call_option.hpp"
#include <cmath>
#include <random>
#include <iostream>

AsianCallOption::AsianCallOption(double S0, double K, double r, double q, double sigma, double T, double m) :
	S0(S0), K(K), r(r), q(q), sigma(sigma), T(T), m(m) {}

double AsianCallOption::price(int N, double& standard_error) {
	double payoff_sum = 0.0;
	double payoff_squared_sum = 0.0;

	for (int i = 0; i < N; i++) { // N: number of sample paths
		double payoff = std::max(0.0, sample_path_generator() - K);
		payoff_sum += payoff;
		payoff_squared_sum += payoff * payoff;
	}
	double discount_factor = std::exp(-r * T);
	double option_price = discount_factor * (payoff_sum / N);

	// calculate standard error
	double variance = (payoff_squared_sum - (payoff_sum * payoff_sum) / N) / (N - 1);
	standard_error = discount_factor * std::sqrt(variance / N);
	return option_price;
}

double AsianCallOption::sample_path_generator() {
	double dt = T / m;
	double drift = (r - q - 0.5 * sigma * sigma) * dt;
	double diffusion = sigma * std::sqrt(dt);
	double sum_S = 0.0;
	double S = S0;


	for (int i = 1; i <= m; i++) {
		double z = gaussian_box_muller();
		S *= std::exp(drift + diffusion * z);
		sum_S += S;
	}
	return sum_S / m;
}

double AsianCallOption::gaussian_box_muller() {
	static std::mt19937 rng(std::random_device{}());
	static std::normal_distribution<> nd(0.0, 1.0);
	return nd(rng);
}