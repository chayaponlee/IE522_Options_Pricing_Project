#include <iostream>
#include <iomanip>
#include <chrono>
#include "asian_opt.hpp"

int main() {
    // Option parameters
    double S0 = 100.0;
    double K = 100.0;
    double r = 0.10;
    double q = 0.00;
    double sigma = 0.20;
    double T = 1.0;
    int m = 50;

    // Sample sizes to investigate convergence
    int sample_sizes[] = { 1000, 4000, 16000, 64000, 128000 };

    // Price the option for different sample sizes
    for (int N : sample_sizes) {
        auto start_time = std::chrono::high_resolution_clock::now();

        AsianCallOption option(S0, K, r, q, sigma, T, m);
        double standard_error = 0.0;
        double price = option.price(N, standard_error);

        auto end_time = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> computational_time = end_time - start_time;

        // 95% confidence interval
        double margin = 1.96 * standard_error;
        double ci_lower = price - margin;
        double ci_upper = price + margin;

        // Output the results
        std::cout << "Sample Size: " << N << std::endl;
        std::cout << "Option Price: " << std::fixed << std::setprecision(4) << price << std::endl;
        std::cout << "Estimated Standard Error: " << standard_error << std::endl;
        std::cout << "95% Confidence Interval: [" << ci_lower << ", " << ci_upper << "]" << std::endl;
        std::cout << "Computational Time (seconds): " << computational_time.count() << std::endl;
        std::cout << "-------------------------------------------" << std::endl;
    }

    return 0;
}
