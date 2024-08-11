#include "Rcpp.h"

// [[Rcpp::export(rng=false)]]
Rcpp::NumericVector compute_density (Rcpp::List distances, double radius) {
    Rcpp::NumericVector output(distances.size());

    for (size_t i = 0, end = distances.size(); i < end; ++i) {
        Rcpp::NumericVector current_distances = distances[i];
        double curdensity = 1; // account for self, i.e., distance of zero.

        for (const auto& d : current_distances) { 
            const double ratio = d / radius;
            const double diffdist = 1 - ratio*ratio*ratio;
            curdensity += diffdist * diffdist * diffdist; // tricube weights.
        }

        output[i] = curdensity;
    }
    return output;
}

