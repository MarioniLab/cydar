#include "Rcpp.h"

// [[Rcpp::export(rng=false)]]
Rcpp::NumericVector compute_density (Rcpp::List distances, double radius) {
    Rcpp::NumericVector output(distances.size());

    for (size_t i=0; i<distances.size(); ++i) {
        Rcpp::NumericVector current_distances=distances[i];
        double& curdensity=(output[i]=0);
        for (const auto& d : current_distances) { 
            const double ratio=d/radius;
            const double diffdist = 1 - ratio*ratio*ratio;
            curdensity += diffdist * diffdist * diffdist; // tricube weights.
        }
    }
    return output;
}

