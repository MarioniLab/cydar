#include "cydar.h"
#include "utils.h"

SEXP compute_density (SEXP distances, SEXP radius) {
    BEGIN_RCPP
    const double rad=check_numeric_scalar(radius, "radius");
    Rcpp::List Distances(distances);
    Rcpp::NumericVector output(Distances.size());

    for (size_t i=0; i<Distances.size(); ++i) {
        Rcpp::NumericVector current_distances=Distances[i];
        double& curdensity=(output[i]=0);
        for (const auto& d : current_distances) { 
            const double ratio=d/rad;
            const double diffdist = 1 - ratio*ratio*ratio;
            curdensity += diffdist * diffdist * diffdist; // tricube weights.
        }
    }
    return output;
    END_RCPP
}

