#include "objects.h"

SEXP compute_density (SEXP coords, SEXP centers, SEXP clust_info, SEXP radius) {
    BEGIN_RCPP

    const double rad=check_numeric_scalar(radius, "radius");
    auto searcher=generate_holder(coords, centers, clust_info);
    const size_t nhyper=searcher->get_ncells();

    Rcpp::NumericVector output(nhyper);
    for (size_t h=0; h<nhyper; ++h) {
        searcher->find_neighbors(h, rad, true);
        const std::deque<double>& distances=searcher->get_distances();

        double& curdensity=(output[h]=0);
        for (const auto& d : distances) { 
            const double ratio=d/rad;
            const double diffdist = 1 - ratio*ratio*ratio;
            curdensity += diffdist * diffdist * diffdist; // tricube weights.
        }
    }
    return output;
    END_RCPP
}

