#include "cydar.h"
#include "utils.h"

SEXP drop_redundant (SEXP intensities, SEXP ordering, SEXP assignments, SEXP threshold) {
    BEGIN_RCPP

    const Rcpp::List Assignments(assignments);
    const int ngroups=Assignments.size();

    const Rcpp::IntegerVector Ordering(ordering);
    if (Ordering.size()!=ngroups) {
        throw std::runtime_error("length of 'ordering' is not equal to the number of groups");
    }

    const Rcpp::NumericMatrix Intensities(intensities);
    if (ngroups!=Intensities.ncol()) {
        throw std::runtime_error("length of 'ordering' is not equal to number of columns in 'intensities'");
    }
    const int nmarkers=Intensities.nrow();

    const double Threshold=check_numeric_scalar(threshold, "threshold");

    // Looking for points that are not redundant to points with lower p-values.
    Rcpp::LogicalVector output(ngroups);
    std::deque<bool> already_seen(ngroups, false);

    for (auto o : Ordering) {
        if (already_seen[o]) { 
            continue; 
        }
        output[o]=1; 
        auto curint=Intensities.column(o);

        const Rcpp::IntegerVector neighbors=Assignments[o];
        for (const auto& neigh : neighbors) {
            auto neighint=Intensities.column(neigh - 1);
            bool is_within=true;

            for (int m=0; m<nmarkers; ++m) {
                if (std::abs(neighint[m] - curint[m]) > Threshold) {
                    is_within=false;
                    break;
                }
            }

            if (is_within) {
                already_seen[neigh-1] = true;
            }
        }
    }

    return output;
    END_RCPP
}

