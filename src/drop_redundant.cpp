#include "Rcpp.h"
#include <cmath>
#include <stdexcept>
#include <deque>

// [[Rcpp::export(rng=false)]]
Rcpp::LogicalVector drop_redundant (Rcpp::NumericMatrix intensities, Rcpp::IntegerVector ordering, 
    Rcpp::List assignments, double threshold) 
{
    const int ngroups=assignments.size();
    if (ordering.size()!=ngroups) {
        throw std::runtime_error("length of 'ordering' is not equal to the number of groups");
    }

    const int nmarkers=intensities.nrow();
    if (ngroups!=intensities.ncol()) {
        throw std::runtime_error("length of 'ordering' is not equal to number of columns in 'intensities'");
    }

    // Looking for points that are not redundant to points with lower p-values.
    Rcpp::LogicalVector output(ngroups);
    std::deque<bool> already_seen(ngroups, false);

    for (auto o : ordering) {
        if (already_seen[o]) { 
            continue; 
        }
        output[o]=1; 
        auto curint=intensities.column(o);

        const Rcpp::IntegerVector neighbors=assignments[o];
        for (const auto& neigh : neighbors) {
            auto neighint=intensities.column(neigh - 1);
            bool is_within=true;

            for (int m=0; m<nmarkers; ++m) {
                if (std::abs(neighint[m] - curint[m]) > threshold) {
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
}

