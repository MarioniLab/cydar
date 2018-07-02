#include "cydar.h"

SEXP drop_redundant (SEXP center_id, SEXP assignments) {
    BEGIN_RCPP

    const Rcpp::List Assignments(assignments);
    const int ngroups=Assignments.size();

    const Rcpp::IntegerVector center_cell(center_id);
    if (center_cell.size()!=ngroups) {
        throw std::runtime_error("length of 'center_id' is not equal to the number of groups");
    }

    // Looking for hypersphers that are not redundant to hyperspheres with lower p-values.
    Rcpp::LogicalVector output(ngroups);
    std::deque<bool> already_seen(ngroups, false);

    for (size_t i=0; i<ngroups; ++i) {
        const int idx=center_cell[i] - 1;
        if (already_seen[idx]) { 
            continue; 
        }
        output[idx]=1;

        const Rcpp::IntegerVector neighbors=Assignments[i];
        for (const auto& neigh : neighbors) { 
            already_seen[neigh-1] = true;
        }
    }

    return output;
    END_RCPP
}

