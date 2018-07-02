#include "cydar.h"

SEXP drop_redundant (SEXP ordering, SEXP center_id, SEXP assignments) {
    BEGIN_RCPP

    const Rcpp::List Assignments(assignments);
    const int ngroups=Assignments.size();

    const Rcpp::IntegerVector center_cell(center_id);
    if (center_cell.size()!=ngroups) {
        throw std::runtime_error("length of 'center_id' is not equal to the number of groups");
    }

    const Rcpp::IntegerVector Ordering(ordering);
    if (Ordering.size()!=ngroups) {
        throw std::runtime_error("length of 'ordering' is not equal to the number of groups");
    }

    // Looking for hypersphers that are not redundant to hyperspheres with lower p-values.
    Rcpp::LogicalVector output(ngroups);
    std::deque<bool> already_seen(ngroups, false);

    for (auto o : Ordering) {
        if (already_seen[center_cell[o]-1]) { 
            continue; 
        }
        output[o]=1; // yes, 'o' is deliberate here; we're referring to the index of 'assignments'.

        const Rcpp::IntegerVector neighbors=Assignments[o];
        for (const auto& neigh : neighbors) { 
            already_seen[neigh-1] = true;
        }
    }

    return output;
    END_RCPP
}

