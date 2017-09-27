#include "cydar.h"
#include "packer.h"
#include "utils.h"

SEXP pack_indices(SEXP assignments, SEXP compact) {
    BEGIN_RCPP

    const Rcpp::List _assignments(assignments);
    const int ngrps=LENGTH(assignments);
    const bool compress=check_logical_scalar(compact, "compact specifier");

    Rcpp::List output(ngrps);
    std::deque<int> sorted_ids, temp;

    for (int g=0; g<ngrps; ++g) {
        const Rcpp::IntegerVector current(_assignments[g]);
            
         if (compress) {
            temp.assign(current.begin(), current.end());
            for (auto& t : temp) { --t; } // getting to zero-indexing.
            pack_index_vector(sorted_ids, temp.begin(), temp.end());
            output[g]=Rcpp::IntegerVector(sorted_ids.begin(), sorted_ids.end());

        } else {
            unpack_index_vector(temp, current.begin(), current.end());
            output[g]=Rcpp::IntegerVector(temp.begin(), temp.end());
        }
    }

    return output;
    END_RCPP
}

