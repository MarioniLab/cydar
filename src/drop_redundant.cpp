#include "objects.h"

SEXP drop_redundant (SEXP actual_order, SEXP coords, SEXP centers, SEXP clust_info, SEXP threshold) {
    BEGIN_RCPP

    Rcpp::NumericMatrix _coords(coords);
    auto searcher=generate_holder(coords, centers, clust_info);
    const size_t nhyper=searcher->get_ncells();
    const size_t nmarkers=searcher->get_nmarkers();

    const double thresh=check_numeric_scalar(threshold, "threshold");
    const double radius=thresh * std::sqrt(nmarkers);

    Rcpp::IntegerVector ordering(actual_order);
    if (ordering.size()!=nhyper) {
        throw std::runtime_error("length of actual_order order vector must be equal to number of hyperspheres");
    }
    
    // Looking for hypersphers that are not redundant to hyperspheres with lower p-values.
    Rcpp::LogicalVector output(nhyper);
    std::deque<bool> already_seen(nhyper, false);

    for (const auto& actual_index : ordering) { 
        if (already_seen[actual_index]) { continue; }
        output[actual_index]=1;

        searcher->find_neighbors(actual_index, radius, false);
        const std::deque<size_t>& neighbors=searcher->get_neighbors();
        auto curcoords=_coords.column(actual_index);

        for (const auto& neigh : neighbors) { 
            bool okay=false;
            auto othercoords=_coords.column(neigh);

            for (size_t mi=0; mi<nmarkers; ++mi) {
                if (std::abs(curcoords[mi] - othercoords[mi]) > thresh) {
                    okay=true;
                    break;
                }
            }
            if (!okay) { 
                already_seen[neigh] = true;
            }
        }
    }

    return output;
    END_RCPP
}

