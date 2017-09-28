#include "objects.h"
#include "packer.h"

SEXP count_cells(SEXP exprs, SEXP distance, SEXP centers, SEXP cluster_info, SEXP curcells) {
    BEGIN_RCPP
    auto searcher=generate_holder(exprs, centers, cluster_info);

    // Checking distances and chosen cells.
    const double threshold=check_numeric_scalar(distance, "distance");
    const Rcpp::IntegerVector _curcells(curcells);
    const int N=_curcells.size();
    const int ncells=searcher->get_ncells();
    const size_t& nmarkers=searcher->get_nmarkers();
    if (nmarkers==0) {
        throw std::runtime_error("number of markers should be positive");
    }
   
    // Setting up output vectors. 
    Rcpp::List outassign(N);
    Rcpp::IntegerVector outtotal(N);
    std::deque<int> sorted_ids;

    // Running through all cells.
    for (int ix=0; ix<N; ++ix) {
        const int& current_cell=_curcells[ix];
        if (current_cell >= ncells || current_cell < 0)  {
            throw std::runtime_error("chosen indices out of range");
        }

        searcher->find_neighbors(current_cell, threshold, false);
        std::deque<size_t>& collected=searcher->get_neighbors();
        if (collected.size()==0) {
            // Check here, otherwise median calculations fail.
            throw std::runtime_error("cell failed to count itself");
        }
            
        // Storing the identities of the cells (compressed).
        pack_index_vector(sorted_ids, collected.begin(), collected.end());
        outassign[ix]=Rcpp::IntegerVector(sorted_ids.begin(), sorted_ids.end());
        outtotal[ix]=collected.size();
    }

    return Rcpp::List::create(outassign, outtotal);
    END_RCPP
}

