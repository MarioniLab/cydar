#include "objects.h"

SEXP get_nndist(SEXP cells, SEXP centers, SEXP clust_info, SEXP nn, SEXP freq) {
    BEGIN_RCPP

    const size_t NN=check_integer_scalar(nn, "number of nearest neighbors");
    const size_t downsample=check_integer_scalar(freq, "downsampling frequency");
    auto searcher=generate_holder(cells, centers, clust_info);
    const size_t ncells=searcher -> get_ncells();

    const size_t true_nrows=(ncells ? 1+int((ncells-1)/downsample) : 0);
    Rcpp::NumericMatrix output(NN, true_nrows);
    auto oIt=output.begin();
        
    for (size_t h=0; h<ncells; h+=downsample) {
        searcher->find_nearest_neighbors(h, NN, true); 
        const std::deque<double>& distances=searcher->get_distances();
        std::copy(distances.begin(), distances.end(), oIt);
        oIt+=NN;
    }

    return output;
    END_RCPP
}

