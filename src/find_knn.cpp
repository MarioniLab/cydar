#include "objects.h"

SEXP find_knn(SEXP cells, SEXP clust_centers, SEXP clust_info, SEXP nn, SEXP mode, SEXP query) {
    BEGIN_RCPP

    const size_t NN=check_integer_scalar(nn, "number of nearest neighbors");
    if (NN<1) { 
        throw std::runtime_error("number of nearest neighbors must be positive");
    }
    auto searcher=generate_holder(cells, clust_centers, clust_info);
    const size_t nmarkers=searcher->get_nmarkers();

    // Just iterating across itself, if there are no query cells.
    const bool self_cycle=(query==R_NilValue);
    Rcpp::NumericMatrix _query;
    if (self_cycle) { 
        _query=Rcpp::NumericMatrix(cells);
    } else {
        _query=Rcpp::NumericMatrix(query);
        if (_query.nrow()!=nmarkers) {
            throw std::runtime_error("host and target intensity matrix do not have same number of markers");
        }
    }
    const size_t ncells=_query.ncol();

    // Getting the output mode.
    int M=check_integer_scalar(mode, "mode");
    const bool keep_last=(M < 0);
    if (keep_last) { M*=-1; }
    const bool store_distances=(M%2==1), store_neighbors=(M>=2);

    Rcpp::List output(2, R_NilValue);
    double* odptr=NULL;
    if (store_distances) { // Record all distances.
        if (keep_last) { 
            Rcpp::NumericVector current(ncells);
            output[1]=current;
            odptr=current.begin();
        } else {
            Rcpp::NumericMatrix current(NN, ncells);
            output[1]=current;
            odptr=current.begin();
        }
    }

    int* onptr=NULL;
    if (store_neighbors) { // Record all neighbours.
        if (keep_last) { 
            Rcpp::IntegerVector current(ncells);
            output[0]=current;
            onptr=current.begin();
        } else {
            Rcpp::IntegerMatrix current(NN, ncells);
            output[0]=current;
            onptr=current.begin();
        }
    }
        
    // Iterating across cells, finding NNs and storing (last) distances or neighbors.
    auto qIt=_query.begin();
    for (size_t h=0; h<ncells; ++h) {
        if (self_cycle) { 
            searcher->find_nearest_neighbors(h, NN, store_distances); 
        } else {
            searcher->find_nearest_neighbors(qIt, NN, store_distances); 
            qIt+=nmarkers;
        }
        const std::deque<double>& distances=searcher->get_distances();
        const std::deque<size_t>& neighbors=searcher->get_neighbors();

        if (store_distances) {
            if (keep_last) {
                (*odptr)=distances.back();
                ++odptr;
            } else {
                std::copy(distances.begin(), distances.end(), odptr);
                odptr+=NN;
            }
        }
        if (store_neighbors) {
            if (keep_last) { 
                (*onptr)=neighbors.back();
                ++onptr;
            } else {
                for (const auto& n : neighbors) { 
                    (*onptr)=n;
                    ++onptr;
                }
            }
        }
    }

    if (store_distances && store_neighbors) { 
        return output;
    } else if (store_distances) {
        return output[1];
    } else if (store_neighbors) {
        return output[0];
    } else {
        return R_NilValue;
    }
    END_RCPP
}


