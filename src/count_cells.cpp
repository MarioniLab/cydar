#include "Rcpp.h"
#include <stdexcept>

// [[Rcpp::export(rng=false)]]
Rcpp::IntegerMatrix count_cells(Rcpp::List assignments, Rcpp::NumericVector sample_id, int nsamples) {
    const int ngroups=assignments.size();
    const int ncells=sample_id.size();
    for (const auto& cursample : sample_id) { 
        if (cursample < 0 || cursample >= nsamples) {
            throw std::runtime_error("sample IDs out of range");
        }
    }

    Rcpp::IntegerMatrix outcounts(nsamples, ngroups);
    for (int g=0; g<ngroups; ++g) {
        const Rcpp::IntegerVector curass=assignments[g];
        for (auto curdex : curass) {
            if (curdex <= 0 || curdex > ncells) {
                throw std::runtime_error("cell assignment indices out of range");
            }
        }
            
        // Computing counts and total weights.
        auto curcounts=outcounts.column(g);
        for (const auto& c : curass) { 
            const int& cursample=sample_id[c-1];
            ++(curcounts[cursample]);
        }
    }

    return outcounts;
}
