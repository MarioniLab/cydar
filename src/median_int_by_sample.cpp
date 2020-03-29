#include "Rcpp.h"
#include <algorithm>
#include <stdexcept>
#include <vector>
#include <deque>

SEXP weighted_median_int(
    Rcpp::NumericMatrix exprs,
    Rcpp::List assignments,
    Rcpp::NumericVector sample_id,
    int nsamples) 
{
    const size_t nmarkers=exprs.ncol();
    const size_t ncells=exprs.nrow();
    const int ngroups=assignments.size();

    // Checking samples.
    if (nsamples <= 0) { 
        throw std::runtime_error("number of samples must be positive"); 
    }
    if (sample_id.size()!=ncells) { 
        throw std::runtime_error("sample IDs should be an integer vector of length equal to the number of cells"); 
    }
    for (const auto& cursample : sample_id) { 
        if (cursample <= 0 || cursample > nsamples) {
            throw std::runtime_error("sample IDs out of range");
        }
    }

    // Setting up output matrices (don't use List here, as it gets stored as a funny Proxy type).
    std::vector<Rcpp::NumericMatrix> output;
    output.reserve(nmarkers);
    for (size_t u=0; u<nmarkers; ++u) {
        output.push_back(Rcpp::NumericMatrix(ngroups, nsamples));
    }

    std::vector<std::deque<double> > all_intensities(nsamples);
    for (int g=0; g<ngroups; ++g) {
        Rcpp::IntegerVector curass=assignments[g];
        for (auto curdex : curass) {
            if (curdex <= 0 || curdex > ncells) {
                throw std::runtime_error("cell assignment indices out of range");
            }
        }
                
        // Computing the median intensity.
        for (size_t u=0; u<nmarkers; ++u) {
            auto curexprs=exprs.column(u);
            for (auto curdex : curass) { 
                --curdex; // 1-based indexing.
                all_intensities[sample_id[curdex] - 1].push_back(curexprs[curdex]);
            }

            auto curout=output[u].row(g);
            for (int s=0; s<nsamples; ++s) {
                std::deque<double>& curint=all_intensities[s];
                if (curint.empty()) {
                    curout[s]=R_NaReal;
                    continue;
                }

                size_t mid=curint.size()/2;
                std::nth_element(curint.begin(), curint.begin() + mid, curint.end());
                if (curint.size()%2==0) { 
                    double medtmp=curint[mid];
                    std::nth_element(curint.begin(), curint.begin()+mid-1, curint.end());
                    curout[s]=(medtmp+curint[mid-1])/2;
                } else {
                    curout[s]=curint[mid];
                }
                curint.clear();
            }
        }
    }

    // Converting to a list.
    return Rcpp::List(output.begin(), output.end());
}
