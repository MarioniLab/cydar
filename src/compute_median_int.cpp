#include "cydar.h"
#include "packer.h"
#include "utils.h"
#include "objects.h"

SEXP compute_median_int(SEXP exprs, SEXP nsamp, SEXP sample_id, SEXP assignments) {
    BEGIN_RCPP

    // Setting up inputs.
    const Rcpp::NumericMatrix _exprs(exprs);
    const size_t& nmarkers=_exprs.nrow();
    const size_t& ncells=_exprs.ncol();

    const Rcpp::List _assignments(assignments);
    const int ngroups=_assignments.size();

    // Checking samples.
    const int nsamples=check_integer_scalar(nsamp, "number of samples");
    if (nsamples <= 0) { throw std::runtime_error("number of samples must be positive"); }

    Rcpp::IntegerVector _sample_id(sample_id);
    if (_sample_id.size()!=ncells) { 
        throw std::runtime_error("sample IDs should be an integer vector of length equal to the number of cells"); 
    }
    for (const auto& cursample : _sample_id) { 
        if (cursample < 0 || cursample >= nsamples) {
            throw std::runtime_error("sample IDs out of range");
        }
    }

    // Setting up output matrices (don't use List here, as it gets stored as a funny Proxy type).
    std::vector<Rcpp::NumericMatrix> output(nmarkers);
    for (size_t u=0; u<nmarkers; ++u) {
        output[u]=Rcpp::NumericMatrix(ngroups, nsamples);
    }

    std::deque<int> collected;
    std::vector<std::deque<double> > all_intensities(nsamples);

    for (int g=0; g<ngroups; ++g) {
        // Unpacking the assignments.
        const Rcpp::IntegerVector curass(_assignments[g]);
        unpack_index_vector(collected, curass.begin(), curass.end());
        for (auto& curdex : collected) {
            --curdex; // Getting to zero-index.
            if (curdex < 0 || curdex >= ncells) {
                throw std::runtime_error("cell assignment indices out of range");
            }
        }
                
        // Computing the median intensity.
        for (size_t u=0; u<nmarkers; ++u) {
            auto curexprs=_exprs.row(u);
            for (const auto& curdex : collected) { 
                all_intensities[_sample_id[curdex]].push_back(curexprs[curdex]);
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
    END_RCPP
}

