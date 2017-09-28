#include "cydar.h"
#include "packer.h"
#include "utils.h"
#include "objects.h"

SEXP compute_hyperstats(SEXP exprs, SEXP nsamp, SEXP sample_id, SEXP assignments) {
    BEGIN_RCPP 

    // Setting up inputs.
    const Rcpp::NumericMatrix _exprs(exprs);
    const size_t nmarkers=_exprs.nrow();
    const size_t ncells=_exprs.ncol();

    const Rcpp::List _assignments(assignments);
    const int ngroups=_assignments.size();

    // Checking samples and computing sample weights.
    const int nsamples=check_integer_scalar(nsamp, "number of samples");
    if (nsamples <= 0) { throw std::runtime_error("number of samples must be positive"); }

    const Rcpp::IntegerVector _sample_id(sample_id);
    if (_sample_id.size()!=ncells) { 
        throw std::runtime_error("sample IDs should be an integer vector of length equal to the number of cells"); 
    }

    std::vector<double> sample_weights(nsamples);
    for (const auto& s : _sample_id) {
        if (s < 0 || s >= nsamples) {
            throw std::runtime_error("sample IDs out of range");
        }
        ++(sample_weights[s]);
    }
    for (auto& w : sample_weights) { // Reciprocal of the total number of cells.
        w=1/w;
    }

    // Setting up output vectors. 
    Rcpp::IntegerMatrix outcounts(ngroups, nsamples);
    Rcpp::NumericMatrix outcoords(ngroups, nmarkers);

    std::deque<std::pair<double, int> > intensities;
    std::deque<int> collected;

    for (int g=0; g<ngroups; ++g) {
        const Rcpp::IntegerVector curass=_assignments[g];
        unpack_index_vector(collected, curass.begin(), curass.end());
        for (size_t icx=0; icx<collected.size(); ++icx) { --(collected[icx]); } // Getting to 1-based indexing.
            
        // Computing counts and total weights.
        auto curcounts=outcounts.row(g);
        double total_weight=0;
        for (const auto& c : collected) { 
            const int& cursample=_sample_id[c];
            ++(curcounts[cursample]);
            total_weight+=sample_weights[cursample];
        }

        // Setting the weighted medians (to avoid large samples from dominating the location).
        intensities.resize(collected.size());
        auto curcoords=outcoords.row(g);
        for (size_t mi=0; mi<nmarkers; ++mi) {
            auto curexprs=_exprs.row(mi);
            for (size_t icx=0; icx<collected.size(); ++icx) {
                const int& curneighbor=collected[icx];
                intensities[icx].first=curexprs[curneighbor];
                intensities[icx].second=_sample_id[curneighbor];
            }

            std::sort(intensities.begin(), intensities.end());
            double cumweight=0;
            size_t midpoint=0;
            for (; midpoint<intensities.size(); ++midpoint) {
                cumweight += sample_weights[intensities[midpoint].second];
                if (cumweight/total_weight >= 0.5) { break; }
            }
   
            if (midpoint==intensities.size()) {
                // Only possible if total_weights is zero. 
                curcoords[mi]=R_NaReal;
            } else {
                if (cumweight/total_weight==0.5) {
                    curcoords[mi]=(intensities[midpoint].first + intensities[midpoint+1].first)/2;
                } else {
                    curcoords[mi]=intensities[midpoint].first;                
                }
            }
        }
    }

    return Rcpp::List::create(outcounts, outcoords);
    END_RCPP
}

