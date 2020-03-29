#include "Rcpp.h"
#include <algorithm>
#include <vector>

// [[Rcpp::export(rng=false)]]
SEXP weighted_median_int(
    Rcpp::NumericMatrix exprs,
    Rcpp::List assignments,
    Rcpp::NumericVector sample_id, 
    Rcpp::NumericVector sample_weight)
{
    const size_t nmarkers=exprs.ncol();
    const size_t ncells=exprs.nrow();
    const int ngroups=assignments.size();

    if (sample_id.size()!=ncells) { 
        throw std::runtime_error("sample IDs should be an integer vector of length equal to the number of cells"); 
    }

    // Setting up output vectors. 
    Rcpp::NumericMatrix outcoords(nmarkers, ngroups);
    std::vector<std::pair<double, int> > intensities;
    intensities.reserve(ncells);

    for (int g=0; g<ngroups; ++g) {
        Rcpp::IntegerVector curass=assignments[g];
            
        // Computing total weights.
        double total_weight=0;
        for (auto c : curass) { 
            int cursample=sample_id[c-1];
            total_weight+=sample_weight[cursample];
        }

        // Setting the weighted medians (to avoid large samples from dominating the location).
        auto curcoords=outcoords.column(g);
        intensities.resize(curass.size());

        for (size_t mi=0; mi<nmarkers; ++mi) {
            auto curexprs=exprs.column(mi);

            for (size_t icx=0; icx<curass.size(); ++icx) {
                const int curneighbor=curass[icx]-1;
                intensities[icx].first=curexprs[curneighbor];
                intensities[icx].second=sample_id[curneighbor];
            }

            std::sort(intensities.begin(), intensities.end());
            double cumweight=0;
            size_t midpoint=0;
            bool exactly_mid=false;

            for (; midpoint<intensities.size(); ++midpoint) {
                cumweight += sample_weight[intensities[midpoint].second];
                double ratio=cumweight/total_weight;
                if (ratio >= 0.49999999) {
                    if (ratio <= 0.50000001) {
                        exactly_mid=true;
                    }
                    break;
                }
            }
   
            if (midpoint==intensities.size()) {
                // Only possible if total_weights is zero. 
                curcoords[mi]=R_NaReal;
            } else {
                if (exactly_mid) {
                    curcoords[mi]=(intensities[midpoint].first + intensities[midpoint+1].first)/2;
                } else {
                    curcoords[mi]=intensities[midpoint].first;                
                }
            }
        }
    }

    return outcoords;
}

