#include "packer.h"
#include "cydar.h"
#include "objects.h"

/* For each original hypersphere, this function identifies the cells within the corresponding nested 
 * hypersphere; then identifies other nested hyperspheres which contain a subset of cells within the 
 * original hypersphere.
 */

SEXP recount_cells(SEXP exprs, SEXP distance, SEXP centers, SEXP assignments) {
    BEGIN_RCPP
    const Rcpp::NumericMatrix _exprs(exprs);
    const size_t& nmarkers=_exprs.nrow();

    // Get centers and assignments for each hypersphere.
    const Rcpp::IntegerVector _centers(centers);
    const int ncenters=_centers.size();

    const Rcpp::List _assignments(assignments);
    if (_assignments.size()!=ncenters) {
        throw std::runtime_error("length of 'assignments' must be equal to length of 'centers'");
    }
    
    const double dist=check_numeric_scalar(distance, "distance");
    const double threshold2=dist*dist;

    // Setting up output objects.
    Rcpp::List outassign(ncenters);
    Rcpp::IntegerVector outcount(ncenters);

    // Going through all centers and finding all cells that belong in the new space.
    std::deque<std::deque<int> > new_assignments(ncenters);
    std::deque<int> temp, packed;

    for (int c=0; c<ncenters; ++c) { 
        const Rcpp::IntegerVector curass(_assignments[c]);
        unpack_index_vector(temp, curass.begin(), curass.end());
            
        std::deque<int>& retained=new_assignments[c];
        const double* centerpoint=_exprs.begin() + nmarkers*(_centers[c]);
        for (size_t ix=0; ix<temp.size(); ++ix) {
            int& curt=temp[ix];
            --curt; // become zero-indexed.
                
            const double* curpoint=_exprs.begin() + nmarkers*curt;
            double dist2=0, tmpdist;
            for (size_t m=0; m<nmarkers; ++m) {
                tmpdist=centerpoint[m] - curpoint[m];
                dist2+=tmpdist*tmpdist;
            }
            if (dist2 <= threshold2) { retained.push_back(curt); }
        }

        // Storing the output.
        outcount[c]=retained.size();
        pack_index_vector(packed, retained.begin(), retained.end());
        outassign[c]=Rcpp::IntegerVector(packed.begin(), packed.end());
    }

    return Rcpp::List::create(outassign, outcount);
    END_RCPP
}
