#include "cydar.h"
#include "packer.h"
#include "utils.h"
#include "objects.h"

SEXP compute_median_int(SEXP exprs, SEXP nsamp, SEXP sample_id, SEXP assignments, SEXP curmarks) try {
    // Setting up inputs.
    finder fx(exprs, curmarks, R_NilValue, R_NilValue);
    const matrix_info& EXPRS=fx.searcher->exprs;
    const size_t& nmarkers=EXPRS.nrow;
    const size_t& ncells=EXPRS.ncol;
    const double* eptr=EXPRS.dptr;
    const std::deque<size_t>& used_markers=fx.searcher->get_used_markers();
    const size_t& nused=used_markers.size();

    if (!isNewList(assignments)) { 
        throw std::runtime_error("'assignments' must be a list of packed assignments");
    }
    const int ngroups=LENGTH(assignments);

    // Checking samples and computing sample weights.
    if (!isInteger(nsamp)|| LENGTH(nsamp)!=1) { throw std::runtime_error("number of samples must be an integer scalar"); }
    const int nsamples=asInteger(nsamp);
    if (nsamples <= 0) { throw std::runtime_error("number of samples must be positive"); }

    if (!isInteger(sample_id) || LENGTH(sample_id)!=ncells) {
        throw std::runtime_error("sample IDs should be an integer vector of length equal to the number of cells"); 
    }
    const int* sample_ids=INTEGER(sample_id);
    for (int i=0; i<LENGTH(sample_id); ++i) {
        const int& cursample=sample_ids[i];
        if (cursample < 0 || cursample >= nsamples) {
            throw std::runtime_error("sample IDs out of range");
        }
    }

    // Setting up output vectors. 
    SEXP output=PROTECT(allocVector(VECSXP, nused));
    try { 
        size_t u, s;
        std::deque<std::deque<double*> > outptrs(nused);
        for (u=0; u<nused; ++u) {
            SET_VECTOR_ELT(output, u, allocMatrix(REALSXP, ngroups, nsamples));
            std::deque<double*>& current=outptrs[u];
            if (nsamples) { current.push_back(REAL(VECTOR_ELT(output, u))); }
            for (s=1; s<nsamples; ++s) {
                current.push_back(current[s-1] + ngroups);
            }
        }

        SEXP curass;
        const int* iptr;
        int ndex;
        std::deque<int> collected;
        size_t icx;

        std::deque<std::deque<double> > all_intensities(nsamples);
        const double* curpoint;
        size_t mid;

        for (int g=0; g<ngroups; ++g) {
            // Unpacking the assignments.
            curass=VECTOR_ELT(assignments, g);
            if (!isInteger(curass)) { 
                throw std::runtime_error("assignment vectors should be integer");
            }
            iptr=INTEGER(curass);
            ndex=LENGTH(curass);
            unpack_index_vector(collected, iptr, iptr+ndex);
            for (icx=0; icx<collected.size(); ++icx) { 
                int& current=--(collected[icx]); // Getting to zero-index.
                if (current < 0 || current >= ncells) {
                    throw std::runtime_error("cell assignment indices out of range");
                }
            }
                
            // Computing the median intensity.
            for (u=0; u<nused; ++u) {
                curpoint = eptr + used_markers[u];

                for (icx=0; icx<collected.size(); ++icx) {
                    const int& curdex=collected[icx];
                    all_intensities[sample_ids[curdex]].push_back(*(curpoint + curdex*nmarkers));
                }

                for (s=0; s<nsamples; ++s) {
                    std::deque<double>& curint=all_intensities[s];
                    if (curint.empty()) {
                        outptrs[u][s][g]=R_NaReal;
                        continue;
                    }
                    std::sort(curint.begin(), curint.end());
                    mid=size_t(curint.size()/2);
                    if (curint.size()%2==0) { 
                        outptrs[u][s][g]=(curint[mid]+curint[mid-1])/2;
                    } else {
                        outptrs[u][s][g]=curint[mid];
                    }
                    curint.clear();
                }
            }
        }
    } catch (std::exception& e) {
        UNPROTECT(1);
        throw;
    }

    UNPROTECT(1);
    return output;
} catch (std::exception& e) {
    return mkString(e.what());
}

