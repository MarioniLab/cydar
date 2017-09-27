#ifndef CYDAR_H
#define CYDAR_H

#include <stdexcept>
#include <algorithm>
#include <map>
#include <deque>
#include <bitset>
#include <cmath>

#include "Rcpp.h"

extern "C" {

SEXP compute_density(SEXP, SEXP, SEXP, SEXP);

SEXP compute_hyperstats(SEXP, SEXP, SEXP, SEXP);

SEXP count_cells(SEXP, SEXP, SEXP, SEXP, SEXP);

SEXP get_nndist(SEXP, SEXP, SEXP, SEXP, SEXP);

SEXP find_knn(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP); 

SEXP drop_redundant(SEXP, SEXP, SEXP, SEXP, SEXP);

SEXP pack_indices(SEXP, SEXP);

SEXP recount_cells(SEXP, SEXP, SEXP, SEXP);

SEXP compute_median_int(SEXP, SEXP, SEXP, SEXP);

}

#endif 
