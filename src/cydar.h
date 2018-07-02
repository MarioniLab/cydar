#ifndef CYDAR_H
#define CYDAR_H

#include <stdexcept>
#include <algorithm>
#include <map>
#include <deque>
#include <cmath>
#include <vector>
#include <queue>
#include <memory>

#include "Rcpp.h"

extern "C" {

SEXP compute_density(SEXP, SEXP);

SEXP compute_hyperstats(SEXP, SEXP, SEXP, SEXP);

SEXP drop_redundant(SEXP, SEXP, SEXP);

SEXP compute_median_int(SEXP, SEXP, SEXP, SEXP);

}

#endif 
