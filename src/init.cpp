#include "cydar.h"
#include "R_ext/Rdynload.h"
#include "R_ext/Visibility.h"

#define REGISTER(x, i) {#x, (DL_FUNC) &x, i}

extern "C" {

static const R_CallMethodDef all_call_entries[] = {
    REGISTER(compute_density, 2),
    REGISTER(compute_hyperstats, 4),
    REGISTER(drop_redundant, 2),
    REGISTER(compute_median_int, 4),
    {NULL, NULL, 0}
};

void attribute_visible R_init_cydar(DllInfo *dll) {
    R_registerRoutines(dll, NULL, all_call_entries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
    R_forceSymbols(dll, TRUE);
}

}

