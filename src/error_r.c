#include <R_ext/Error.h>
#include "error_r.h"

static const char * const INVALID_COVARIANCE_ERROR_MSG =
        "There was a group with invalid covariance matrix (perhaps not positive-definite).";

static const char * const MALLOC_ERROR_MSG = "Memory allocation error.";

static const char * const UNKNOWN_ERROR_MSG = "Unknown error.";

static const char * const ALL_CLUSTERS_REMOVED_ERROR_MSG = 
        "All clusters have been removed before first iteraions.";

static const char * const CENTERS_INIT_ERROR_MSG = "Centers initialization error.";

void noreturn error_r(enum error_code code) {
    switch (code)
    {
        case MALLOC_ERROR:
            error(MALLOC_ERROR_MSG);
        case INVALID_COVARIANCE_ERROR:
            error(INVALID_COVARIANCE_ERROR_MSG);
        case ALL_CLUSTERS_REMOVED_ERROR:
            error(ALL_CLUSTERS_REMOVED_ERROR_MSG);
        case CENTERS_INIT_ERROR:
            error(CENTERS_INIT_ERROR_MSG);
            
        default:
            error(UNKNOWN_ERROR_MSG);
    }
}