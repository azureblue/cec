#include <string.h>
#include <R_ext/Error.h>
#include "error_r.h"

static const char * const INVALID_COVARIANCE_ERROR_MSG =
        "There was a group with invalid covariance matrix (perhaps not positive-definite).";
static const char * const MALLOC_ERROR_MSG = "Memory allocation error.";
static const char * const UNKNOWN_ERROR_MSG = "Unknown error.";
static const char * const ALL_CLUSTERS_REMOVED_ERROR_MSG = 
        "All clusters have been removed before first iteraions.";
static const char * const CENTERS_INIT_ERROR_MSG = "Centers initialization error.";
static const char * const INVALID_CENTERS_INIT_METHOD_ERROR_MSG = "Invalid centers initialization method";
static const char * const LIBRARY_DEFECT_ERROR_MSG = "CEC library error";

static char error_msg_buffer[256];

void noreturn error_r(enum cec_result_code code) {
    switch (code)
    {
        case MEM_ALLOC_ERROR:
            error(MALLOC_ERROR_MSG);
        case INVALID_COVARIANCE_ERROR:
            error(INVALID_COVARIANCE_ERROR_MSG);
        case ALL_CLUSTERS_REMOVED_ERROR:
            error(ALL_CLUSTERS_REMOVED_ERROR_MSG);
        case CENTERS_INIT_ERROR:
            error(CENTERS_INIT_ERROR_MSG);
        case INVALID_CENTERS_INIT_METHOD_ERROR:
            error(INVALID_CENTERS_INIT_METHOD_ERROR_MSG);
        case LIBRARY_DEFECT_ERROR:
            error(LIBRARY_DEFECT_ERROR_MSG);
            
        case UNKNOWN_ERROR:
        case NO_ERROR:
            error(UNKNOWN_ERROR_MSG);
    }
}

void noreturn defect_error_r(const char *details) {
    strcpy(error_msg_buffer, LIBRARY_DEFECT_ERROR_MSG);
    strcat(error_msg_buffer, ": ");
    strcat(error_msg_buffer, details);
    error(error_msg_buffer);
}

void noreturn missing_param_error_r(const char *param_name) {
    strcpy(error_msg_buffer, LIBRARY_DEFECT_ERROR_MSG);
    strcat(error_msg_buffer, ": ");
    strcat(error_msg_buffer, "missing param: ");
    strcat(error_msg_buffer, param_name);
    error(error_msg_buffer);
}
