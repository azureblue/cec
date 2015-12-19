#ifndef ERRORS_H
#define	ERRORS_H
#include "errors.h"

extern const char * const INVALID_COVARIANCE_ERROR_MSG;

extern const char * const MALLOC_ERROR_MSG;

extern const char * const UNKNOWN_ERROR_MSG;

extern const char * const ALL_CLUSTERS_REMOVED_MSG;

enum error_code
{
    NO_ERROR = 0,
    MALLOC_ERROR = 1,
    INVALID_COVARIANCE_ERROR = 2,
    ALL_CLUSTERS_REMOVED_ERROR = 3,
    UNKNOWN_ERROR = 4
};

#endif	/* ERRORS_H */
