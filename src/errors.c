#include "errors.h"

const char * const INVALID_COVARIANCE_ERROR_MSG =
	"There was a group with invalid covariance matrix (perhaps not positive-definite).";

const char * const MALLOC_ERROR_MSG = "Memory allocation error.";

const char * const UNKNOWN_ERROR_MSG = "Internal library error.";

const char * const ALL_CLUSTERS_REMOVED_MSG =
	"All clusters have been removed before first iteraions.";
