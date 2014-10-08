#include "errors.h"

const char * const POSITIVE_DEFINITE_ERROR_MSG =
	"There was a group with not positive-definite covariance matrix during computation.";

const char * const MALLOC_ERROR_MSG = "Memory allocation error.";

const char * const UNKNOWN_ERROR_MSG = "Internal library error.";

const char * const ALL_CLUSTERS_REMOVED_MSG =
	"All clusters have been removed before first iteraions.";
