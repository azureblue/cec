#ifndef ERRORS_H
#define	ERRORS_H

enum error_code
{
    NO_ERROR = 0,
            
    MALLOC_ERROR = 1,
    INVALID_COVARIANCE_ERROR = 2,
    ALL_CLUSTERS_REMOVED_ERROR = 3,
    CENTERS_INIT_ERROR = 4,

    UNKNOWN_ERROR = -1
};

#endif	/* ERRORS_H */
