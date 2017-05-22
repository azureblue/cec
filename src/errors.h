#ifndef ERRORS_H
#define	ERRORS_H

enum cec_result_code
{
    NO_ERROR = 0,
            
    MEM_ALLOC_ERROR = 1,
    INVALID_COVARIANCE_ERROR = 2,
    ALL_CLUSTERS_REMOVED_ERROR = 3,
    CENTERS_INIT_ERROR = 4,
    INVALID_CENTERS_INIT_METHOD = 5,

    UNKNOWN_ERROR = -1
};

typedef enum cec_result_code cec_res;

#endif	/* ERRORS_H */
