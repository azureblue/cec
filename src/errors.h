#ifndef ERRORS_H
#define	ERRORS_H

enum cec_result_code
{
    NO_ERROR = 0,
            
    MEM_ALLOC_ERROR = 1,
    INVALID_COVARIANCE_ERROR = 2,
    ALL_CLUSTERS_REMOVED_ERROR = 3,
    CENTERS_INIT_ERROR = 4,
    INVALID_CENTERS_INIT_METHOD_ERROR = 5,
    LIBRARY_DEFECT_ERROR = 6,

    UNKNOWN_ERROR = -1
};

typedef enum cec_result_code res_code;

#endif	/* ERRORS_H */
