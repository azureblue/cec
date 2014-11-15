#ifndef KMEANSPP_H
#define	KMEANSPP_H

#include "matrix.h"

/*
 * K-means++ algorithm.
 */
int kmeanspp(struct cec_matrix * X, struct cec_matrix * C);

#endif	/* KMEANSPP_H */
