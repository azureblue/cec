#ifndef CEC_H
#define	CEC_H

#include "cec_context.h"

#define BIG_DOUBLE 100000000000000.0

/*
 * Performs the CEC algorithm on the cec_context structure.
 */
int cec(struct cec_context * cec_context);

#endif	/* CEC_H */

