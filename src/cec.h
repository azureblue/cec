#ifndef CEC_H
#define	CEC_H

#include "cec_context.h"

#define BIG_DOUBLE 100000000000000.0

static inline double compute_energy(const int m, const double hx,
	const int card)
{
    double p = (card / (double) m);
    return p * (-log(p) + hx);
}

int cec(struct cec_context * cec_context);

#endif	/* CEC_H */

