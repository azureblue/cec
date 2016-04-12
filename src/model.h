#ifndef MODEL_H
#define MODEL_H

#include "energy.h"

struct cec_model
{
    enum density_family family;
    cross_entropy_function cross_entropy;    
    struct cross_entropy_context * cross_entropy_context; 
};

double cluster_energy(struct cec_model * model, struct cec_matrix * cov, int card, int m);

int cluster_energy_get_last_error(struct cec_model * model);

#endif /* MODEL_H */

