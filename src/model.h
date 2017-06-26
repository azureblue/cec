#ifndef MODEL_H
#define MODEL_H

#include <Rdefines.h>
#include "models/energy.h"

struct cec_model
{
    cross_entropy_function cross_entropy;    
    struct cross_entropy_context * cross_entropy_context; 
};

double cluster_energy(struct cec_model * model, struct cec_matrix * cov, int card, int m);

res_code cluster_energy_get_last_error(struct cec_model * model);

#endif /* MODEL_H */

