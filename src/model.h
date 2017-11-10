#ifndef MODEL_H
#define MODEL_H

#include "models/energy.h"

typedef struct
{
    cross_entropy_function cross_entropy;    
    struct cross_entropy_context * cross_entropy_context; 
} cec_model;

double cluster_energy(const cec_model * model, cec_mat * cov, int card, int m);

res_code cluster_energy_get_last_error(const cec_model * model);

#endif /* MODEL_H */
