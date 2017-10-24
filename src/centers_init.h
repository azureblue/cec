#ifndef CENTERS_INIT_H
#define CENTERS_INIT_H

#include "matrix.h"
#include "cec_params.h"

typedef memptr_t centers_init_ctx;
typedef void (*centers_init_function)(centers_init_ctx ctx, const cec_mat * x, cec_mat * c);

typedef struct {
    centers_init_function init;
    centers_init_ctx ctx;
} centers_initializer;

void cec_init_centers(centers_initializer * ci, const cec_mat * x, cec_mat * c);
centers_initializer * create_centers_init(cec_centers_par *centers_par, int m_max);
centers_initializer *create_copy_initializer(cec_mat *centers);

#endif //CENTERS_INIT_H
