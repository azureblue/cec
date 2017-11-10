#ifndef ARRAY_H
#define	ARRAY_H

struct cec_array_int {
    const int len;
    int ar[];
};

struct cec_array_double {
    const int len;
    double ar[];
};

typedef struct cec_array_int vec_i;
typedef struct cec_array_double vec_d;

vec_i * cec_array_int_create(int len);
vec_i * cec_array_int_create_from(int len, const int * src);
vec_d * cec_array_double_create_from(int len,const double * src);
vec_d * cec_array_double_create(int len);

const int * cec_array_int_const_data(const vec_i *src);
const double * cec_array_double_const_data(const vec_d *src);
int * cec_array_int_data(vec_i *src);
double * cec_array_double_data(vec_d *src);

void cec_array_int_copy(const vec_i *src, vec_i *dst);
void cec_array_double_copy(const vec_d *src, vec_d *dst);
void cec_array_int_copy_to(const vec_i *src, int *dst);
void cec_array_double_copy_to(const vec_d *src, double *dst);

#define vec_i_create cec_array_int_create
#define vec_d_create cec_array_double_create
#define vec_i_create_from cec_array_int_create_from
#define vec_d_create_from cec_array_double_create_from
#define vec_i_copy cec_array_int_copy
#define vec_d_copy cec_array_double_copy
#define vec_i_copy_to cec_array_int_copy_to
#define vec_d_copy_to cec_array_double_copy_to
#define vec_i_cdata cec_array_int_const_data
#define vec_i_data cec_array_int_data
#define vec_d_cdata cec_array_double_const_data
#define vec_d_data cec_array_double_data

int max_i(const vec_i *ar);

void array_fill(double * x, double val, int n);

void array_copy(const double *restrict from, double *restrict to, int n);

void array_add(double *restrict x1, const double *restrict x2, int n);

void array_sub(double *restrict x1, const double *restrict x2, int n);

void array_mul(double * x, double val, int n);

double dist_sq(const double *x, const double *x2, int n);

void array_sum_mul(double *restrict dest, const double *restrict x1, double m1,
        const double *restrict x2, double m2, int n);

void mean_add_point(double *restrict new_mean, const double *restrict mean, 
        const double *restrict point, int card, int n);

void mean_remove_point(double *restrict new_mean, const double *restrict mean,
        const double *restrict point, int card, int n);

#endif	/* ARRAY_H */
