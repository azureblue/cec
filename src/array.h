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

vec_i * cec_array_int_create_from(int len, int * src);

vec_d * cec_array_double_create_from(int len, double * src);

vec_d * cec_array_double_create(int len);

void cec_array_int_copy(vec_i *src, vec_i *dst);
void cec_array_double_copy(vec_d *src, vec_d *dst);
void cec_array_int_copy_to(vec_i *src, int *dst);
void cec_array_double_copy_to(vec_d *src, double *dst);

#define vec_i_create cec_array_int_create
#define vec_d_create cec_array_double_create
#define vec_i_create_from cec_array_int_create_from
#define vec_d_create_from cec_array_double_create_from
#define vec_i_copy cec_array_int_copy
#define vec_d_copy cec_array_double_copy
#define vec_i_copy_to cec_array_int_copy_to
#define vec_d_copy_to cec_array_double_copy_to

void array_fill(double * x, const double val, const int n);

void array_copy(const double *restrict from, double *restrict to, const int n);

void array_add(double *restrict x1, const double *restrict x2, const int n);

void array_sub(double *restrict x1, const double *restrict x2, const int n);

void array_mul(double * x, const double val, const int n);

double dist_sq(const double *x, const double *x2, const int n);

void array_sum_mul(double *restrict dest, const double *restrict x1, const double m1,
        const double *restrict x2, const double m2, const int n);

void mean_add_point(double *restrict new_mean, const double *restrict mean, 
        const double *restrict point, const int card, const int n);

void mean_remove_point(double *restrict new_mean, const double *restrict mean,
        const double *restrict point, const int card, const int n);


#endif	/* ARRAY_H */
