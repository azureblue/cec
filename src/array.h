#ifndef ARRAY_H
#define	ARRAY_H

/*
 * Simple operations on C arrays.
 */

void array_fill(double * x, const double val, const int n);

void array_copy(const double *restrict from, double *restrict to, const int n);

void array_add(double *restrict x1, const double *restrict x2, const int n);

void array_sub(double *restrict x1, const double *restrict x2, const int n);

void array_mul(double * x, const double val, const int n);

double dist2(const double * x, const double * x2, const int n);

void array_sum_mul(double *restrict dest, const double *restrict x1, const double m1,
        const double *restrict x2, const double m2, const int n);

void mean_add_point(double *restrict new_mean, const double *restrict mean, 
        const double *restrict point, const int card, const int n);

void mean_remove_point(double *restrict new_mean, const double *restrict mean,
        const double *restrict point, const int card, const int n);

#endif	/* ARRAY_H */
