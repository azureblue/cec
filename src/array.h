#ifndef ARRAY_H
#define	ARRAY_H

/*
 * Simple operations on C arrays.
 */

void array_fill(double * x, const double val, const int n);

void array_copy(const double * from, double * to, const int n);

void array_add(double * x1, const double *x2, const int n);

void array_sub(double * x1, const double * x2, const int n);

void array_sum_multiplied(double * dest, const double * x1, const double m1,
	const double *x2, const double m2, const int n);

void array_mul(double * x1, const double val, const int n);

double dist2(const double * x1, const double * x2, const int n);

void mean_add_point(const double * mean, double * new_mean,
	const double * point, const int card, const int n);

void mean_remove_point(const double * mean, double * new_mean,
	const double * point, const int card, const int n);

#endif	/* ARRAY_H */

