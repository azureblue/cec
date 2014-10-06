#include "array.h"

void array_fill(double * x, const double val, const int n)
{
    for (int i = 0; i < n; i++)
    {
	x[i] = val;
    }
}

void array_copy(const double * from, double * to, const int n)
{
    for (int i = 0; i < n; i++)
    {
	to[i] = from[i];
    }
}

void array_add(double * x1, const double *x2, const int n)
{
    for (int i = 0; i < n; i++)
    {
	x1[i] += x2[i];
    }
}

void array_sub(double * x1, const double * x2, const int n)
{
    for (int i = 0; i < n; i++)
    {
	x1[i] -= x2[i];
    }
}

void array_sum_multiplied(double * dest, const double * x1, const double m1,
	const double *x2, const double m2, const int n)
{
    for (int i = 0; i < n; i++)
    {
	dest[i] = x1[i] * m1 + x2[i] * m2;
    }
}

void array_mul(double * x1, const double val, const int n)
{
    for (int i = 0; i < n; i++)
    {
	x1[i] *= val;
    }
}

double dist2(const double * x1, const double * x2, const int n)
{
    double res = 0;
    for (int i = 0; i < n; i++)
    {
	res += (x1[i] - x2[i]) * (x1[i] - x2[i]);
    }
    return res;
}

void mean_add_point(const double * mean, double * new_mean,
	const double * point, const int card, const int n)
{
    double d = (double) card;
    array_sum_multiplied(new_mean, mean, d / (d + 1.0), point, 1.0 / (d + 1.0),
	    n);
}

void mean_remove_point(const double * mean, double * new_mean,
	const double * point, const int card, const int n)
{
    double d = (double) card;
    array_sum_multiplied(new_mean, mean, d / (d - 1.0), point, -1.0 / (d - 1.0),
	    n);
}

