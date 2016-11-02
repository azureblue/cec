#include "array.h"

void array_fill(double * x, const double val, const int n)
{
    for (int i = 0; i < n; i++)
        x[i] = val;
}

void array_copy(const double *restrict from, double *restrict to, const int n)
{
    for (int i = 0; i < n; i++)
        to[i] = from[i];
}

void array_add(double *restrict x1, const double *restrict x2, const int n)
{
    for (int i = 0; i < n; i++)
        x1[i] += x2[i];
}

void array_sub(double *restrict x1, const double *restrict x2, const int n)
{
    for (int i = 0; i < n; i++)
        x1[i] -= x2[i];
}

void array_sum_multiplied(double *restrict dest, const double *restrict x1, const double m1,
        const double *restrict x2, const double m2, const int n)
{
    for (int i = 0; i < n; i++)
        dest[i] = x1[i] * m1 + x2[i] * m2;
}

void array_mul(double * x, const double val, const int n)
{
    for (int i = 0; i < n; i++)
        x[i] *= val;
}

double dist2(const double * x1, const double * x2, const int n)
{
    double dist = 0;
    for (int i = 0; i < n; i++)
        dist += (x1[i] - x2[i]) * (x1[i] - x2[i]);
    return dist;
}

void mean_add_point(double *restrict dest_mean, const double *restrict mean,
        const double *restrict point, const int card, const int n)
{
    double d = (double) card;
    array_sum_multiplied(dest_mean, mean, d / (d + 1.0), point, 1.0 / (d + 1.0), n);
}

void mean_remove_point(double *restrict dest_mean, const double *restrict mean,
        const double *restrict point, const int card, const int n)
{
    double d = (double) card;
    array_sum_multiplied(dest_mean, mean, d / (d - 1.0), point, -1.0 / (d - 1.0), n);
}
