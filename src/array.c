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

inline void array_sum_mul(double *restrict dest, const double *restrict x1, const double m1,
        const double *restrict x2, const double m2, const int n)
{
    for (int i = 0; i < n; i++)
        dest[i] = x1[i] * m1 + x2[i] * m2;
}

void mean_add_point(double *restrict dest_mean, const double *restrict mean,
        const double *restrict point, const int card, const int n)
{
    const double card_n = card + 1;
    array_sum_mul(dest_mean, mean, card / card_n, point, 1 / card_n, n);
}

void mean_remove_point(double *restrict dest_mean, const double *restrict mean,
        const double *restrict point, const int card, const int n)
{
    const double card_n = card - 1;
    array_sum_mul(dest_mean, mean, card / card_n, point, -1 / card_n, n);
}
