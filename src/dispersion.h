#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

#define EPS 1e-6

double compare(double *a, double *b);
double nanmax(const double *x, size_t n);
double nanmin(const double *x, size_t n);
double nanmean(const double *x, size_t n);
double nanstd(const double *x, size_t n);
double studentized_range(const double *x, size_t n);
double coefficient_of_lvariation(const double *x, size_t n);
double coefficient_of_variation(const double *x, size_t n);
double robust_coefficient_of_variation(const double *x, size_t n);
double quartile_coefficient_of_dispersion(const double *x, size_t n);
// double dispersion_ratio(const double *x, size_t n);
double fisher_index_of_dispersion(const double *x, size_t n);
double morisita_index_of_dispersion(const double *x, size_t n);
double standard_quantile_absolute_deviation(const double *x, size_t n);
double shamos_estimator(const double *x, size_t n);
double coefficient_of_range(const double *x, size_t n);
double cole_index_of_dispersion(const double *x, size_t n);
double gini_mean_difference(const double *x, size_t n);
