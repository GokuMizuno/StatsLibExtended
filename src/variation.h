#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int unique_counts(const double *x, int n, double **unique_vals);
double mod_vr(const double *x, int n);
double range_vr(const double *x, int n);
double b_index(const double *x, int n);
double avdev(const double *x, int n);
double renyi_entropy(const double *x, int n, double alpha);
double negative_extropy(const double *x, int n);
double mcintosh_d(const double *x, int n);
double gibbs_m1(const double *x, int n);
double gibbs_m2(const double *x, int n);
