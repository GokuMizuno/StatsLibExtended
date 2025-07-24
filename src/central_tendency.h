#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

double mean(double* x, int n);
double median(double* x, int n);
int compare(const void* a, const void* b);
double* mode(double* x, int n, int* mode_count);
double midrange(double* x, int n);
double midhinge(double* x, int n);
double trimean(double* x, int n);
double contraharmonic_mean(double* x, int n);
double midmean(double* x, int n);
double hodges_lehmann_sen_location(double* x, int n);
double standard_trimmed_harrell_davis_quantile(double* x, int n, double q);
double half_sample_mode(double* x, int n);
