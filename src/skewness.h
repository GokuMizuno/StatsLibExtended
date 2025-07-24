#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

double compare(double *a, double *b);
double l_skew(double *x, size_t n);
double pearson_mode_skew(double *x, size_t n);
double bickel_mode_skew(double *x, size_t n);
double pearson_median_skew(double *x, size_t n);
double medeen_skew(double *x, size_t n);
double bowley_skew(double *x, size_t n);
double groeneveld_skew(double *x, size_t n);
double kelly_skew(double *x, size_t n);
double hossain_adnan_skew(double *x, size_t n);
double forhad_shorna_rank_skew(double *x, size_t n);
double _auc_skew_gamma(double *x, size_t n, double dp, double *w);
double auc_skew_gamma(double *x, size_t n, double dp);
double wauc_skew_gamma(double *x, size_t n, double dp);
double cumulative_skew(double *x, size_t n);
