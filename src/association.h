#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>

bool check_arrays(double *x, double *y, int len);
void prep_arrays(double *x, double *y, double **_x, double **_y, int *new_len, int len);
double chatterjeexi(double *x, double *y, int len);
double concordance_corrcoef(double *x, double *y, int len);
double concordance_rate(double *x, double *y, size_t n);
double symmetric_chatterjeexi(double *x, double *y, size_t n);
double zhangi(double *x, double *y, size_t n);
double tanimoto_similarity(double *x, double *y, size_t n);
double blomqvistbeta(double *x, double *y, size_t n);
double winsorized_correlation(double *x, double *y, size_t n, double k);
double rank_minrelation_coefficient(double *x, double *y, size_t n);
double tukey_correlation(double *x, double *y, size_t n);
double gaussian_rank_correlation(double *x, double *y, int n);