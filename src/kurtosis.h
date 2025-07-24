#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double compare(double *a, double *b);
double comb(int n, int k);
double nanmean(double *arr, int size);
double l_kurt(double *x, int n);
double moors_kurt(double *x, int n);
double moors_octile_kurt(double *x, int n);
double hogg_kurt(double *x, int n);
double crow_siddiqui_kurt(double *x, int n);
double reza_ma_kurt(double *x, int n);
