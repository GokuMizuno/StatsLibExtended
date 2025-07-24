#include "kurtosis.h"

double compare(double *a, double *b) {
    return (*a - *b);
}

double comb(int n, int k) {
    if (k > n) return 0;
    if (k == 0 || k == n) return 1;
    double res = 1;
    for (int i = 1; i <= k; i++) {
        res *= (double)(n - i + 1) / i;
    }
    return res;
}

double nanmean(double *arr, int size) {
    double sum = 0;
    int count = 0;
    for (int i = 0; i < size; i++) {
        if (!isnan(arr[i])) {
            sum += arr[i];
            count++;
        }
    }
    return count > 0 ? sum / count : NAN;
}

double l_kurt(double *x, int n) {
    double *_x = malloc(n * sizeof(double));
    for (int i = 0; i < n; i++) _x[i] = x[i];
    qsort(_x, n, sizeof(double), (int (*)(const void *, const void *))compare);
    
    double common[4];
    for (int i = 0; i < 4; i++) {
        common[i] = 1 / (comb(n - 1, i) * n);
    }
    
    double betas[4] = {0};
    for (int i = 0; i < 4; i++) {
        for (int j = i; j < n; j++) {
            betas[i] += comb(j, i) * _x[j];
        }
        betas[i] *= common[i];
    }
    
    double l4 = 20 * betas[3] - 30 * betas[2] + 12 * betas[1] - betas[0];
    double l2 = 2 * betas[1] - betas[0];
    free(_x);
    return l4 / l2;
}

double moors_kurt(double *x, int n) {
    double *z_scores = malloc(n * sizeof(double));
    double mean = nanmean(x, n);
    double variance = 0;
    for (int i = 0; i < n; i++) {
        if (!isnan(x[i])) {
            z_scores[i] = (x[i] - mean) / sqrt(nanmean(x, n));
            variance += z_scores[i] * z_scores[i];
        } else {
            z_scores[i] = NAN;
        }
    }
    free(z_scores);
    return variance / n + 1;
}

double moors_octile_kurt(double *x, int n) {
    double o1 = quantile(x, n, 0.125);
    double o2 = quantile(x, n, 0.25);
    double o3 = quantile(x, n, 0.375);
    double o5 = quantile(x, n, 0.625);
    double o6 = quantile(x, n, 0.75);
    double o7 = quantile(x, n, 0.875);
    return ((o7 - o5) + (o3 - o1)) / (o6 - o2);
}

double hogg_kurt(double *x, int n) {
    double p05 = quantile(x, n, 0.05);
    double p50 = quantile(x, n, 0.5);
    double p95 = quantile(x, n, 0.95);
    
    double masked_p95[n], masked_p05[n], masked_p50g[n], masked_p50l[n];
    for (int i = 0; i < n; i++) {
        masked_p95[i] = (x[i] >= p95) ? x[i] : NAN;
        masked_p05[i] = (x[i] <= p05) ? x[i] : NAN;
        masked_p50g[i] = (x[i] >= p50) ? x[i] : NAN;
        masked_p50l[i] = (x[i] <= p50) ? x[i] : NAN;
    }
    
    return (nanmean(masked_p95, n) - nanmean(masked_p05, n)) / 
           (nanmean(masked_p50g, n) - nanmean(masked_p50l, n));
}

double crow_siddiqui_kurt(double *x, int n) {
    double p025 = quantile(x, n, 0.025);
    double p25 = quantile(x, n, 0.25);
    double p75 = quantile(x, n, 0.75);
    double p975 = quantile(x, n, 0.975);
    return (p975 + p025) / (p75 - p25);
}

double reza_ma_kurt(double *x, int n) {
    double h1 = quantile(x, n, 0.0625);
    double h7 = quantile(x, n, 0.4375);
    double h9 = quantile(x, n, 0.5625);
    double h15 = quantile(x, n, 0.9375);
    return ((h15 - h9) + (h7 - h1)) / (h15 - h1);
}
