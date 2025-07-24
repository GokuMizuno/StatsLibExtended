#include "skewness.h"
#include "central_tendency.h"

double compare(double *a, double *b) {
    return (*a - *b);
}

double l_skew(double *x, size_t n) {
    double *_x = malloc(n * sizeof(double));
    for (size_t i = 0; i < n; i++) {
        _x[i] = x[i];
    }
    qsort(_x, n, sizeof(double), compare);
    
    double common[3];
    for (size_t i = 0; i < 3; i++) {
        common[i] = 1.0 / comb(n - 1, i) / n;
    }
    
    double betas[3] = {0};
    for (size_t i = 0; i < 3; i++) {
        for (size_t j = i; j < n; j++) {
            betas[i] += comb(j, i) * _x[j];
        }
        betas[i] *= common[i];
    }
    
    double l3 = 6 * betas[2] - 6 * betas[1] + betas[0];
    double l2 = 2 * betas[1] - betas[0];
    
    free(_x);
    return l3 / l2;
}

double pearson_mode_skew(double *x, size_t n) {
    double mean = nanmean(x, n);
    double mode = mode(x, n);
    double std = nanstd(x, n);
    return (mean - mode) / std;
}

double bickel_mode_skew(double *x, size_t n) {
    double mode = half_sample_mode(x, n);
    return nanmean(x, n) / fabs(x - mode);
}

double pearson_median_skew(double *x, size_t n) {
    double mean = nanmean(x, n);
    double median = nanmedian(x, n);
    double std = nanstd(x, n);
    return 3 * (mean - median) / std;
}

double medeen_skew(double *x, size_t n) {
    double median = nanmedian(x, n);
    double mean = nanmean(x, n);
    return (mean - median) / nanmean(fabs(x - median), n);
}

double bowley_skew(double *x, size_t n) {
    double q1 = nanquantile(x, n, 0.25);
    double q2 = nanquantile(x, n, 0.5);
    double q3 = nanquantile(x, n, 0.75);
    return (q3 + q1 - 2 * q2) / (q3 - q1);
}

double groeneveld_skew(double *x, size_t n) {
    double q1 = nanquantile(x, n, 0.25);
    double q2 = nanquantile(x, n, 0.5);
    double q3 = nanquantile(x, n, 0.75);
    double rs = (q3 + q1 - 2 * q2) / (q2 - q1);
    double ls = (q3 + q1 - 2 * q2) / (q3 - q2);
    return fabs(rs) > fabs(ls) ? rs : ls;
}

double kelly_skew(double *x, size_t n) {
    double d1 = nanquantile(x, n, 0.1);
    double d5 = nanquantile(x, n, 0.5);
    double d9 = nanquantile(x, n, 0.9);
    return (d9 + d1 - 2 * d5) / (d9 - d1);
}

double hossain_adnan_skew(double *x, size_t n) {
    double median = nanmedian(x, n);
    double *diff = malloc(n * sizeof(double));
    for (size_t i = 0; i < n; i++) {
        diff[i] = x[i] - median;
    }
    return nanmean(diff, n) / nanmean(fabs(diff), n);
}

double forhad_shorna_rank_skew(double *x, size_t n) {
    double mr = (nanmin(x, n) + nanmax(x, n)) * 0.5;
    double *arr = malloc((n + 1) * sizeof(double));
    for (size_t i = 0; i < n; i++) {
        arr[i] = x[i];
    }
    arr[n] = mr;
    
    double *arr_ranked = rankdata(arr, n + 1);
    double *diff = malloc(n * sizeof(double));
    for (size_t i = 0; i < n; i++) {
        diff[i] = arr_ranked[n] - arr_ranked[i];
    }
    
    double sum_diff = nansum(diff, n);
    double sum_abs_diff = nansum(fabs(diff), n);
    
    free(arr);
    free(arr_ranked);
    free(diff);
    
    return sum_diff / sum_abs_diff;
}

double _auc_skew_gamma(double *x, size_t n, double dp, double *w) {
    size_t half_n = n / 2;
    double *qs = malloc((n + 1) * sizeof(double));
    for (size_t i = 0; i < n; i++) {
        qs[i] = nanquantile(x, n, (double)i / n);
    }
    double med = qs[n - 1];
    double *qs_low = malloc(half_n * sizeof(double));
    double *qs_high = malloc(half_n * sizeof(double));
    
    for (size_t i = 0; i < half_n; i++) {
        qs_low[i] = qs[i];
        qs_high[i] = qs[n - 1 - i];
    }
    
    double *skews = malloc(half_n * sizeof(double));
    for (size_t i = 0; i < half_n; i++) {
        skews[i] = (qs_low[i] + qs_high[i] - 2 * med) / (qs_high[i] - qs_low[i]) * w[i];
    }
    
    double result = trapezoid(skews, half_n, dp);
    
    free(qs);
    free(qs_low);
    free(qs_high);
    free(skews);
    
    return result;
}

double auc_skew_gamma(double *x, size_t n, double dp) {
    double w = 1.0;
    return _auc_skew_gamma(x, n, dp, &w);
}

double wauc_skew_gamma(double *x, size_t n, double dp) {
    size_t half_n = n / 2;
    double *w = malloc(half_n * sizeof(double));
    for (size_t i = 0; i < half_n; i++) {
        w[i] = (double)(half_n - i) / half_n;
    }
    return _auc_skew_gamma(x, n, dp, w);
}

double cumulative_skew(double *x, size_t n) {
    double *p = malloc(n * sizeof(double));
    for (size_t i = 0; i < n; i++) {
        p[i] = x[i];
    }
    nancumsum(p, n);
    double sum_p = p[n - 1];
    for (size_t i = 0; i < n; i++) {
        p[i] /= sum_p;
    }
    double *r = malloc(n * sizeof(double));
    for (size_t i = 0; i < n; i++) {
        r[i] = (double)i / n;
    }
    double *d = malloc(n * sizeof(double));
    for (size_t i = 0; i < n; i++) {
        d[i] = r[i] - p[i];
    }
    double *w = malloc(n * sizeof(double));
    for (size_t i = 0; i < n; i++) {
        w[i] = (2 * i - n) * 3.0 / n;
    }
    double result = nansum(d, n) / nansum(d, n);
    
    free(p);
    free(r);
    free(d);
    free(w);
    
    return result;
}
