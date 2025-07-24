#include "dispersion.h"

double compare(double *a, double *b) {
    return (*a - *b);
}

double nanmax(const double *x, size_t n) {
    double max_val = -DBL_MAX;
    for (size_t i = 0; i < n; i++) {
        if (!isnan(x[i]) && x[i] > max_val) {
            max_val = x[i];
        }
    }
    return max_val;
}

double nanmin(const double *x, size_t n) {
    double min_val = DBL_MAX;
    for (size_t i = 0; i < n; i++) {
        if (!isnan(x[i]) && x[i] < min_val) {
            min_val = x[i];
        }
    }
    return min_val;
}

double nanmean(const double *x, size_t n) {
    double sum = 0.0;
    int count = 0;
    for (size_t i = 0; i < n; i++) {
        if (!isnan(x[i])) {
            sum += x[i];
            count++;
        }
    }
    return count > 0 ? sum / count : NAN;
}

double nanstd(const double *x, size_t n) {
    double mean = nanmean(x, n);
    double sum_sq_diff = 0.0;
    int count = 0;
    for (size_t i = 0; i < n; i++) {
        if (!isnan(x[i])) {
            sum_sq_diff += (x[i] - mean) * (x[i] - mean);
            count++;
        }
    }
    return count > 1 ? sqrt(sum_sq_diff / (count - 1)) : NAN;
}

double studentized_range(const double *x, size_t n) {
    double maximum = nanmax(x, n);
    double minimum = nanmin(x, n);
    double std = nanstd(x, n);
    return (maximum - minimum) / std;
}

double coefficient_of_lvariation(const double *x, size_t n) {
    double l1 = nanmean(x, n);
    if (fabs(l1) <= EPS) {
        fprintf(stderr, "Mean is close to 0. Statistic is undefined.\n");
        return INFINITY;
    }
    double *sorted_x = malloc(n * sizeof(double));
    memcpy(sorted_x, x, n * sizeof(double));
    qsort(sorted_x, n, sizeof(double), compare);
    double common = 1.0 / (n - 1) / n;
    double beta_1 = common * nansum(sorted_x + 1, n - 1);
    double l2 = 2 * beta_1 - l1;
    free(sorted_x);
    return l2 / l1;
}

double coefficient_of_variation(const double *x, size_t n) {
    double mean = nanmean(x, n);
    if (fabs(mean) <= EPS) {
        fprintf(stderr, "Mean is close to 0. Statistic is undefined.\n");
        return INFINITY;
    }
    return nanstd(x, n) / mean;
}

double robust_coefficient_of_variation(const double *x, size_t n) {
    double med = nanmedian(x, n);
    if (fabs(med) <= EPS) {
        fprintf(stderr, "Median is close to 0. Statistic is undefined.\n");
        return INFINITY;
    }
    double med_abs_dev = nanmedian(abs_diff(x, med, n), n);
    return med_abs_dev / med;
}

double quartile_coefficient_of_dispersion(const double *x, size_t n) {
    double q1 = nanquantile(x, n, 0.25);
    double q3 = nanquantile(x, n, 0.75);
    if (fabs(q3 + q1) <= EPS) {
        fprintf(stderr, "Midhinge is close to 0. Statistic is undefined.\n");
        return INFINITY;
    }
    return (q3 - q1) / (q3 + q1);
}

// double dispersion_ratio(const double *x, size_t n) {
//     double *non_zero_x = malloc(n * sizeof(double));
//     for (size_t i = 0; i < n; i++) {
//         non_zero_x[i] = (x[i] == 0) ? NAN : x[i];
//     }
//     double mean = nanmean(x, n);
//     double gmean = gmean(non_zero_x, n);
//     free(non_zero_x);
//     return mean / gmean;
// }

double fisher_index_of_dispersion(const double *x, size_t n) {
    double mean = nanmean(x, n);
    if (fabs(mean) <= EPS) {
        fprintf(stderr, "Mean is close to 0. Statistic is undefined.\n");
        return INFINITY;
    }
    return (n - 1) * nanvar(x, n) / mean;
}

double morisita_index_of_dispersion(const double *x, size_t n) {
    double x_sum = nansum(x, n);
    return n * (nansum_square(x, n) - x_sum) / (x_sum * x_sum - x_sum);
}

double standard_quantile_absolute_deviation(const double *x, size_t n) {
    double med = nanmedian(x, n);
    double k = 1.0 + 0.762 / n + 0.967 / (n * n);
    double q = 0.6826894921370850;  // stats.norm.cdf(1) - stats.norm.cdf(-1)
    return k * nanquantile(abs_diff(x, med, n), n, q);
}

double shamos_estimator(const double *x, size_t n) {
    double *product = malloc(n * n * sizeof(double));
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            product[i * n + j] = fabs(x[i] - x[j]);
        }
    }
    double result = nanmedian(product, n * n);
    free(product);
    return result;
}

double coefficient_of_range(const double *x, size_t n) {
    double min_ = nanmin(x, n);
    double max_ = nanmax(x, n);
    if (fabs(min_ + max_) <= EPS) {
        fprintf(stderr, "Midrange is close to 0. Statistic is undefined.\n");
        return INFINITY;
    }
    return (max_ - min_) / (max_ + min_);
}

double cole_index_of_dispersion(const double *x, size_t n) {
    return nansum_square(x, n) / pow(nansum(x, n), 2);
}

double gini_mean_difference(const double *x, size_t n) {
    double *product = malloc(n * n * sizeof(double));
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            product[i * n + j] = fabs(x[i] - x[j]);
        }
    }
    double result = nansum(product, n * n) / (n * (n - 1));
    free(product);
    return result;
}
