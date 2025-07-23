#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int unique_counts(const double *x, int n, double **unique_vals) {
    int *counts = (int *)calloc(n, sizeof(int));
    int unique_count = 0;

    for (int i = 0; i < n; i++) {
        int found = 0;
        for (int j = 0; j < unique_count; j++) {
            if (x[i] == (*unique_vals)[j]) {
                counts[j]++;
                found = 1;
                break;
            }
        }
        if (!found) {
            (*unique_vals)[unique_count] = x[i];
            counts[unique_count] = 1;
            unique_count++;
        }
    }

    return unique_count;
}

double mod_vr(const double *x, int n) {
    double *unique_vals = (double *)malloc(n * sizeof(double));
    int *counts = (int *)calloc(n, sizeof(int));
    int unique_count = unique_counts(x, n, &unique_vals);
    
    int max_count = 0;
    for (int i = 0; i < unique_count; i++) {
        if (counts[i] > max_count) {
            max_count = counts[i];
        }
    }

    free(unique_vals);
    free(counts);
    return 1 - (double)max_count / n;
}

double range_vr(const double *x, int n) {
    double *unique_vals = (double *)malloc(n * sizeof(double));
    int *counts = (int *)calloc(n, sizeof(int));
    int unique_count = unique_counts(x, n, &unique_vals);
    
    int min_count = counts[0];
    int max_count = counts[0];
    for (int i = 0; i < unique_count; i++) {
        if (counts[i] < min_count) {
            min_count = counts[i];
        }
        if (counts[i] > max_count) {
            max_count = counts[i];
        }
    }

    free(unique_vals);
    free(counts);
    return (double)min_count / max_count;
}

double gibbs_m1(const double *x, int n) {
    double *unique_vals = (double *)malloc(n * sizeof(double));
    int *counts = (int *)calloc(n, sizeof(int));
    int unique_count = unique_counts(x, n, &unique_vals);
    
    double sum_freq_squared = 0.0;
    for (int i = 0; i < unique_count; i++) {
        double freq = (double)counts[i] / n;
        sum_freq_squared += freq * freq;
    }
double b_index(const double *x, int n) {
    double *unique_vals = (double *)malloc(n * sizeof(double));
    int *counts = (int *)calloc(n, sizeof(int));
    int unique_count = unique_counts(x, n, &unique_vals);
    for (int i = 0; i < unique_count; i++) {
        counts[i] = (int)(counts[i] / (double)n);
    }
    double gmean = 1.0;
    for (int i = 0; i < unique_count; i++) {
        gmean *= counts[i];
    }
    gmean = pow(gmean, 1.0 / unique_count);
    free(unique_vals);
    free(counts);
    return 1 - sqrt(1 - gmean * gmean);
}

double avdev(const double *x, int n) {
    double *unique_vals = (double *)malloc(n * sizeof(double));
    int *counts = (int *)calloc(n, sizeof(int));
    int unique_count = unique_counts(x, n, &unique_vals);
    double mean = (double)n / unique_count;
    double sum_abs_diff = 0.0;
    for (int i = 0; i < unique_count; i++) {
        sum_abs_diff += fabs(counts[i] - mean);
    }
    free(unique_vals);
    free(counts);
    return 1 - (sum_abs_diff / (2 * mean * fmax(unique_count - 1, 1)));
}

double renyi_entropy(const double *x, int n, double alpha) {
    if (alpha < 0) {
        fprintf(stderr, "Parameter alpha should be positive!\n");
        exit(EXIT_FAILURE);
    }
    double *unique_vals = (double *)malloc(n * sizeof(double));
    int *counts = (int *)calloc(n, sizeof(int));
    int unique_count = unique_counts(x, n, &unique_vals);
    double sum = 0.0;
    for (int i = 0; i < unique_count; i++) {
        sum += pow(counts[i], alpha);
    }
    free(unique_vals);
    free(counts);
    if (alpha == 1) {
        return -sum * log2(sum);
    }
    return 1 / (1 - alpha) * log2(sum);
}

double negative_extropy(const double *x, int n) {
    double *unique_vals = (double *)malloc(n * sizeof(double));
    int *counts = (int *)calloc(n, sizeof(int));
    int unique_count = unique_counts(x, n, &unique_vals);
    double *p_inv = (double *)malloc(unique_count * sizeof(double));
    for (int i = 0; i < unique_count; i++) {
        p_inv[i] = 1.0 - counts[i];
    }
    double sum = 0.0;
    for (int i = 0; i < unique_count; i++) {
        sum += p_inv[i] * log2(p_inv[i]);
    }
    free(unique_vals);
    free(counts);
    free(p_inv);
    return -sum;
}

double mcintosh_d(const double *x, int n) {
    double *unique_vals = (double *)malloc(n * sizeof(double));
    int *counts = (int *)calloc(n, sizeof(int));
    int unique_count = unique_counts(x, n, &unique_vals);
    double sum_squares = 0.0;
    for (int i = 0; i < unique_count; i++) {
        sum_squares += counts[i] * counts[i];
    }
    free(unique_vals);
    free(counts);
    return (n - sqrt(sum_squares)) / (n - sqrt(n));
}
    free(unique_vals);
    free(counts);
    return 1 - sum_freq_squared;
}

double gibbs_m2(const double *x, int n) {
    double *unique_vals = (double *)malloc(n * sizeof(double));
    int *counts = (int *)calloc(n, sizeof(int));
    int unique_count = unique_counts(x, n, &unique_vals);
    
    double sum_freq_squared = 0.0;
    for (int i = 0; i < unique_count; i++) {
        double freq = (double)counts[i] / n;
        sum_freq_squared += freq * freq;
    }

    double k = (double)unique_count;
    free(unique_vals);
    free(counts);
    return (k / (k - 1)) * (1 - sum_freq_squared);
}
