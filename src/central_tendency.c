#include "central_tendency.h"

double mean(double* x, int n) {
    double sum = 0.0;
    for (int i = 0; i < n; i++) {
        sum += x[i];
    }
    return sum / n;
}

double median(double* x, int n) {
    qsort(x, n, sizeof(double), compare);
    if (n % 2 == 0) {
        return (x[n / 2] + x[n / 2 - 1]) / 2.0;
    } else {
        return x[n / 2];
    }
}

int compare(const void* a, const void* b) {
    return (*(double*)a - *(double*)b);
}

double* mode(double* x, int n, int* mode_count) {
    int* frequency = (int*)calloc(n, sizeof(int));
    double* high_list = (double*)malloc(n * sizeof(double));
    int hf = 0, count = 0;

    for (int i = 0; i < n; i++) {
        frequency[(int)x[i]]++;
        if (frequency[(int)x[i]] > hf) {
            hf = frequency[(int)x[i]];
        }
    }

    for (int i = 0; i < n; i++) {
        if (frequency[(int)x[i]] == hf) {
            high_list[count++] = x[i];
        }
    }

    free(frequency);
    *mode_count = count;
    return high_list;
}

double midrange(double* x, int n) {
    double maximum = -DBL_MAX;
    double minimum = DBL_MAX;
    for (int i = 0; i < n; i++) {
        if (x[i] > maximum) maximum = x[i];
        if (x[i] < minimum) minimum = x[i];
    }
    return (maximum + minimum) / 2.0;
}

double midhinge(double* x, int n) {
    double q1 = quantile(x, n, 0.25);
    double q3 = quantile(x, n, 0.75);
    return (q1 + q3) / 2.0;
}

double trimean(double* x, int n) {
    double q1 = quantile(x, n, 0.25);
    double q2 = quantile(x, n, 0.5);
    double q3 = quantile(x, n, 0.75);
    return 0.5 * q2 + 0.25 * q1 + 0.25 * q3;
}

double contraharmonic_mean(double* x, int n) {
    double sum_sq = 0.0;
    double sum = 0.0;
    for (int i = 0; i < n; i++) {
        sum_sq += x[i] * x[i];
        sum += x[i];
    }
    return sum_sq / sum;
}

double midmean(double* x, int n) {
    double q1 = quantile(x, n, 0.25);
    double q3 = quantile(x, n, 0.75);
    double sum = 0.0;
    int count = 0;
    for (int i = 0; i < n; i++) {
        if (x[i] >= q1 && x[i] <= q3) {
            sum += x[i];
            count++;
        }
    }
    return sum / count;
}

double hodges_lehmann_sen_location(double* x, int n) {
    double* product = (double*)malloc(n * n * sizeof(double));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            product[i * n + j] = x[i] + x[j];
        }
    }
    double result = median(product, n * n);
    free(product);
    return result;
}

double standard_trimmed_harrell_davis_quantile(double* x, int n, double q) {
    if (q <= 0 || q >= 1) {
        fprintf(stderr, "Parameter q should be in range (0, 1).\n");
        exit(EXIT_FAILURE);
    }
    qsort(x, n, sizeof(double), compare);
    int n_calculated = (int)(1.0 / sqrt(n));
    double a = (n + 1) * q;
    double b = (n + 1) * (1.0 - q);
    // Further implementation needed for beta.cdf
    return 0.0; // Placeholder
}

double half_sample_mode(double* x, int n) {
    qsort(x, n, sizeof(double), compare);
    int corner_cases[] = {4, 3};
    while (n >= corner_cases[0]) {
        int half_y = n / 2;
        double w_min = x[n - 1] - x[0];
        int j = 0;
        for (int i = 0; i < n - half_y; i++) {
            double w = x[i + half_y - 1] - x[i];
            if (w <= w_min) {
                w_min = w;
                j = i;
            }
        }
        if (w_min == 0) {
            return x[j];
        }
        n = half_y - 1;
    }
    if (n == corner_cases[1]) {
        double z = 2 * x[1] - x[0] - x[2];
        if (z < 0) {
            return mean(x, 1);
        }
        if (z > 0) {
            return mean(x + 1, 1);
        }
        return x[1];
    }
    return mean(x, n);
}
