/*
 * simulate.c — driver application for the probsim library.
 *
 * Walked through in the companion paper, paper/probsim.tex §6 ("The driver
 * and the tests").  Usage:
 *
 *     ./simulate [seed] [n]      (defaults: seed 20260718, n 200000)
 *
 * The program does four things:
 *
 *   1. Simulates n draws from each of the sixteen distributions
 *      and prints the sample mean/variance next to the theoretical
 *      values, so agreement is visible at a glance (law of large
 *      numbers in action, paper §6.1).
 *
 *   1b. Runs the balls-into-bins experiment of paper §4.6: Pearson's
 *      discrete Q statistic on uniform bin counts, repeated, so its
 *      convergence onto the continuous chi-squared law is visible —
 *      and prints even-df chi-squared tails next to the finite Poisson
 *      sums they equal (the count/wait duality).
 *
 *   2. Runs the two t-tests of src/stats.c on freshly simulated data:
 *      a one-sample test where H0 is TRUE (the p-value should look
 *      uniform-ish, i.e. usually > 0.05) and a Welch test where H0 is
 *      FALSE (the p-value should be tiny) — paper §6.2.
 *
 *   3. Runs Tarone's overdispersion test on simulated sibships, where H0
 *      (binomial: one coin for everyone) and H1 (beta-binomial: a coin
 *      per family) differ by a correlation of 0.002 — and shows that
 *      only registry-scale data can tell them apart (paper §3.6, §6.2).
 */
#include "probsim.h"

#include <stdio.h>
#include <stdlib.h>

/* One row of the moments table: name, theoretical moments, and n draws
 * produced by the callback.  Keeping the callback signature uniform lets
 * the table itself stay declarative (see rows[] in main). */
typedef double (*sampler_fn)(ps_rng *rng);

typedef struct {
    const char *name;
    double      mean;   /* theoretical E[X]   */
    double      var;    /* theoretical Var[X] */
    sampler_fn  draw;
} dist_row;

/* --- fixed-parameter wrappers so each table row is a plain function --- */
/* Parameters here deliberately match the worked examples in paper §3–§4. */
static double d_bernoulli(ps_rng *r)  { return ps_bernoulli(r, 0.3); }
static double d_binomial(ps_rng *r)   { return (double)ps_binomial(r, 20, 0.25); }
static double d_geometric(ps_rng *r)  { return (double)ps_geometric(r, 0.2); }
static double d_poisson(ps_rng *r)    { return (double)ps_poisson(r, 4.0); }
static double d_uniform(ps_rng *r)    { return ps_uniform(r, -1.0, 3.0); }
static double d_exponential(ps_rng *r){ return ps_exponential(r, 0.5); }
static double d_normal(ps_rng *r)     { return ps_normal(r, 10.0, 2.0); }
static double d_gamma(ps_rng *r)      { return ps_gamma(r, 3.0, 2.0); }
static double d_beta(ps_rng *r)       { return ps_beta(r, 2.0, 5.0); }
static double d_chi_squared(ps_rng *r){ return ps_chi_squared(r, 6.0); }
static double d_student_t(ps_rng *r)  { return ps_student_t(r, 8.0); }
static double d_f(ps_rng *r)          { return ps_f(r, 5.0, 12.0); }
static double d_negbin(ps_rng *r)     { return (double)ps_negative_binomial(r, 3, 0.4); }
static double d_betabin(ps_rng *r)    { return (double)ps_beta_binomial(r, 20, 2.0, 5.0); }
static double d_rayleigh(ps_rng *r)   { return ps_rayleigh(r, 2.0); }
static double d_gumbel(ps_rng *r)     { return ps_gumbel(r, 0.5, 2.0); }

static void print_moments_table(ps_rng *rng, size_t n, double *buf)
{
    /* Theoretical moments follow the formulas tabulated in paper §3–§4. */
    const dist_row rows[] = {
        { "Bernoulli(0.3)",       0.3,        0.21,       d_bernoulli   },
        { "Binomial(20, 0.25)",   5.0,        3.75,       d_binomial    },
        { "Geometric(0.2)",       4.0,        20.0,       d_geometric   },
        { "Poisson(4)",           4.0,        4.0,        d_poisson     },
        { "Uniform(-1, 3)",       1.0,        4.0 / 3.0,  d_uniform     },
        { "Exponential(0.5)",     2.0,        4.0,        d_exponential },
        { "Normal(10, 2)",        10.0,       4.0,        d_normal      },
        { "Gamma(3, 2)",          6.0,        12.0,       d_gamma       },
        { "Beta(2, 5)",           2.0 / 7.0,  10.0/392.0, d_beta        },
        /* Same mean as Binomial(20, 2/7), 3.4x the variance (paper §3.6). */
        { "BetaBinom(20, 2, 5)",  40.0 / 7.0, 5400.0/392.0, d_betabin   },
        { "ChiSquared(6)",        6.0,        12.0,       d_chi_squared },
        { "StudentT(8)",          0.0,        8.0 / 6.0,  d_student_t   },
        { "F(5, 12)",             1.2,        1.08,       d_f           },
        /* The hash modeler's annex (paper §3.5, §4.9–§4.10). */
        { "NegBinomial(3, 0.4)",  4.5,        11.25,      d_negbin      },
        { "Rayleigh(2)",          2.506628,   1.716815,   d_rayleigh    },
        { "Gumbel(0.5, 2)",       1.654431,   6.579736,   d_gumbel      },
    };
    const size_t n_rows = sizeof rows / sizeof rows[0];

    printf("Moments over %zu draws per distribution "
           "(sample vs. theory)\n\n", n);
    printf("%-20s %10s %10s   %10s %10s\n",
           "distribution", "mean", "E[X]", "variance", "Var[X]");
    printf("%-20s %10s %10s   %10s %10s\n",
           "------------", "----", "----", "--------", "------");

    for (size_t i = 0; i < n_rows; i++) {
        for (size_t j = 0; j < n; j++)
            buf[j] = rows[i].draw(rng);
        ps_summary s = ps_summarize(buf, n);
        printf("%-20s %10.4f %10.4f   %10.4f %10.4f\n",
               rows[i].name, s.mean, rows[i].mean, s.var, rows[i].var);
    }
}

/*
 * Balls into bins (paper §4.6): Pearson's Q on m equal bins is a
 * DISCRETE statistic — the counts are integers, so Q ranges over a
 * lattice of values — whose large-n law is the continuous chi-squared
 * with m-1 degrees of freedom, exactly as the discrete binomial's
 * large-n law is the continuous normal.  Repeat the experiment and the
 * discrete statistic is seen sitting on its continuous asymptote:
 * mean near m-1, variance near 2(m-1), 5% of runs beyond the 95% point.
 */
#define MAX_BINS 64

static void balls_into_bins_demo(ps_rng *rng, size_t experiments,
                                 size_t balls, size_t bins)
{
    long counts[MAX_BINS];
    if (bins < 2 || bins > MAX_BINS) return;

    const double expect = (double)balls / (double)bins;
    const double df = (double)bins - 1.0;
    double sum_q = 0.0, sum_q2 = 0.0;
    size_t rejections = 0;

    for (size_t e = 0; e < experiments; e++) {
        for (size_t i = 0; i < bins; i++) counts[i] = 0;
        for (size_t b = 0; b < balls; b++) {
            size_t i = (size_t)(ps_runif(rng) * (double)bins);
            counts[i >= bins ? bins - 1 : i]++;
        }
        double q = 0.0;
        for (size_t i = 0; i < bins; i++) {
            const double d = (double)counts[i] - expect;
            q += d * d / expect;
        }
        sum_q += q;
        sum_q2 += q * q;
        if (1.0 - ps_chi_squared_cdf(q, df) < 0.05) rejections++;
    }

    const double mean = sum_q / (double)experiments;
    const double var  = sum_q2 / (double)experiments - mean * mean;
    printf("  %4zu balls, %2zu bins   mean Q = %6.3f (df = %g)   "
           "var Q = %7.3f (2 df = %g)   reject at 5%%: %.1f%%\n",
           balls, bins, mean, df, var, 2.0 * df,
           100.0 * (double)rejections / (double)experiments);
}

/*
 * Sibship demo (paper §3.6, §5.4): simulate `families` sibships of
 * `size` children each, either from one common coin (rho = 0) or from a
 * per-family coin drawn from Beta(a, b), then ask Tarone's test whether
 * it can tell.  This is the Harvard-vs-Scandinavia disagreement in
 * miniature: at rho = 0.002 the answer depends entirely on how many
 * families you have.
 */
static void sibship_demo(ps_rng *rng, const char *label, size_t families,
                         long size, double p, double rho,
                         long *ks, long *ns)
{
    double a = 0.0, b = 0.0;
    const bool mixed = ps_beta_binomial_ab(p, rho, &a, &b);
    for (size_t i = 0; i < families; i++) {
        ns[i] = size;
        ks[i] = mixed ? ps_beta_binomial(rng, size, a, b)
                      : ps_binomial(rng, size, p);
    }
    ps_overdispersion od = ps_tarone_z(ks, ns, families);
    if (!od.ok) {
        printf("  %-38s (test undefined)\n", label);
        return;
    }
    printf("  %-38s rho_hat = %+8.5f   Z = %6.2f   p = %8.4g  %s\n",
           label, od.rho, od.z, od.p,
           od.p < 0.05 ? "<- detected" : "");
}

static void print_t_test(const char *title, ps_t_test t)
{
    printf("%s\n", title);
    if (!t.ok) {
        printf("  (test could not be run)\n\n");
        return;
    }
    printf("  estimate = %8.4f   se = %.4f\n", t.estimate, t.se);
    printf("  t = %8.4f   df = %8.2f   two-sided p = %.4g\n",
           t.t, t.df, t.p);
    printf("  => %s H0 at the 5%% level\n\n",
           t.p < 0.05 ? "REJECT" : "do not reject");
}

int main(int argc, char **argv)
{
    uint64_t seed = argc > 1 ? strtoull(argv[1], NULL, 10) : 20260718ULL;
    size_t   n    = argc > 2 ? (size_t)strtoull(argv[2], NULL, 10) : 200000;
    if (n < 2) {
        fprintf(stderr, "usage: %s [seed] [n>=2]\n", argv[0]);
        return EXIT_FAILURE;
    }

    ps_rng rng;
    ps_rng_seed(&rng, seed);

    double *buf = malloc(n * sizeof *buf);
    double *buf2 = malloc(n * sizeof *buf2);
    if (buf == NULL || buf2 == NULL) {
        fprintf(stderr, "out of memory\n");
        free(buf);
        free(buf2);
        return EXIT_FAILURE;
    }

    printf("probsim driver — seed %llu\n", (unsigned long long)seed);
    printf("======================================================"
           "==============\n\n");
    print_moments_table(&rng, n, buf);

    /* ---- discrete statistic, continuous law (paper §4.6) ---- */
    printf("\nBalls into bins: discrete Pearson Q vs its continuous "
           "chi-squared limit\n");
    printf("------------------------------------------------------"
           "--------------\n");
    printf("  10000 experiments each; Q is integer-built, "
           "chi-squared is its large-n law.\n\n");
    balls_into_bins_demo(&rng, 10000, 100, 2);   /* Q = squared std binomial */
    balls_into_bins_demo(&rng, 10000, 600, 6);   /* a fair die, 100 per face */
    balls_into_bins_demo(&rng, 10000, 2000, 20);

    /* The duality read the other way: even-df chi-squared tails ARE
     * finite Poisson sums (paper §4.6), so the continuous tail area and
     * the discrete CDF must agree to machine precision. */
    printf("\n  Even-df chi-squared tails as finite Poisson sums "
           "(x = 12.5916, the chi2_6 95%% point):\n");
    for (long k = 1; k <= 3; k++) {
        const double x = 12.591587243743977;
        printf("    1 - F_chi2(x, %ld df) = %.15f   "
               "Pr[Poisson(x/2) <= %ld] = %.15f\n",
               2 * (k + 1), 1.0 - ps_chi_squared_cdf(x, 2.0 * (double)(k + 1)),
               k, ps_poisson_cdf(k, x / 2.0));
    }

    /* ---- t-test demos on simulated data (paper §6.2) ---- */
    printf("\nt-tests on freshly simulated data (m = 30 per sample)\n");
    printf("-----------------------------------------------------\n\n");
    const size_t m = 30;

    /* (a) H0 true: X ~ N(10, 2), testing E[X] = 10. */
    for (size_t i = 0; i < m; i++)
        buf[i] = ps_normal(&rng, 10.0, 2.0);
    print_t_test("One-sample t-test: X ~ N(10, 2), H0: mean = 10 (H0 true)",
                 ps_t_test_one_sample(buf, m, 10.0));

    /* (b) H0 false: X ~ N(10, 2) vs Y ~ N(11.5, 3) — a real difference,
     * unequal variances, so Welch is the right tool. */
    for (size_t i = 0; i < m; i++)
        buf2[i] = ps_normal(&rng, 11.5, 3.0);
    print_t_test("Welch t-test: N(10, 2) vs N(11.5, 3), H0: equal means "
                 "(H0 false)",
                 ps_t_test_welch(buf, m, buf2, m));

    /* ---- one coin or one coin per family?  (paper §3.6, §5.4) ---- */
    printf("Tarone's test on simulated sibships of 4, p = 0.5122\n");
    printf("-----------------------------------------------------\n");
    printf("  H0: every child is a flip of the SAME coin (binomial).\n"
           "  H1: each family gets its own coin (beta-binomial).\n\n");
    {
        const size_t big = 1000000, small = 5000;
        long *ks = malloc(big * sizeof *ks);
        long *ns = malloc(big * sizeof *ns);
        if (ks != NULL && ns != NULL) {
            sibship_demo(&rng, "one coin, 1000000 families",
                         big, 4, 0.5122, 0.0, ks, ns);
            sibship_demo(&rng, "rho = 0.002, 5000 families",
                         small, 4, 0.5122, 0.002, ks, ns);
            sibship_demo(&rng, "rho = 0.002, 1000000 families",
                         big, 4, 0.5122, 0.002, ks, ns);
            sibship_demo(&rng, "rho = 0.05 (strong), 5000 families",
                         small, 4, 0.5122, 0.05, ks, ns);
            printf("\n  A correlation of 0.002 is real but nearly invisible:"
                   " Z grows as\n  sqrt(pairs), so it takes registry-scale"
                   " data to see it at all.\n");
        }
        free(ks);
        free(ns);
    }

    free(buf);
    free(buf2);
    return EXIT_SUCCESS;
}
