/*
 * probsim.h — public API of the probsim library.
 *
 * A small C99 library for simulating the major classical probability
 * distributions and running the basic inference that connects them
 * (one-sample and Welch two-sample t-tests).
 *
 * The library is documented, distribution by distribution, in the companion
 * paper paper/probsim.tex (built to paper/probsim.pdf).  Every declaration
 * below names the paper section that explains it; the paper in turn cites
 * the implementing function, so the two can be read side by side.
 *
 * Conventions
 *   - All samplers take a ps_rng * first, so streams are explicit and the
 *     library is reentrant; there is no hidden global state.
 *   - Continuous parameters follow the usual textbook conventions noted at
 *     each declaration (e.g. the exponential and gamma use RATE lambda and
 *     SHAPE/SCALE respectively).
 *   - Functions never abort: invalid parameters yield NaN (continuous),
 *     -1 (discrete counts), or a false `ok` flag (tests).
 */
#ifndef PROBSIM_H
#define PROBSIM_H

#include <stddef.h>
#include <stdint.h>
#include <stdbool.h>

/* ------------------------------------------------------------------ *
 *  Random number generation                       (paper §2, src/rng.c)
 * ------------------------------------------------------------------ */

/*
 * xoshiro256++ state (Blackman & Vigna 2019, refs/blackman-vigna-xoshiro.pdf).
 * 256 bits of state, period 2^256 - 1, passes BigCrush.
 */
typedef struct {
    uint64_t s[4];
} ps_rng;

/* Seed deterministically from any 64-bit value via splitmix64. */
void     ps_rng_seed(ps_rng *rng, uint64_t seed);
/* Next raw 64-bit output. */
uint64_t ps_rng_next(ps_rng *rng);
/* Uniform double on [0,1) with 53 random bits — the basis of every sampler. */
double   ps_runif(ps_rng *rng);

/* ------------------------------------------------------------------ *
 *  Discrete distributions — samplers and pmf/cdf  (paper §3, src/dist.c)
 *  Samplers return -1 on invalid parameters.
 * ------------------------------------------------------------------ */

/* Bernoulli(p): 1 with probability p, else 0.                 (paper §3.1) */
int    ps_bernoulli(ps_rng *rng, double p);

/* Binomial(n, p): number of successes in n Bernoulli trials.  (paper §3.2) */
long   ps_binomial(ps_rng *rng, long n, double p);
double ps_binomial_pmf(long k, long n, double p);

/* Geometric(p): number of FAILURES before the first success,
 * support {0, 1, 2, ...} (the convention used by R's rgeom).   (paper §3.3) */
long   ps_geometric(ps_rng *rng, double p);
double ps_geometric_pmf(long k, double p);

/* Poisson(lambda): event counts at constant rate.             (paper §3.4) */
long   ps_poisson(ps_rng *rng, double lambda);
double ps_poisson_pmf(long k, double lambda);

/* Negative binomial(r, p): failures before the r-th success
 * (generalizes the geometric, which is r = 1).                (paper §3.5) */
long   ps_negative_binomial(ps_rng *rng, long r, double p);
double ps_negative_binomial_pmf(long k, long r, double p);
double ps_negative_binomial_cdf(long k, long r, double p);

/* ------------------------------------------------------------------ *
 *  Continuous distributions — samplers, pdf, cdf  (paper §4, src/dist.c)
 *  Samplers/pdf/cdf return NaN on invalid parameters.
 * ------------------------------------------------------------------ */

/* Uniform(a, b) on [a, b).                                    (paper §4.1) */
double ps_uniform(ps_rng *rng, double a, double b);

/* Exponential(lambda), RATE parameterisation: mean 1/lambda.  (paper §4.2) */
double ps_exponential(ps_rng *rng, double lambda);
double ps_exponential_cdf(double x, double lambda);

/* Normal(mu, sigma) via the Marsaglia polar variant of the
 * Box–Muller transform (refs/box-muller-1958.pdf).            (paper §4.3) */
double ps_normal(ps_rng *rng, double mu, double sigma);
double ps_normal_pdf(double x, double mu, double sigma);
double ps_normal_cdf(double x, double mu, double sigma);

/* Gamma(shape k, scale theta), mean k*theta, via the
 * Marsaglia–Tsang (2000) squeeze method.                      (paper §4.4) */
double ps_gamma(ps_rng *rng, double shape, double scale);
double ps_gamma_cdf(double x, double shape, double scale);

/* Beta(a, b) as G1/(G1+G2) with independent gammas.           (paper §4.5) */
double ps_beta(ps_rng *rng, double a, double b);
double ps_beta_cdf(double x, double a, double b);

/* Chi-squared(k) = Gamma(k/2, 2).                             (paper §4.6) */
double ps_chi_squared(ps_rng *rng, double k);
double ps_chi_squared_cdf(double x, double k);

/* Student's t(nu) = Z / sqrt(V/nu), Z~N(0,1), V~chi2(nu).     (paper §4.7) */
double ps_student_t(ps_rng *rng, double nu);
double ps_student_t_pdf(double x, double nu);
double ps_student_t_cdf(double x, double nu);

/* Fisher's F(d1, d2) = (U/d1)/(V/d2), U,V independent chi2.   (paper §4.8) */
double ps_f(ps_rng *rng, double d1, double d2);
double ps_f_cdf(double x, double d1, double d2);

/* Rayleigh(sigma): length of a 2-D Gaussian vector; governs the
 * birthday-collision waiting time of an ideal hash.           (paper §4.9) */
double ps_rayleigh(ps_rng *rng, double sigma);
double ps_rayleigh_pdf(double x, double sigma);
double ps_rayleigh_cdf(double x, double sigma);

/* Gumbel(mu, beta): limit law of maxima of light-tailed draws;
 * longest runs and fullest hash buckets.                     (paper §4.10) */
double ps_gumbel(ps_rng *rng, double mu, double beta);
double ps_gumbel_pdf(double x, double mu, double beta);
double ps_gumbel_cdf(double x, double mu, double beta);

/* ------------------------------------------------------------------ *
 *  Special functions                            (paper §5.1, src/special.c)
 *  These power the t, F, beta, gamma and chi-squared CDFs above.
 * ------------------------------------------------------------------ */

/* Regularized incomplete beta I_x(a, b), by Lentz continued fraction. */
double ps_incbeta(double x, double a, double b);
/* Regularized lower incomplete gamma P(a, x), series + continued fraction. */
double ps_incgamma_lower(double a, double x);

/* ------------------------------------------------------------------ *
 *  Descriptive statistics and t-tests            (paper §5, src/stats.c)
 * ------------------------------------------------------------------ */

typedef struct {
    size_t n;
    double mean;
    double var;      /* unbiased sample variance (divisor n-1) */
    double sd;
    double min;
    double max;
} ps_summary;

/* One-pass Welford summary of x[0..n-1]; n==0 gives NaN fields. */
ps_summary ps_summarize(const double *x, size_t n);

typedef struct {
    bool   ok;        /* false if inputs were unusable (e.g. n < 2)   */
    double t;         /* the t statistic                              */
    double df;        /* degrees of freedom (fractional under Welch)  */
    double p;         /* two-sided p-value from ps_student_t_cdf      */
    double estimate;  /* mean (one-sample) or mean difference (Welch) */
    double se;        /* standard error of the estimate               */
} ps_t_test;

/* One-sample t-test of H0: E[X] = mu0.                        (paper §5.2) */
ps_t_test ps_t_test_one_sample(const double *x, size_t n, double mu0);

/* Welch's two-sample t-test of H0: E[X] = E[Y], unequal
 * variances, Welch–Satterthwaite df (Welch 1947).             (paper §5.3) */
ps_t_test ps_t_test_welch(const double *x, size_t nx,
                          const double *y, size_t ny);

#endif /* PROBSIM_H */
