/*
 * dist.c — samplers, pmfs/pdfs and CDFs for the twelve classical
 * distributions.
 *
 * Explained distribution-by-distribution in the companion paper,
 * paper/probsim.tex §3 (discrete) and §4 (continuous); each function below
 * names its paper section.  Public interface: include/probsim.h.
 *
 * Design notes
 *   - Every sampler reduces to ps_runif() (src/rng.c), so a single seed
 *     reproduces an entire simulation.
 *   - Samplers favour the clearest correct classical algorithm over the
 *     fastest one; where a cleverer method is standard (Box–Muller polar,
 *     Marsaglia–Tsang gamma) the paper explains it and the bibliography
 *     holds the original article (refs/).
 *   - CDFs delegate to the special functions in src/special.c; that
 *     reduction (t, F, chi-squared -> incomplete beta/gamma) is the
 *     subject of paper §5.1.
 */
#include "probsim.h"

#include <math.h>

/* M_PI is POSIX, not C99, so -std=c99 -pedantic builds need our own copy. */
static const double PS_PI = 3.14159265358979323846;

/* ============================ Discrete ============================ */

/* Bernoulli(p) — paper §3.1. */
int ps_bernoulli(ps_rng *rng, double p)
{
    if (!(p >= 0.0 && p <= 1.0)) return -1;
    return ps_runif(rng) < p ? 1 : 0;
}

/*
 * Binomial(n, p) — paper §3.2.
 *
 * Simulated as what it *is*: n independent Bernoulli(p) trials.  O(n) per
 * draw, which is the right trade-off for a reference implementation; see
 * paper §3.2 for pointers to sub-linear samplers.
 */
long ps_binomial(ps_rng *rng, long n, double p)
{
    if (n < 0 || !(p >= 0.0 && p <= 1.0)) return -1;
    long k = 0;
    for (long i = 0; i < n; i++)
        k += (ps_runif(rng) < p);
    return k;
}

double ps_binomial_pmf(long k, long n, double p)
{
    if (n < 0 || !(p >= 0.0 && p <= 1.0)) return NAN;
    if (k < 0 || k > n) return 0.0;
    if (p == 0.0) return k == 0 ? 1.0 : 0.0;
    if (p == 1.0) return k == n ? 1.0 : 0.0;
    double lc = lgamma((double)n + 1.0) - lgamma((double)k + 1.0)
              - lgamma((double)(n - k) + 1.0);
    return exp(lc + (double)k * log(p) + (double)(n - k) * log1p(-p));
}

/*
 * Geometric(p) — paper §3.3.  Counts FAILURES before the first success
 * (support {0,1,2,...}).  Sampled by inversion of the CDF:
 * K = floor( log(1-U) / log(1-p) ).
 */
long ps_geometric(ps_rng *rng, double p)
{
    if (!(p > 0.0 && p <= 1.0)) return -1;
    if (p == 1.0) return 0;
    return (long)floor(log1p(-ps_runif(rng)) / log1p(-p));
}

double ps_geometric_pmf(long k, double p)
{
    if (!(p > 0.0 && p <= 1.0)) return NAN;
    if (k < 0) return 0.0;
    return p * exp((double)k * log1p(-p));
}

/*
 * Poisson(lambda) — paper §3.4.
 *
 * Knuth's product-of-uniforms method: count uniforms until their product
 * drops below exp(-lambda).  To keep that product away from underflow for
 * large lambda, we exploit additivity — Poisson(a+b) is the sum of
 * independent Poisson(a) and Poisson(b) — and draw in chunks of rate <= 30.
 * Cost is O(lambda) uniforms per draw, fine for a reference implementation.
 */
static long poisson_knuth(ps_rng *rng, double lambda)
{
    const double limit = exp(-lambda);
    long k = 0;
    double prod = ps_runif(rng);
    while (prod > limit) {
        k++;
        prod *= ps_runif(rng);
    }
    return k;
}

long ps_poisson(ps_rng *rng, double lambda)
{
    if (!(lambda >= 0.0)) return -1;
    const double chunk = 30.0;
    long k = 0;
    while (lambda > chunk) {
        k += poisson_knuth(rng, chunk);
        lambda -= chunk;
    }
    return k + poisson_knuth(rng, lambda);
}

double ps_poisson_pmf(long k, double lambda)
{
    if (!(lambda >= 0.0)) return NAN;
    if (k < 0) return 0.0;
    if (lambda == 0.0) return k == 0 ? 1.0 : 0.0;
    return exp((double)k * log(lambda) - lambda - lgamma((double)k + 1.0));
}

/*
 * Negative binomial(r, p) — paper §3.5.  Counts FAILURES before the
 * r-th success; the geometric of §3.3 is the case r = 1, and the
 * sampler is that identity: a sum of r independent geometric waits.
 * Governs the expected work of multi-collision searches against an
 * ideal hash (paper §3.5) and overdispersed counts everywhere else.
 */
long ps_negative_binomial(ps_rng *rng, long r, double p)
{
    if (r < 1 || !(p > 0.0 && p <= 1.0)) return -1;
    long k = 0;
    for (long i = 0; i < r; i++)
        k += ps_geometric(rng, p);
    return k;
}

double ps_negative_binomial_pmf(long k, long r, double p)
{
    if (r < 1 || !(p > 0.0 && p <= 1.0)) return NAN;
    if (k < 0) return 0.0;
    double lc = lgamma((double)(k + r)) - lgamma((double)k + 1.0)
              - lgamma((double)r);
    return exp(lc + (double)r * log(p) + (double)k * log1p(-p));
}

/* P(X <= k) = I_p(r, k+1) — the incomplete-beta reduction again. */
double ps_negative_binomial_cdf(long k, long r, double p)
{
    if (r < 1 || !(p > 0.0 && p <= 1.0)) return NAN;
    if (k < 0) return 0.0;
    return ps_incbeta(p, (double)r, (double)k + 1.0);
}

/* =========================== Continuous =========================== */

/* Uniform(a, b) — paper §4.1. */
double ps_uniform(ps_rng *rng, double a, double b)
{
    if (!(a <= b)) return NAN;
    return a + (b - a) * ps_runif(rng);
}

/* Exponential(lambda), rate parameterisation — paper §4.2.
 * Inversion: X = -log(1-U)/lambda. */
double ps_exponential(ps_rng *rng, double lambda)
{
    if (!(lambda > 0.0)) return NAN;
    return -log1p(-ps_runif(rng)) / lambda;
}

double ps_exponential_cdf(double x, double lambda)
{
    if (!(lambda > 0.0)) return NAN;
    return x <= 0.0 ? 0.0 : -expm1(-lambda * x);
}

/*
 * Normal(mu, sigma) — paper §4.3.
 *
 * Marsaglia's polar variant of the Box–Muller transform (Box & Muller
 * 1958, refs/box-muller-1958.pdf): draw (v1, v2) uniform on the unit
 * disc, then v1 * sqrt(-2 log s / s) is exactly standard normal.  The
 * method produces two independent normals; we keep one and discard the
 * spare so the sampler stays stateless (see paper §4.3, footnote).
 */
double ps_normal(ps_rng *rng, double mu, double sigma)
{
    if (!(sigma >= 0.0)) return NAN;
    double v1, v2, s;
    do {
        v1 = 2.0 * ps_runif(rng) - 1.0;
        v2 = 2.0 * ps_runif(rng) - 1.0;
        s = v1 * v1 + v2 * v2;
    } while (s >= 1.0 || s == 0.0);
    return mu + sigma * v1 * sqrt(-2.0 * log(s) / s);
}

double ps_normal_pdf(double x, double mu, double sigma)
{
    if (!(sigma > 0.0)) return NAN;
    double z = (x - mu) / sigma;
    return exp(-0.5 * z * z) / (sigma * sqrt(2.0 * PS_PI));
}

double ps_normal_cdf(double x, double mu, double sigma)
{
    if (!(sigma > 0.0)) return NAN;
    return 0.5 * erfc(-(x - mu) / (sigma * sqrt(2.0)));
}

/*
 * Gamma(shape k, scale theta) — paper §4.4.
 *
 * Marsaglia & Tsang's (2000) squeeze method for shape >= 1: with
 * d = k - 1/3 and c = 1/sqrt(9d), accept d*(1 + c*Z)^3 for standard
 * normal Z with a cheap squeeze then an exact log test.  For shape < 1,
 * boost via  Gamma(k) = Gamma(k+1) * U^(1/k).
 */
static double gamma_marsaglia_tsang(ps_rng *rng, double shape)
{
    const double d = shape - 1.0 / 3.0;
    const double c = 1.0 / sqrt(9.0 * d);
    for (;;) {
        double z, v;
        do {
            z = ps_normal(rng, 0.0, 1.0);
            v = 1.0 + c * z;
        } while (v <= 0.0);
        v = v * v * v;
        double u = ps_runif(rng);
        if (u < 1.0 - 0.0331 * z * z * z * z)               /* squeeze  */
            return d * v;
        if (log(u) < 0.5 * z * z + d * (1.0 - v + log(v)))  /* exact    */
            return d * v;
    }
}

double ps_gamma(ps_rng *rng, double shape, double scale)
{
    if (!(shape > 0.0) || !(scale > 0.0)) return NAN;
    if (shape < 1.0) {
        /* Boost: draw at shape+1 and multiply by U^(1/shape). */
        double u = ps_runif(rng);
        while (u == 0.0)  /* avoid log(0) / 0^(1/shape) edge */
            u = ps_runif(rng);
        return scale * gamma_marsaglia_tsang(rng, shape + 1.0)
                     * pow(u, 1.0 / shape);
    }
    return scale * gamma_marsaglia_tsang(rng, shape);
}

double ps_gamma_cdf(double x, double shape, double scale)
{
    if (!(shape > 0.0) || !(scale > 0.0)) return NAN;
    return x <= 0.0 ? 0.0 : ps_incgamma_lower(shape, x / scale);
}

/* Beta(a, b) — paper §4.5.  X = G1 / (G1 + G2) with independent
 * G1 ~ Gamma(a, 1), G2 ~ Gamma(b, 1). */
double ps_beta(ps_rng *rng, double a, double b)
{
    if (!(a > 0.0) || !(b > 0.0)) return NAN;
    double g1 = ps_gamma(rng, a, 1.0);
    double g2 = ps_gamma(rng, b, 1.0);
    return g1 / (g1 + g2);
}

double ps_beta_cdf(double x, double a, double b)
{
    return ps_incbeta(x, a, b);
}

/* Chi-squared(k) — paper §4.6.  Exactly Gamma(k/2, 2). */
double ps_chi_squared(ps_rng *rng, double k)
{
    if (!(k > 0.0)) return NAN;
    return ps_gamma(rng, k / 2.0, 2.0);
}

double ps_chi_squared_cdf(double x, double k)
{
    if (!(k > 0.0)) return NAN;
    return x <= 0.0 ? 0.0 : ps_incgamma_lower(k / 2.0, x / 2.0);
}

/* Student's t(nu) — paper §4.7.  T = Z / sqrt(V/nu) from its very
 * definition (Student 1908, refs/student-1908.pdf). */
double ps_student_t(ps_rng *rng, double nu)
{
    if (!(nu > 0.0)) return NAN;
    double z = ps_normal(rng, 0.0, 1.0);
    double v = ps_chi_squared(rng, nu);
    return z / sqrt(v / nu);
}

double ps_student_t_pdf(double x, double nu)
{
    if (!(nu > 0.0)) return NAN;
    double lc = lgamma((nu + 1.0) / 2.0) - lgamma(nu / 2.0)
              - 0.5 * log(nu * PS_PI);
    return exp(lc - 0.5 * (nu + 1.0) * log1p(x * x / nu));
}

/* P(T <= x) via the incomplete beta reduction — paper §5.1, eq. (5.2). */
double ps_student_t_cdf(double x, double nu)
{
    if (!(nu > 0.0)) return NAN;
    double tail = 0.5 * ps_incbeta(nu / (nu + x * x), nu / 2.0, 0.5);
    return x >= 0.0 ? 1.0 - tail : tail;
}

/* Fisher's F(d1, d2) — paper §4.8. */
double ps_f(ps_rng *rng, double d1, double d2)
{
    if (!(d1 > 0.0) || !(d2 > 0.0)) return NAN;
    double u = ps_chi_squared(rng, d1);
    double v = ps_chi_squared(rng, d2);
    return (u / d1) / (v / d2);
}

double ps_f_cdf(double x, double d1, double d2)
{
    if (!(d1 > 0.0) || !(d2 > 0.0)) return NAN;
    if (x <= 0.0) return 0.0;
    return ps_incbeta(d1 * x / (d1 * x + d2), d1 / 2.0, d2 / 2.0);
}

/*
 * Rayleigh(sigma) — paper §4.9.  Length of a 2-D Gaussian vector;
 * equivalently sigma * sqrt(2E) for E ~ Exp(1), which is the inversion
 * sampler below.  Scaled by sqrt(N), it is the limit law of the
 * birthday-collision waiting time of an ideal N-bin hash.
 */
double ps_rayleigh(ps_rng *rng, double sigma)
{
    if (!(sigma > 0.0)) return NAN;
    return sigma * sqrt(-2.0 * log1p(-ps_runif(rng)));
}

double ps_rayleigh_pdf(double x, double sigma)
{
    if (!(sigma > 0.0)) return NAN;
    if (x <= 0.0) return 0.0;
    double s2 = sigma * sigma;
    return (x / s2) * exp(-0.5 * x * x / s2);
}

double ps_rayleigh_cdf(double x, double sigma)
{
    if (!(sigma > 0.0)) return NAN;
    return x <= 0.0 ? 0.0 : -expm1(-0.5 * x * x / (sigma * sigma));
}

/*
 * Gumbel(mu, beta) — paper §4.10.  The extreme-value limit for maxima
 * of light-tailed samples; sampled by inversion of the double
 * exponential, X = mu - beta * ln(-ln U).
 */
double ps_gumbel(ps_rng *rng, double mu, double beta)
{
    if (!(beta > 0.0)) return NAN;
    double u = ps_runif(rng);
    while (u == 0.0)          /* log(0) guard; U = 0 has measure zero */
        u = ps_runif(rng);
    return mu - beta * log(-log(u));
}

double ps_gumbel_pdf(double x, double mu, double beta)
{
    if (!(beta > 0.0)) return NAN;
    double z = (x - mu) / beta;
    return exp(-(z + exp(-z))) / beta;
}

double ps_gumbel_cdf(double x, double mu, double beta)
{
    if (!(beta > 0.0)) return NAN;
    return exp(-exp(-(x - mu) / beta));
}
