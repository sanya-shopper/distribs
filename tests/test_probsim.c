/*
 * test_probsim.c — unit tests for the probsim library.
 *
 * Described in the companion paper, docs/probsim.tex §6.3.  Run with
 * `make test`.  Three layers of checking:
 *
 *   1. EXACT-ISH fixtures: CDF/pmf/pdf values and t-test results are
 *      compared against reference numbers computed independently with
 *      SciPy 1.x (each fixture notes the generating call), to ~1e-9.
 *
 *   2. IDENTITIES: symmetry of the incomplete beta, the erf form of the
 *      half-shape incomplete gamma, CDF edge cases, RNG determinism.
 *
 *   3. STATISTICAL sanity: with a fixed seed, sample moments of every
 *      distribution must land within 5 standard errors of theory —
 *      deterministic in practice, yet sensitive to real sampler bugs.
 *
 * The test harness is deliberately tiny: CHECK(cond) and
 * CHECK_CLOSE(got, want, tol) macros that count failures and print
 * context on mismatch.
 */
#include "probsim.h"

#include <math.h>
#include <stdio.h>
#include <string.h>

static int checks = 0;
static int failures = 0;

#define CHECK(cond)                                                        \
    do {                                                                   \
        checks++;                                                          \
        if (!(cond)) {                                                     \
            failures++;                                                    \
            printf("FAIL %s:%d  %s\n", __FILE__, __LINE__, #cond);         \
        }                                                                  \
    } while (0)

#define CHECK_CLOSE(got, want, tol)                                        \
    do {                                                                   \
        checks++;                                                          \
        double g_ = (got), w_ = (want);                                    \
        if (!(fabs(g_ - w_) <= (tol))) {                                   \
            failures++;                                                    \
            printf("FAIL %s:%d  %s = %.12g, want %.12g (tol %g)\n",        \
                   __FILE__, __LINE__, #got, g_, w_, (double)(tol));       \
        }                                                                  \
    } while (0)

/* ---- 1. RNG basics: determinism and uniform range (paper §2) ---- */
static void test_rng(void)
{
    ps_rng a, b;
    ps_rng_seed(&a, 42);
    ps_rng_seed(&b, 42);
    for (int i = 0; i < 1000; i++)
        CHECK(ps_rng_next(&a) == ps_rng_next(&b));

    ps_rng_seed(&a, 7);
    for (int i = 0; i < 100000; i++) {
        double u = ps_runif(&a);
        CHECK(u >= 0.0 && u < 1.0);
    }
}

/* ---- 2. Fixtures: pmf/pdf values vs SciPy (paper §3–§4) ---- */
static void test_pmf_pdf_fixtures(void)
{
    /* scipy.stats.binom.pmf(5, 20, 0.25) */
    CHECK_CLOSE(ps_binomial_pmf(5, 20, 0.25), 0.202331151856924, 1e-12);
    /* scipy.stats.poisson.pmf(3, 4) */
    CHECK_CLOSE(ps_poisson_pmf(3, 4.0), 0.195366814813165, 1e-12);
    /* geometric: p (1-p)^k = 0.2 * 0.8^2 */
    CHECK_CLOSE(ps_geometric_pmf(2, 0.2), 0.128, 1e-12);
    /* scipy.stats.norm.pdf(1.0) */
    CHECK_CLOSE(ps_normal_pdf(1.0, 0.0, 1.0), 0.241970724519143, 1e-12);
    /* scipy.stats.t.pdf(1.0, 8) */
    CHECK_CLOSE(ps_student_t_pdf(1.0, 8.0), 0.227607580145303, 1e-12);

    /* Support edges. */
    CHECK(ps_binomial_pmf(-1, 20, 0.25) == 0.0);
    CHECK(ps_binomial_pmf(21, 20, 0.25) == 0.0);
    CHECK(ps_poisson_pmf(-1, 4.0) == 0.0);
}

/* ---- 3. Fixtures: CDF values vs SciPy (paper §5.1) ---- */
static void test_cdf_fixtures(void)
{
    /* scipy.stats.norm.cdf(1.96) */
    CHECK_CLOSE(ps_normal_cdf(1.96, 0.0, 1.0), 0.975002104851780, 1e-12);
    /* scipy.stats.t.cdf(2.2281388519649385, 10) — the 97.5% point */
    CHECK_CLOSE(ps_student_t_cdf(2.2281388519649385, 10.0), 0.975, 1e-9);
    /* scipy.stats.t.cdf(-1.5, 5) */
    CHECK_CLOSE(ps_student_t_cdf(-1.5, 5.0), 0.096951840121237, 1e-12);
    /* scipy.stats.chi2.cdf(3.8414588206941236, 1) — the 95% point */
    CHECK_CLOSE(ps_chi_squared_cdf(3.8414588206941236, 1.0), 0.95, 1e-9);
    /* scipy.stats.chi2.cdf(12.591587243743977, 6) — the 95% point */
    CHECK_CLOSE(ps_chi_squared_cdf(12.591587243743977, 6.0), 0.95, 1e-9);
    /* scipy.stats.f.cdf(4.964602743730711, 1, 10) — the 95% point */
    CHECK_CLOSE(ps_f_cdf(4.964602743730711, 1.0, 10.0), 0.95, 1e-9);
    /* scipy.stats.f.cdf(2.5, 5, 12) */
    CHECK_CLOSE(ps_f_cdf(2.5, 5.0, 12.0), 0.910175846395065, 1e-12);
    /* scipy.stats.beta.cdf(0.4, 2, 5) — also exact: 1-(1-x)^5(1+5x+15x^2)... */
    CHECK_CLOSE(ps_beta_cdf(0.4, 2.0, 5.0), 0.766720000000000, 1e-12);
    /* scipy.stats.gamma.cdf(5.0, 3, scale=2) */
    CHECK_CLOSE(ps_gamma_cdf(5.0, 3.0, 2.0), 0.456186884116670, 1e-12);
    /* scipy.stats.expon.cdf(1.2, scale=2) — rate 0.5 */
    CHECK_CLOSE(ps_exponential_cdf(1.2, 0.5), 0.451188363905974, 1e-12);

    /* Edges and invalid parameters. */
    CHECK(ps_exponential_cdf(-1.0, 0.5) == 0.0);
    CHECK(ps_chi_squared_cdf(0.0, 6.0) == 0.0);
    CHECK(isnan(ps_normal_cdf(0.0, 0.0, -1.0)));
    CHECK(isnan(ps_gamma_cdf(1.0, -2.0, 1.0)));
}

/* ---- 4. Identities among the special functions (paper §5.1) ---- */
static void test_special_identities(void)
{
    /* Symmetry: I_x(a,b) + I_{1-x}(b,a) = 1. */
    const double xs[] = { 0.05, 0.3, 0.5, 0.7, 0.95 };
    for (size_t i = 0; i < sizeof xs / sizeof xs[0]; i++)
        CHECK_CLOSE(ps_incbeta(xs[i], 2.5, 4.0)
                    + ps_incbeta(1.0 - xs[i], 4.0, 2.5), 1.0, 1e-12);

    /* Half-shape incomplete gamma is erf: P(1/2, x) = erf(sqrt(x)). */
    for (double x = 0.25; x <= 4.0; x += 0.75)
        CHECK_CLOSE(ps_incgamma_lower(0.5, x), erf(sqrt(x)), 1e-12);

    /* Shape-1 incomplete gamma is the exponential CDF: P(1,x) = 1-e^-x. */
    CHECK_CLOSE(ps_incgamma_lower(1.0, 2.0), -expm1(-2.0), 1e-12);

    /* CDF bounds. */
    CHECK(ps_incbeta(0.0, 2.0, 3.0) == 0.0);
    CHECK(ps_incbeta(1.0, 2.0, 3.0) == 1.0);
}

/* ---- 5. t-tests vs SciPy fixtures (paper §5.2–§5.3) ---- */
static void test_t_tests(void)
{
    const double x[] = { 5.1, 4.9, 6.2, 5.7, 5.5, 4.8, 5.9, 6.0, 5.2, 5.4 };
    const double y[] = { 6.3, 5.8, 6.1, 6.5, 5.9, 6.4,
                         6.2, 5.7, 6.0, 6.6, 6.1, 5.9 };
    const size_t nx = sizeof x / sizeof x[0];
    const size_t ny = sizeof y / sizeof y[0];

    /* scipy.stats.ttest_1samp(x, 5.0):
     * t = 3.121027607504, p = 0.012298238977 */
    ps_t_test one = ps_t_test_one_sample(x, nx, 5.0);
    CHECK(one.ok);
    CHECK_CLOSE(one.estimate, 5.47, 1e-12);
    CHECK_CLOSE(one.se, 0.150591426641, 1e-9);
    CHECK_CLOSE(one.t, 3.121027607504, 1e-9);
    CHECK_CLOSE(one.df, 9.0, 0.0);
    CHECK_CLOSE(one.p, 0.012298238977, 1e-9);

    /* scipy.stats.ttest_ind(x, y, equal_var=False):
     * t = -3.822418627849, p = 0.001847149818, Welch df = 14.087037411612 */
    ps_t_test welch = ps_t_test_welch(x, nx, y, ny);
    CHECK(welch.ok);
    CHECK_CLOSE(welch.estimate, -0.655, 1e-12);
    CHECK_CLOSE(welch.se, 0.171357473833, 1e-9);
    CHECK_CLOSE(welch.t, -3.822418627849, 1e-9);
    CHECK_CLOSE(welch.df, 14.087037411612, 1e-9);
    CHECK_CLOSE(welch.p, 0.001847149818, 1e-9);

    /* Degenerate inputs must be flagged, not crash. */
    CHECK(!ps_t_test_one_sample(x, 1, 5.0).ok);
    CHECK(!ps_t_test_one_sample(NULL, 10, 5.0).ok);
    const double flat[] = { 2.0, 2.0, 2.0 };
    CHECK(!ps_t_test_one_sample(flat, 3, 2.0).ok);
    CHECK(!ps_t_test_welch(x, nx, flat, 1).ok);
}

/* ---- 6. Sample moments within 5 standard errors of theory (paper §6.3) --- */

/* Draw n samples with `draw`, then require the sample mean to be within
 * 5 * sqrt(var/n) of `mean` and the sample variance within 12% of `var`
 * (a loose but bug-catching bound for these sample sizes). */
typedef double (*draw_fn)(ps_rng *rng);

static void check_moments(const char *name, draw_fn draw,
                          double mean, double var, ps_rng *rng)
{
    enum { N = 200000 };
    static double buf[N];
    for (size_t i = 0; i < N; i++)
        buf[i] = draw(rng);
    ps_summary s = ps_summarize(buf, N);

    checks += 2;
    double se = sqrt(var / (double)N);
    if (fabs(s.mean - mean) > 5.0 * se) {
        failures++;
        printf("FAIL moments(%s): mean %.6g, want %.6g +/- %.2g\n",
               name, s.mean, mean, 5.0 * se);
    }
    if (fabs(s.var - var) > 0.12 * var) {
        failures++;
        printf("FAIL moments(%s): var %.6g, want %.6g +/- 12%%\n",
               name, s.var, var);
    }
}

/* Same fixed parameter choices as the driver's table (app/simulate.c). */
static double m_bernoulli(ps_rng *r)  { return ps_bernoulli(r, 0.3); }
static double m_binomial(ps_rng *r)   { return (double)ps_binomial(r, 20, 0.25); }
static double m_geometric(ps_rng *r)  { return (double)ps_geometric(r, 0.2); }
static double m_poisson(ps_rng *r)    { return (double)ps_poisson(r, 4.0); }
static double m_poisson_big(ps_rng *r){ return (double)ps_poisson(r, 100.0); }
static double m_uniform(ps_rng *r)    { return ps_uniform(r, -1.0, 3.0); }
static double m_exponential(ps_rng *r){ return ps_exponential(r, 0.5); }
static double m_normal(ps_rng *r)     { return ps_normal(r, 10.0, 2.0); }
static double m_gamma(ps_rng *r)      { return ps_gamma(r, 3.0, 2.0); }
static double m_gamma_small(ps_rng *r){ return ps_gamma(r, 0.5, 1.0); }
static double m_beta(ps_rng *r)       { return ps_beta(r, 2.0, 5.0); }
static double m_chi_squared(ps_rng *r){ return ps_chi_squared(r, 6.0); }
static double m_student_t(ps_rng *r)  { return ps_student_t(r, 8.0); }
static double m_f(ps_rng *r)          { return ps_f(r, 5.0, 12.0); }

static void test_sampler_moments(void)
{
    ps_rng rng;
    ps_rng_seed(&rng, 20260718);

    check_moments("Bernoulli(0.3)",     m_bernoulli,   0.3,     0.21,    &rng);
    check_moments("Binomial(20,0.25)",  m_binomial,    5.0,     3.75,    &rng);
    check_moments("Geometric(0.2)",     m_geometric,   4.0,     20.0,    &rng);
    check_moments("Poisson(4)",         m_poisson,     4.0,     4.0,     &rng);
    check_moments("Poisson(100)",       m_poisson_big, 100.0,   100.0,   &rng);
    check_moments("Uniform(-1,3)",      m_uniform,     1.0,     4.0/3.0, &rng);
    check_moments("Exponential(0.5)",   m_exponential, 2.0,     4.0,     &rng);
    check_moments("Normal(10,2)",       m_normal,      10.0,    4.0,     &rng);
    check_moments("Gamma(3,2)",         m_gamma,       6.0,     12.0,    &rng);
    check_moments("Gamma(0.5,1)",       m_gamma_small, 0.5,     0.5,     &rng);
    check_moments("Beta(2,5)",          m_beta,        2.0/7.0, 10.0/392.0, &rng);
    check_moments("ChiSquared(6)",      m_chi_squared, 6.0,     12.0,    &rng);
    check_moments("StudentT(8)",        m_student_t,   0.0,     8.0/6.0, &rng);
    /* Var[F(d1,d2)] = 2 d2^2 (d1+d2-2) / (d1 (d2-2)^2 (d2-4)) = 1.08 here. */
    check_moments("F(5,12)",            m_f,           1.2,     1.08,    &rng);
}

/* ---- 7. Simulated type-I error rate of the t-test (paper §6.3) ---- */
static void test_type_one_error_rate(void)
{
    /* Under H0, p-values are Uniform(0,1), so about 5% of tests should
     * reject at alpha = 0.05.  Binomial(2000, 0.05) has sd ~ 9.7, so
     * [60, 140] is a ~4-sigma acceptance band. */
    ps_rng rng;
    ps_rng_seed(&rng, 1234);
    enum { TRIALS = 2000, M = 25 };
    double x[M];
    int rejections = 0;
    for (int t = 0; t < TRIALS; t++) {
        for (int i = 0; i < M; i++)
            x[i] = ps_normal(&rng, 3.0, 1.5);
        ps_t_test r = ps_t_test_one_sample(x, M, 3.0);
        if (r.ok && r.p < 0.05)
            rejections++;
    }
    CHECK(rejections >= 60 && rejections <= 140);
}

int main(void)
{
    test_rng();
    test_pmf_pdf_fixtures();
    test_cdf_fixtures();
    test_special_identities();
    test_t_tests();
    test_sampler_moments();
    test_type_one_error_rate();

    printf("%d checks, %d failure(s)\n", checks, failures);
    return failures == 0 ? 0 : 1;
}
