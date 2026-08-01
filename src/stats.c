/*
 * stats.c — descriptive statistics and the t-tests that tie the
 * distribution family together.
 *
 * Explained in the companion paper, paper/probsim.tex §5 ("From
 * distributions to inference: the t-test").  Public interface:
 * include/probsim.h.
 *
 * Also here: Tarone's C(alpha) test for binomial overdispersion (paper
 * §5.4), which is how the binomial and the beta-binomial of §3.6 are
 * told apart in practice.
 *
 * The two t-tests here are the payoff of the whole library: the normal
 * (src/dist.c: ps_normal) models the data, the chi-squared distribution
 * describes the sample variance, their ratio is Student's t (Student
 * 1908, refs/student-1908.pdf), and ps_student_t_cdf turns the observed
 * statistic into a p-value.  Welch's unequal-variance test (Welch 1947)
 * follows the same route with an approximate df.
 */
#include "probsim.h"

#include <math.h>

/*
 * One-pass mean and variance by Welford's method — numerically stable
 * where the naive sum-of-squares formula cancels catastrophically.
 * (Paper §5.2 explains why the n-1 divisor is what makes the t statistic
 * work out.)
 */
ps_summary ps_summarize(const double *x, size_t n)
{
    ps_summary s = { n, NAN, NAN, NAN, NAN, NAN };
    if (x == NULL || n == 0)
        return s;

    double mean = 0.0, m2 = 0.0, min = x[0], max = x[0];
    for (size_t i = 0; i < n; i++) {
        double delta = x[i] - mean;
        mean += delta / (double)(i + 1);
        m2 += delta * (x[i] - mean);
        if (x[i] < min) min = x[i];
        if (x[i] > max) max = x[i];
    }
    s.mean = mean;
    s.min = min;
    s.max = max;
    if (n > 1) {
        s.var = m2 / (double)(n - 1);
        s.sd = sqrt(s.var);
    }
    return s;
}

/*
 * Tarone's C(alpha) test for binomial overdispersion — paper §5.4.
 *
 * Given m clusters (k[i] successes in n[i] trials), H0 says every trial
 * everywhere is a flip of the SAME coin p; the beta-binomial alternative
 * (src/dist.c: ps_beta_binomial, paper §3.6) says each cluster gets its
 * own coin.  Write phat for the pooled proportion and
 *
 *     S = sum_i (k[i] - n[i] phat)^2 / (phat (1 - phat)),
 *     T = sum_i n[i] (n[i] - 1)                    ("ordered pairs").
 *
 * Under H0, E[S] = sum n[i]; under the mixture, E[S] = sum n[i] + rho T.
 * Hence the moment estimate rho = (S - N)/T and Tarone's statistic
 *
 *     Z = (S - N) / sqrt(2 T)  =  rho * sqrt(T / 2),
 *
 * standard normal under H0 and large-positive under overdispersion.  The
 * second form is the useful one for planning: detecting a correlation of
 * rho needs about T = 2 (z/rho)^2 ordered within-cluster pairs, which is
 * why tiny rho only becomes visible in registry-sized data (paper §3.6).
 * Note T = 0 when every cluster has one trial: no pairs, nothing to
 * compare, and the test is correctly undefined.
 */
ps_overdispersion ps_tarone_z(const long *k, const long *n, size_t m)
{
    ps_overdispersion r = { false, NAN, NAN, NAN, NAN };
    if (k == NULL || n == NULL || m == 0)
        return r;

    double sum_k = 0.0, sum_n = 0.0, pairs = 0.0;
    for (size_t i = 0; i < m; i++) {
        if (n[i] < 0 || k[i] < 0 || k[i] > n[i])
            return r;                      /* not a (successes, trials) pair */
        sum_k += (double)k[i];
        sum_n += (double)n[i];
        pairs += (double)n[i] * (double)(n[i] - 1);
    }
    if (!(pairs > 0.0))
        return r;  /* every cluster of size <= 1: no within-cluster pairs */

    double phat = sum_k / sum_n;
    if (!(phat > 0.0 && phat < 1.0))
        return r;  /* all successes or all failures: dispersion undefined */

    double s = 0.0;
    for (size_t i = 0; i < m; i++) {
        double d = (double)k[i] - (double)n[i] * phat;
        s += d * d;
    }
    s /= phat * (1.0 - phat);

    r.phat = phat;
    r.rho  = (s - sum_n) / pairs;
    r.z    = (s - sum_n) / sqrt(2.0 * pairs);
    r.p    = ps_normal_cdf(-r.z, 0.0, 1.0);   /* one-sided: overdispersion */
    r.ok   = true;
    return r;
}

/* Shared tail: turn (t, df) into a two-sided p-value via the t CDF. */
static double two_sided_p(double t, double df)
{
    return 2.0 * ps_student_t_cdf(-fabs(t), df);
}

/*
 * One-sample t-test of H0: E[X] = mu0 — paper §5.2.
 *
 *   t = (xbar - mu0) / (s / sqrt(n)),   df = n - 1.
 */
ps_t_test ps_t_test_one_sample(const double *x, size_t n, double mu0)
{
    ps_t_test r = { false, NAN, NAN, NAN, NAN, NAN };
    if (x == NULL || n < 2)
        return r;

    ps_summary s = ps_summarize(x, n);
    if (!(s.sd > 0.0))
        return r;  /* zero-variance sample: t undefined */

    r.estimate = s.mean;
    r.se = s.sd / sqrt((double)n);
    r.t = (s.mean - mu0) / r.se;
    r.df = (double)(n - 1);
    r.p = two_sided_p(r.t, r.df);
    r.ok = true;
    return r;
}

/*
 * Welch's two-sample t-test of H0: E[X] = E[Y] — paper §5.3.
 *
 *   t  = (xbar - ybar) / sqrt(sx^2/nx + sy^2/ny)
 *   df = (sx^2/nx + sy^2/ny)^2
 *        / ( (sx^2/nx)^2/(nx-1) + (sy^2/ny)^2/(ny-1) )     (Welch 1947)
 *
 * No equal-variance assumption, which is why it is the modern default.
 */
ps_t_test ps_t_test_welch(const double *x, size_t nx,
                          const double *y, size_t ny)
{
    ps_t_test r = { false, NAN, NAN, NAN, NAN, NAN };
    if (x == NULL || y == NULL || nx < 2 || ny < 2)
        return r;

    ps_summary sx = ps_summarize(x, nx);
    ps_summary sy = ps_summarize(y, ny);
    double vx = sx.var / (double)nx;
    double vy = sy.var / (double)ny;
    if (!(vx + vy > 0.0))
        return r;  /* both samples constant: t undefined */

    r.estimate = sx.mean - sy.mean;
    r.se = sqrt(vx + vy);
    r.t = r.estimate / r.se;
    r.df = (vx + vy) * (vx + vy)
         / (vx * vx / (double)(nx - 1) + vy * vy / (double)(ny - 1));
    r.p = two_sided_p(r.t, r.df);
    r.ok = true;
    return r;
}
