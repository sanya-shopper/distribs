/*
 * stats.c — descriptive statistics and the t-tests that tie the
 * distribution family together.
 *
 * Explained in the companion paper, docs/probsim.tex §5 ("From
 * distributions to inference: the t-test").  Public interface:
 * include/probsim.h.
 *
 * The two tests here are the payoff of the whole library: the normal
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
