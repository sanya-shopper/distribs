/*
 * special.c — the two special functions behind every classical test CDF.
 *
 * Explained in the companion paper, docs/probsim.tex §5.1 ("Two special
 * functions carry all the tables").  Public interface: include/probsim.h.
 *
 * Almost every CDF a classical statistician needs — beta, gamma,
 * chi-squared, Student's t, Fisher's F — reduces to one of:
 *
 *   ps_incbeta(x, a, b)       regularized incomplete beta  I_x(a, b)
 *   ps_incgamma_lower(a, x)   regularized lower incomplete gamma  P(a, x)
 *
 * Both are evaluated with the standard series / continued-fraction pairs
 * (Lentz's algorithm for the fractions), following the presentation in
 * Press et al., "Numerical Recipes in C", 2nd ed., §6.2 and §6.4.
 * Accuracy is ample for p-values (relative error near 1e-14 over the
 * ranges exercised in tests/test_probsim.c).
 */
#include "probsim.h"

#include <math.h>

/* Convergence controls shared by both continued fractions. */
#define PS_CF_EPS      1e-15   /* target relative error                    */
#define PS_CF_TINY     1e-300  /* Lentz floor to avoid division by zero    */
#define PS_CF_MAX_ITER 300

/* Continued-fraction core of the incomplete beta (Numerical Recipes betacf),
 * evaluated by the modified Lentz method. */
static double incbeta_cf(double x, double a, double b)
{
    double c = 1.0;
    double d = 1.0 - (a + b) * x / (a + 1.0);
    if (fabs(d) < PS_CF_TINY) d = PS_CF_TINY;
    d = 1.0 / d;
    double h = d;

    for (int m = 1; m <= PS_CF_MAX_ITER; m++) {
        int m2 = 2 * m;

        /* Even step. */
        double aa = m * (b - m) * x / ((a + m2 - 1.0) * (a + m2));
        d = 1.0 + aa * d;
        if (fabs(d) < PS_CF_TINY) d = PS_CF_TINY;
        c = 1.0 + aa / c;
        if (fabs(c) < PS_CF_TINY) c = PS_CF_TINY;
        d = 1.0 / d;
        h *= d * c;

        /* Odd step. */
        aa = -(a + m) * (a + b + m) * x / ((a + m2) * (a + m2 + 1.0));
        d = 1.0 + aa * d;
        if (fabs(d) < PS_CF_TINY) d = PS_CF_TINY;
        c = 1.0 + aa / c;
        if (fabs(c) < PS_CF_TINY) c = PS_CF_TINY;
        d = 1.0 / d;
        double del = d * c;
        h *= del;
        if (fabs(del - 1.0) < PS_CF_EPS)
            break;
    }
    return h;
}

double ps_incbeta(double x, double a, double b)
{
    if (isnan(x) || !(a > 0.0) || !(b > 0.0)) return NAN;
    if (x <= 0.0) return 0.0;
    if (x >= 1.0) return 1.0;

    /* Prefactor  x^a (1-x)^b / (a B(a,b)),  computed in log space. */
    double lbeta = lgamma(a) + lgamma(b) - lgamma(a + b);
    double front = exp(a * log(x) + b * log1p(-x) - lbeta);

    /* Use the continued fraction on whichever side converges fast,
     * exploiting the symmetry I_x(a,b) = 1 - I_{1-x}(b,a). */
    if (x < (a + 1.0) / (a + b + 2.0))
        return front * incbeta_cf(x, a, b) / a;
    return 1.0 - front * incbeta_cf(1.0 - x, b, a) / b;
}

/* Series expansion for P(a, x), good for x < a + 1. */
static double incgamma_series(double a, double x)
{
    double ap = a;
    double sum = 1.0 / a;
    double del = sum;
    for (int n = 1; n <= PS_CF_MAX_ITER; n++) {
        ap += 1.0;
        del *= x / ap;
        sum += del;
        if (fabs(del) < fabs(sum) * PS_CF_EPS)
            break;
    }
    return sum * exp(-x + a * log(x) - lgamma(a));
}

/* Continued fraction for the UPPER tail Q(a, x), good for x >= a + 1
 * (Numerical Recipes gcf), by the modified Lentz method. */
static double incgamma_cf(double a, double x)
{
    double b = x + 1.0 - a;
    double c = 1.0 / PS_CF_TINY;
    double d = 1.0 / b;
    double h = d;
    for (int i = 1; i <= PS_CF_MAX_ITER; i++) {
        double an = -i * (i - a);
        b += 2.0;
        d = an * d + b;
        if (fabs(d) < PS_CF_TINY) d = PS_CF_TINY;
        c = b + an / c;
        if (fabs(c) < PS_CF_TINY) c = PS_CF_TINY;
        d = 1.0 / d;
        double del = d * c;
        h *= del;
        if (fabs(del - 1.0) < PS_CF_EPS)
            break;
    }
    return h * exp(-x + a * log(x) - lgamma(a));
}

double ps_incgamma_lower(double a, double x)
{
    if (!(a > 0.0) || isnan(x)) return NAN;
    if (x <= 0.0) return 0.0;
    if (x < a + 1.0)
        return incgamma_series(a, x);
    return 1.0 - incgamma_cf(a, x);
}
