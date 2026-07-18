/*
 * rng.c — pseudo-random number generation for probsim.
 *
 * Explained in the companion paper, paper/probsim.tex §2 ("Where the
 * randomness comes from"); §2.2 examines this generator in detail —
 * what is proved (period 2^256-1, exact jump functions, the forbidden
 * all-zero state), what is empirical (BigCrush, PractRand, the
 * Hamming-weight dependency test), and what remains open.  The engine's
 * primitive characteristic polynomial, on which the period theorem
 * rests, is derived from this very code and certified by
 * scripts/charpoly.py (printed in paper §2.2).
 * Public interface: include/probsim.h.
 *
 * The generator is xoshiro256++ (Blackman & Vigna 2021; local copy at
 * refs/blackman-vigna-xoshiro.pdf).  It is small, fast, portable, and
 * shows no systematic failures on BigCrush — more than adequate for
 * statistical simulation, though NOT for cryptography (see paper §2.2).
 *
 * Seeding uses splitmix64, as the xoshiro authors recommend, so that
 * even trivially-different seeds (0, 1, 2, ...) yield well-separated
 * initial states.
 */
#include "probsim.h"

/* splitmix64 step: a 64-bit finalizer-style generator used only for seeding. */
static uint64_t splitmix64(uint64_t *state)
{
    uint64_t z = (*state += 0x9E3779B97F4A7C15ULL);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}

void ps_rng_seed(ps_rng *rng, uint64_t seed)
{
    for (int i = 0; i < 4; i++)
        rng->s[i] = splitmix64(&seed);
}

static uint64_t rotl(uint64_t x, int k)
{
    return (x << k) | (x >> (64 - k));
}

uint64_t ps_rng_next(ps_rng *rng)
{
    uint64_t *s = rng->s;
    const uint64_t result = rotl(s[0] + s[3], 23) + s[0];
    const uint64_t t = s[1] << 17;

    s[2] ^= s[0];
    s[3] ^= s[1];
    s[1] ^= s[2];
    s[0] ^= s[3];
    s[2] ^= t;
    s[3] = rotl(s[3], 45);

    return result;
}

/*
 * Uniform double on [0,1): take the top 53 bits (a double's full mantissa
 * precision) and scale by 2^-53.  Every sampler in src/dist.c reduces to
 * calls of this function — see paper §2.1 on the "uniform is enough" idea.
 */
double ps_runif(ps_rng *rng)
{
    return (double)(ps_rng_next(rng) >> 11) * 0x1.0p-53;
}
