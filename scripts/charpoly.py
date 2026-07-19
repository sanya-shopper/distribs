#!/usr/bin/env python3
"""
charpoly.py -- derive and certify the characteristic polynomial of the
xoshiro256 linear engine used by probsim (src/rng.c).

Referenced by the companion paper, paper/probsim.tex §2.2 ("What is
proved"), which prints the polynomial this script outputs.  The paper's
period theorem rests on this polynomial being PRIMITIVE over GF(2); this
script re-establishes that from the code itself, with no trust in any
external table:

  1. Simulate the engine's linear state transition (the xoshiro256 update
     WITHOUT the ++ scrambler -- the scrambler does not touch the state)
     bit-exactly, and record one state bit per step.
  2. Berlekamp--Massey over GF(2) recovers the minimal recurrence of that
     bit sequence; for an engine with irreducible characteristic
     polynomial the minimal polynomial IS the characteristic polynomial
     P(x), of degree 256.
  3. Certify the recurrence against thousands of further steps and
     against a second, independent state bit.
  4. Certify primitivity: x^(2^256-1) = 1 (mod P) and
     x^((2^256-1)/q) != 1 (mod P) for every prime q | 2^256-1
     (the complete factorization of 2^256-1 is classical; the script
     asserts the product reconstitutes 2^256-1).  Order exactly 2^256-1
     forces P primitive, which is precisely the single-full-cycle
     theorem of paper §2.2.

Usage:  python3 scripts/charpoly.py        (runs in under a second)

The answer, for the reader who does not run it (and asserted against the
derivation below, so this docstring cannot silently rot): P(x) is a
degree-256, weight-115 polynomial over GF(2) -- an odd weight, as every
irreducible binary polynomial must have -- whose coefficients, packed
with bit i holding the coefficient of x^i, are

    0x1_0003c03c_3f3ecb19_04b4edcf_26259f85
      _0280002b_cefd1a5e_9d116f2b_b0f0f001

i.e. P(x) = x^256 + x^241 + x^240 + ... + x^12 + ... + 1 (115 terms; the
leading x^256 is the top hex digit's 1, the trailing +1 the bottom
digit's).  This is the polynomial printed in paper §2.2.

Theory background: Lidl & Niederreiter, "Finite Fields"; the engine and
its design are from Blackman & Vigna (2021) -- see paper/references.bib.
"""

# The documented answer (see docstring).  main() derives P independently
# and asserts equality, so a change to the engine in src/rng.c that is
# not mirrored here -- or vice versa -- fails loudly.
EXPECTED_P = int(
    "1" "0003c03c" "3f3ecb19" "04b4edcf" "26259f85"
        "0280002b" "cefd1a5e" "9d116f2b" "b0f0f001", 16)

M64 = (1 << 64) - 1


def rotl(x: int, k: int) -> int:
    return ((x << k) | (x >> (64 - k))) & M64


def engine_step(s: list) -> None:
    """One step of the xoshiro256 LINEAR engine (mirrors src/rng.c
    ps_rng_next() with the nonlinear output scrambler omitted)."""
    t = (s[1] << 17) & M64
    s[2] ^= s[0]
    s[3] ^= s[1]
    s[1] ^= s[2]
    s[0] ^= s[3]
    s[2] ^= t
    s[3] = rotl(s[3], 45)


def bit_sequence(bit_word: int, bit_index: int, n: int) -> list:
    """n successive values of one state bit, from a fixed nonzero seed."""
    s = [1, 0, 0, 0]
    out = []
    for _ in range(n):
        out.append((s[bit_word] >> bit_index) & 1)
        engine_step(s)
    return out


def berlekamp_massey(seq: list):
    """Minimal LFSR (connection polynomial C, length L) over GF(2).
    Polynomials are ints: bit i holds the coefficient of x^i; c_0 = 1 and
    sum_{i=0..L} c_i * seq[n-i] = 0 for all valid n."""
    C, B, L, m = 1, 1, 0, 1
    for n in range(len(seq)):
        d = 0
        for i in range(L + 1):
            if (C >> i) & 1:
                d ^= seq[n - i]
        if d:
            T = C
            C ^= B << m
            if 2 * L <= n:
                L, B, m = n + 1 - L, T, 1
            else:
                m += 1
        else:
            m += 1
    return C, L


# ---- GF(2)[x] arithmetic on int-encoded polynomials ----

def pmulmod(a: int, b: int, P: int, deg: int) -> int:
    """(a * b) mod P, carry-less."""
    r = 0
    while b:
        if b & 1:
            r ^= a
        a <<= 1
        b >>= 1
    while r.bit_length() > deg:
        r ^= P << (r.bit_length() - 1 - deg)
    return r


def ppowmod(e: int, P: int, deg: int) -> int:
    """x^e mod P by square-and-multiply."""
    result, base = 1, 2          # 2 encodes the polynomial x
    while e:
        if e & 1:
            result = pmulmod(result, base, P, deg)
        base = pmulmod(base, base, P, deg)
        e >>= 1
    return result


def main() -> None:
    DEG = 256
    N = 4 * DEG                  # plenty for BM (needs 2*DEG) + certification

    # 1-2. Recover the recurrence.
    seq = bit_sequence(0, 0, N)
    C, L = berlekamp_massey(seq)
    assert L == DEG, f"expected degree {DEG}, Berlekamp-Massey found {L}"
    # Characteristic-polynomial form: P(x) = x^L * C(1/x)  (bit reversal).
    P = 0
    for i in range(L + 1):
        if (C >> i) & 1:
            P |= 1 << (L - i)

    # 3. Certify the recurrence on the tail of this sequence...
    for n in range(DEG, N):
        acc = 0
        for i in range(DEG):
            if (P >> i) & 1:
                acc ^= seq[n - DEG + i]
        assert acc == seq[n], f"recurrence fails at step {n}"
    # ...and on an independent state bit (word 2, bit 63).
    seq2 = bit_sequence(2, 63, N)
    for n in range(DEG, N):
        acc = 0
        for i in range(DEG):
            if (P >> i) & 1:
                acc ^= seq2[n - DEG + i]
        assert acc == seq2[n], f"recurrence fails on bit 2:63 at step {n}"

    # 4. Certify primitivity.  Complete prime factorization of 2^256 - 1.
    order = (1 << 256) - 1
    primes = [3, 5, 17, 257, 641, 65537, 274177, 6700417,
              67280421310721, 59649589127497217, 5704689200685129054721]
    prod = 1
    for q in primes:
        prod *= q
    assert prod == order, "factor list does not reconstitute 2^256 - 1"

    assert ppowmod(order, P, DEG) == 1, "x^(2^256-1) != 1 (mod P)"
    for q in primes:
        assert ppowmod(order // q, P, DEG) != 1, \
            f"order divides (2^256-1)/{q}: P is not primitive"

    # 5. The derivation must agree with the documented answer up top.
    assert P == EXPECTED_P, \
        "derived polynomial differs from EXPECTED_P in this file's " \
        "docstring — engine mirror and documentation are out of sync"

    weight = bin(P).count("1")
    print("characteristic polynomial of the xoshiro256 engine over GF(2):")
    print(f"  degree : {DEG}")
    print(f"  weight : {weight} nonzero coefficients")
    print(f"  P (hex, coefficient of x^i at bit i):")
    h = format(P, "x")
    for i in range(0, len(h), 32):
        print(f"    {h[i:i+32]}")
    print("certified: recurrence holds on 2 independent state bits over "
          f"{N} steps;")
    print("certified: ord(x) = 2^256 - 1 exactly => P is primitive "
          "(single full cycle).")


if __name__ == "__main__":
    main()
