# probsim — classical probability distributions in C99

A small, portable C99 library and driver application that simulates the major
classical probability distributions, together with a typeset companion paper
(`docs/probsim.tex` → `docs/probsim.pdf`) explaining what each distribution is
for in statistical inference and modeling, and how the whole family funnels
into the *t*-test.

## Layout

| Path                    | Contents                                                        | Documented in |
|-------------------------|-----------------------------------------------------------------|---------------|
| `include/probsim.h`     | Public API: RNG, samplers, densities/CDFs, stats, *t*-tests     | paper §2–§5   |
| `src/rng.c`             | xoshiro256++ generator + splitmix64 seeding                     | paper §2      |
| `src/special.c`         | Regularized incomplete beta/gamma functions                     | paper §5.1    |
| `src/dist.c`            | Samplers, PDFs/PMFs and CDFs for the 12 distributions           | paper §3–§4   |
| `src/stats.c`           | Descriptive statistics (Welford), one-sample & Welch *t*-tests  | paper §5      |
| `app/simulate.c`        | Driver: theoretical vs. empirical moments, *t*-test demos       | paper §6      |
| `tests/test_probsim.c`  | Unit tests: CDF spot values, moment checks, *t*-test fixtures   | paper §6      |
| `docs/probsim.tex`      | Companion paper (LaTeX)                                         | —             |
| `docs/references.bib`   | Bibliography (BibTeX)                                           | —             |
| `refs/`                 | Local copies of open-access reference PDFs (`make fetch-refs`)  | —             |
| `scripts/fetch_refs.sh` | Downloads the open-access reference PDFs into `refs/`           | —             |

Each source file carries comments pointing back to the section of the paper
that explains it, and the paper cites the implementing function for every
distribution it discusses, so you can read the two side by side.

## Build, test, run (Ubuntu)

```sh
make            # builds libprobsim.a and the `simulate` driver
make test       # builds and runs the unit tests
make run        # runs the driver (simulation report + t-test demos)
make docs       # builds docs/probsim.pdf (needs texlive + latexmk)
make fetch-refs # downloads open-access reference PDFs into refs/
make clean
```

Requirements: any C99 compiler (`gcc`/`clang`), `make`; for the paper,
`texlive-latex-base`, `texlive-latex-recommended` and `latexmk`.

The library uses only the C99 standard library (`math.h` incl. `lgamma`,
`erf`); no external dependencies.

## Distributions covered

Discrete: Bernoulli, binomial, geometric, Poisson.
Continuous: uniform, exponential, normal, gamma, beta, chi-squared,
Student's *t*, Fisher's *F*.
