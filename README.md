# probsim — classical probability distributions in C99

*Last updated 2026-08-09.*

A small, portable C99 library and driver application that simulates the major
classical probability distributions, together with a typeset companion paper
(`paper/probsim.tex` → `paper/probsim.pdf`) explaining what each distribution is
for in statistical inference and modeling, and how the whole family funnels
into the *t*-test.

## Layout

| Path                    | Contents                                                        | Documented in |
|-------------------------|-----------------------------------------------------------------|---------------|
| `include/probsim.h`     | Public API: RNG, samplers, densities/CDFs, stats, *t*-tests     | paper §2–§5   |
| `src/rng.c`             | xoshiro256++ generator + splitmix64 seeding                     | paper §2      |
| `src/special.c`         | Regularized incomplete beta/gamma functions                     | paper §5.1    |
| `src/dist.c`            | Samplers, PDFs/PMFs and CDFs for the 16 distributions           | paper §3–§4   |
| `src/stats.c`           | Descriptive statistics (Welford), *t*-tests, Tarone's test      | paper §5      |
| `app/simulate.c`        | Driver: moments, balls-into-bins, *t*-test and Tarone demos     | paper §6      |
| `tests/test_probsim.c`  | Unit tests: CDF spot values, moment checks, *t*-test fixtures   | paper §6      |
| `paper/probsim.tex`     | Companion paper (LaTeX)                                         | —             |
| `paper/references.bib`  | Bibliography (BibTeX)                                           | —             |
| `docs/index.html`       | Interactive web version of the paper (GitHub Pages)             | —             |
| `docs/probsim.pdf`      | Built copy of the paper, served alongside the site              | —             |
| `bibsrc/`               | Provenance map of the local reference PDFs (see `bibsrc/README.md`) | —         |
| `scripts/fetch_refs.sh` | Downloads the open-access reference PDFs into `../_refs/distribs/` | —            |
| `scripts/charpoly.py`   | Derives & certifies the RNG's primitive characteristic polynomial | paper §2.2  |

Each source file carries comments pointing back to the section of the paper
that explains it, and the paper cites the implementing function for every
distribution it discusses, so you can read the two side by side.

## Build, test, run (Ubuntu)

```sh
make            # builds libprobsim.a and the `simulate` driver
make test       # builds and runs the unit tests
make run        # runs the driver (simulation report + t-test demos)
make docs       # builds paper/probsim.pdf (needs texlive + latexmk)
make fetch-refs # downloads open-access reference PDFs into ../_refs/distribs/
make clean
```

Requirements: any C99 compiler (`gcc`/`clang`), `make`; for the paper,
`texlive-latex-base`, `texlive-latex-recommended` and `latexmk`.

The library uses only the C99 standard library (`math.h` incl. `lgamma`,
`erf`); no external dependencies.

## Web version (GitHub Pages)

`docs/index.html` is a self-contained interactive edition of the paper —
same sections and worked examples, but every figure is animated (family
tree, inverse-transform sampling, polar Box–Muller, density explorers,
and a live *t*-test playground).  It needs no build step and no network
resources.  To publish it:

1. Push this repository to GitHub.
2. In the repo: *Settings → Pages → Deploy from a branch → `main` /
   `docs`*.
3. The `REPO_URL` constant in `docs/index.html` and `\siteurl` in
   `paper/probsim.tex` are set for `sanya-shopper/distribs`
   (site: <https://sanya-shopper.github.io/distribs/>); update both if
   the repository moves.

The site links to `docs/probsim.pdf` (a committed copy of the built
paper) so the full narrative is always one click away.  This copy is kept
in sync mechanically: `make docs` rebuilds the paper **and** refreshes
`docs/probsim.pdf`, `make check-sync` verifies the two are identical (and
that the PDF is not older than its LaTeX sources), and `make
install-hooks` installs a git pre-commit hook (`scripts/pre-commit`) that
refuses any commit touching `paper/*.tex` or `*.bib` with a stale or
unstaged site PDF.  Run `make install-hooks` once per clone.  The paper
links back to the site from each figure it animates.

## Distributions covered

Discrete: Bernoulli, binomial, geometric, Poisson, beta-binomial,
negative binomial.  Continuous: uniform, exponential, normal, gamma,
beta, chi-squared, Student's *t*, Fisher's *F*, Rayleigh, Gumbel.

The beta-binomial (paper §3.6) is the binomial with its *p* drawn afresh
from a beta for each unit rather than fixed once for all — the standard
model for counts that cluster.  It comes with Tarone's C(α) test for
binomial overdispersion (`ps_tarone_z()`, paper §5.4), which asks of real
data whether one coin or many were at work.

The chi-squared works both sides of the discrete/continuous divide
(paper §4.6): Pearson's balls-into-bins Q is a discrete statistic whose
continuous limit is the χ² — the multinomial's de Moivre–Laplace — while
every even-df χ² tail is a finite Poisson sum, which is how
`ps_poisson_cdf()` is evaluated (one incomplete-gamma call) and what the
driver prints, to machine precision, next to the continuous tail areas.

The negative binomial, Rayleigh and Gumbel form the "hash modeler's
annex" (paper §3.5, §4.9–§4.10): less common in introductory courses, but
exactly the laws that govern an ideal cryptographic hash — multi-collision
search cost (negative binomial), the birthday bound (Rayleigh), and
worst-case bucket loads and longest runs (Gumbel).

The negative binomial also carries a worked example of what using the
wrong count law costs (paper §3.5): the attacker's block count in a
Bitcoin double-spend race is exactly NB(*z*, *p*), and the white paper's
§11 substitutes the mean elapsed time to get a Poisson with the right
mean but 1/*p* times too little variance — understating the merchant's
risk by a factor that grows with every confirmation waited for (6× at ten
confirmations against a 10% attacker, 84× at twenty).  The exact answer is
one more incomplete-beta call, *I*<sub>4pq</sub>(*z*, ½), which
`ps_incbeta()` evaluates to full precision where the textbook series has
cancelled itself into noise.

## Sibling project: randtests

probsim has a sibling project, **randtests**
([repo](https://github.com/sanya-shopper/randtests) ·
[site](https://sanya-shopper.github.io/randtests/)), which surveys the
statistical test suites for PRNGs and hash functions — Diehard, Dieharder,
NIST SP 800-22, TestU01, PractRand, ENT, SMHasher — with the same
conventions (site in `docs/`, cross-referenced sources, fetchable refs,
sync-enforcing hooks).  The two are complementary: where probsim derives
the distributions (the uniform null of §4.1, the chi-squared of §4.6, the
Rayleigh birthday bound of §4.9, the Gumbel longest-run law of §4.10),
randtests covers the batteries that use them against real generators —
including the TestU01 and PractRand runs cited for xoshiro256++ in
paper §2.2.

## Sibling experiment: shastats

A second sibling, **shastats** (a local repository with no public remote
yet), hashed 8.6 billion 256-bit messages to test whether a prefixed byte
moves SHA-256's leading-zero distribution, and whether the inputs sharing
a leading-zero count have any common structure (both answers: no).  Its
measurements appear throughout the paper as *field applications* of the
distributions derived here: the binomial (avalanche, §3.2), the geometric
(the leading-zero law itself, §3.3), the uniform (256 p-values on trial,
§4.1), the chi-squared (one homogeneity test on 5100 degrees of freedom,
§4.6), the Gumbel (per-population champion digests, §4.10), and the
dispersion machinery (an *under*-dispersion curiosity, §5.4).  A PDF
snapshot of its README and working notes lands in
`../_refs/distribs/shastats-2026.pdf` via `make fetch-refs` when the
sibling checkout (directory `../256-shastats`) sits beside this
repository.
