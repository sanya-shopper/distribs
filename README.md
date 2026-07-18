# probsim — classical probability distributions in C99

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
| `src/dist.c`            | Samplers, PDFs/PMFs and CDFs for the 12 distributions           | paper §3–§4   |
| `src/stats.c`           | Descriptive statistics (Welford), one-sample & Welch *t*-tests  | paper §5      |
| `app/simulate.c`        | Driver: theoretical vs. empirical moments, *t*-test demos       | paper §6      |
| `tests/test_probsim.c`  | Unit tests: CDF spot values, moment checks, *t*-test fixtures   | paper §6      |
| `paper/probsim.tex`     | Companion paper (LaTeX)                                         | —             |
| `paper/references.bib`  | Bibliography (BibTeX)                                           | —             |
| `docs/index.html`       | Interactive web version of the paper (GitHub Pages)             | —             |
| `docs/probsim.pdf`      | Built copy of the paper, served alongside the site              | —             |
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
make docs       # builds paper/probsim.pdf (needs texlive + latexmk)
make fetch-refs # downloads open-access reference PDFs into refs/
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
3. Edit the `REPO_URL` constant at the top of the `<script>` block in
   `docs/index.html` so the site's "code" links point at your repository.

The site links to `docs/probsim.pdf` (a committed copy of the built
paper) so the full narrative is always one click away; after editing the
paper, run `make docs` and copy `paper/probsim.pdf` over
`docs/probsim.pdf` to keep the two in sync.  The paper links back to the
site from each figure it animates.

## Distributions covered

Discrete: Bernoulli, binomial, geometric, Poisson.
Continuous: uniform, exponential, normal, gamma, beta, chi-squared,
Student's *t*, Fisher's *F*.
