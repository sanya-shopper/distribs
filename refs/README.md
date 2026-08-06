# refs/ — local copies of the bibliography

*Last updated 2026-08-06.*

This directory holds local PDF copies of the open-access references cited in
`paper/references.bib`.  Populate it with:

```sh
sh scripts/fetch_refs.sh     # or: make fetch-refs
```

(The PDFs are not committed: the cloud sandbox in which this repository was
authored blocks arbitrary outbound downloads, and re-fetching from the
canonical sources keeps the repo small and the provenance clean.)

## What lands here, and what cannot

| File                               | Reference (bib key)            | Status |
|------------------------------------|--------------------------------|--------|
| `student-1908.pdf`                 | Student 1908 (`student1908`)   | public-domain scan, U. York history-of-statistics archive |
| `box-muller-1958.pdf`              | Box & Muller 1958 (`boxmuller1958`) | open access, Project Euclid |
| `marsaglia-tsang-ziggurat-2000.pdf`| Marsaglia & Tsang 2000 (`marsaglia2000ziggurat`) | open access, Journal of Statistical Software |
| `blackman-vigna-xoshiro.pdf`       | Blackman & Vigna 2021 (`blackman2021xoshiro`) | open access, arXiv:1805.01407 |
| `oneill-pcg-2014.pdf`              | O'Neill 2014 (`oneill2014pcg`) | free technical report, Harvey Mudd College |
| `marsaglia-2003-xorshift.pdf`      | Marsaglia 2003 (`marsaglia2003xorshift`) | open access, Journal of Statistical Software |
| `bayes-1763.pdf`                   | Bayes 1763 (`bayes1763`)       | open access, Royal Society |
| `wilks-1938.pdf`                   | Wilks 1938 (`wilks1938`)       | open access, Project Euclid |
| `nist-sp800-22r1a.pdf`             | NIST SP 800-22 (`nist2010`)    | open access, NIST |
| `nakamoto-bitcoin-2008.pdf`        | Nakamoto 2008 (`nakamoto2008`) | freely published, bitcoin.org |
| `rosenfeld-2014-doublespend.pdf`   | Rosenfeld 2014 (`rosenfeld2014`) | open access, arXiv:1402.2009 |
| `grunspan-perez-marco-2018.pdf`    | Grunspan & Pérez-Marco 2018 (`grunspan2018`) | open access, arXiv:1702.02867 |
| —                                  | L'Ecuyer & Simard 2007 (`lecuyer2007testu01`) | paywalled (ACM TOMS); cited by DOI |
| —                                  | Vigna 2016 (`vigna2016xorshift`) | paywalled (ACM TOMS); preprint arXiv:1402.6246 |
| —                                  | PractRand (`dotyhumphrey_practrand`) | software; cited by URL |
| —                                  | Welch 1947 (`welch1947`)       | paywalled (Biometrika); cited by DOI |
| —                                  | Marsaglia & Tsang 2000 gamma (`marsaglia2000gamma`) | paywalled (ACM TOMS); cited by DOI |
| —                                  | Welford 1962 (`welford1962`)   | paywalled (Technometrics); cited by DOI |
| —                                  | Devroye 1986 (`devroye1986`)   | whole book free chapter-by-chapter at <http://luc.devroye.org/rnbookindex.html> |
| —                                  | Skellam 1948 (`skellam1948`)   | paywalled (JRSS-B); cited by DOI |
| —                                  | Tarone 1979 (`tarone1979`)     | paywalled (Biometrika); cited by DOI |
| —                                  | Edwards 1958 (`edwards1958`)   | paywalled (Annals of Human Genetics); cited by DOI |
| —                                  | Lindsey & Altham 1998 (`lindsey1998`) | paywalled (JRSS-C); cited by DOI |
| —                                  | Jacobsen et al. 1999 (`jacobsen1999`) | Human Reproduction; cited by DOI |
