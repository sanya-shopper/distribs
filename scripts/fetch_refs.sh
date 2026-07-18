#!/bin/sh
# fetch_refs.sh — download local PDF copies of the open-access references
# cited in paper/references.bib into refs/.
#
# Run from the repository root:   sh scripts/fetch_refs.sh
#
# Only references with a legitimately free PDF are fetched; paywalled items
# (Welch 1947, Marsaglia & Tsang's TOMS gamma paper, Welford 1962) are cited
# by DOI in the bibliography instead — see refs/README.md for the full map.
#
# Note: the cloud sandbox this project was authored in blocks arbitrary
# outbound downloads, which is why the PDFs are fetched by script rather
# than committed.  On an ordinary Ubuntu machine this just works.

set -eu
cd "$(dirname "$0")/.."
mkdir -p refs

fetch() {
    out="refs/$1"
    url="$2"
    if [ -s "$out" ]; then
        echo "have    $out"
        return
    fi
    echo "fetch   $out"
    curl -fsSL --retry 2 -o "$out" "$url" || echo "FAILED  $out  ($url)"
}

# Student (1908), "The probable error of a mean" — public-domain scan
# hosted by the University of York's history-of-statistics archive.
fetch student-1908.pdf \
    "https://www.york.ac.uk/depts/maths/histstat/student.pdf"

# Box & Muller (1958) — open access at Project Euclid.
fetch box-muller-1958.pdf \
    "https://projecteuclid.org/journals/annals-of-mathematical-statistics/volume-29/issue-2/A-Note-on-the-Generation-of-Random-Normal-Deviates/10.1214/aoms/1177706645.pdf"

# Marsaglia & Tsang (2000), "The Ziggurat Method" — open access, JSS.
fetch marsaglia-tsang-ziggurat-2000.pdf \
    "https://www.jstatsoft.org/index.php/jss/article/view/v005i08/ziggurat.pdf"

# Blackman & Vigna (2019), xoshiro/xoroshiro generators — arXiv.
fetch blackman-vigna-xoshiro.pdf \
    "https://arxiv.org/pdf/1805.01407"

# O'Neill (2014), PCG technical report — free from the author/HMC.
fetch oneill-pcg-2014.pdf \
    "https://www.cs.hmc.edu/tr/hmc-cs-2014-0905.pdf"

# Marsaglia (2003), "Xorshift RNGs" — open access, JSS (cited in the
# paper's §2.2 on the ancestry of xoshiro256++).
fetch marsaglia-2003-xorshift.pdf \
    "https://www.jstatsoft.org/index.php/jss/article/view/v008i14/xorshift.pdf"

# Bayes (1763) — open access at the Royal Society (paper §4.5 history).
fetch bayes-1763.pdf \
    "https://royalsocietypublishing.org/doi/pdf/10.1098/rstl.1763.0053"

# Wilks (1938) — open access at Project Euclid (paper §4.6 history).
fetch wilks-1938.pdf \
    "https://projecteuclid.org/journals/annals-of-mathematical-statistics/volume-9/issue-1/The-Large-Sample-Distribution-of-the-Likelihood-Ratio-for-Testing/10.1214/aoms/1177732360.pdf"

# NIST SP 800-22 rev 1a — open access (paper §4.10, randomness testing).
fetch nist-sp800-22r1a.pdf \
    "https://nvlpubs.nist.gov/nistpubs/legacy/sp/nistspecialpublication800-22r1a.pdf"

echo "done."
