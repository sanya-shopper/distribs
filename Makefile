# probsim — build the library, driver, tests, and the companion paper.
# See README.md for an overview and paper/probsim.tex for the paper itself.

CC      ?= gcc
CFLAGS  ?= -std=c99 -Wall -Wextra -pedantic -O2
CPPFLAGS = -Iinclude
LDLIBS   = -lm

LIB_SRCS = src/rng.c src/special.c src/dist.c src/stats.c
LIB_OBJS = $(LIB_SRCS:.c=.o)
LIB      = libprobsim.a

.PHONY: all test run docs check-sync install-hooks fetch-refs clean

all: simulate test_probsim

$(LIB): $(LIB_OBJS)
	$(AR) rcs $@ $^

%.o: %.c include/probsim.h
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

simulate: app/simulate.o $(LIB)
	$(CC) $(CFLAGS) -o $@ $^ $(LDLIBS)

test_probsim: tests/test_probsim.o $(LIB)
	$(CC) $(CFLAGS) -o $@ $^ $(LDLIBS)

test: test_probsim
	./test_probsim

run: simulate
	./simulate

# The paper: LaTeX + BibTeX, driven by latexmk.  Bibliography sources live in
# paper/references.bib; local PDFs of open-access references live in refs/.
# The built PDF is ALWAYS copied to docs/probsim.pdf so the GitHub Pages
# site serves the current paper — never edit docs/probsim.pdf by hand.
docs:
	cd paper && latexmk -pdf -bibtex -interaction=nonstopmode probsim.tex
	cp paper/probsim.pdf docs/probsim.pdf

# Rebuild-if-needed (latexmk is checksum-based, so this is fast when the
# sources are unchanged), re-sync the served copy, then fail if the PDFs
# now differ from what is staged — i.e. if a content change would land
# without its rebuilt paper.  The pre-commit hook runs this; see
# install-hooks.
check-sync:
	@$(MAKE) -s docs >/dev/null 2>&1
	@if git rev-parse --git-dir >/dev/null 2>&1 \
	  && ! git diff --quiet -- paper/probsim.pdf docs/probsim.pdf; then \
	  echo "OUT OF SYNC: the rebuilt paper differs from the staged PDFs —"; \
	  echo "run: git add paper/probsim.pdf docs/probsim.pdf"; exit 1; fi
	@echo "paper and site PDF are in sync"

# One-time per clone: install the pre-commit hook that enforces check-sync.
install-hooks:
	cp scripts/pre-commit .git/hooks/pre-commit
	chmod +x .git/hooks/pre-commit
	@echo "pre-commit hook installed"

# Local PDFs of the open-access references (see refs/README.md).
fetch-refs:
	sh scripts/fetch_refs.sh

clean:
	rm -f $(LIB_OBJS) app/simulate.o tests/test_probsim.o $(LIB) simulate test_probsim
	cd paper && latexmk -C probsim.tex 2>/dev/null || true
