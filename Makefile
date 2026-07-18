# probsim — build the library, driver, tests, and the companion paper.
# See README.md for an overview and paper/probsim.tex for the paper itself.

CC      ?= gcc
CFLAGS  ?= -std=c99 -Wall -Wextra -pedantic -O2
CPPFLAGS = -Iinclude
LDLIBS   = -lm

LIB_SRCS = src/rng.c src/special.c src/dist.c src/stats.c
LIB_OBJS = $(LIB_SRCS:.c=.o)
LIB      = libprobsim.a

.PHONY: all test run docs fetch-refs clean

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
docs:
	cd paper && latexmk -pdf -bibtex -interaction=nonstopmode probsim.tex

# Local PDFs of the open-access references (see refs/README.md).
fetch-refs:
	sh scripts/fetch_refs.sh

clean:
	rm -f $(LIB_OBJS) app/simulate.o tests/test_probsim.o $(LIB) simulate test_probsim
	cd paper && latexmk -C probsim.tex 2>/dev/null || true
