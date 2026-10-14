#!/bin/bash
# Compile the harness against the libfac.a of a built FAC repository.
# Usage: [HARNESS_DEFS=-DNO_EESCREENMU] build_harness.sh REPO_DIR OUTPUT_EXECUTABLE
set -e
REPO="$(cd "$1" && pwd)"; OUT="$2"
HERE="$(cd "$(dirname "$0")" && pwd)"
CC_BIN="$(grep -m1 '^CC = ' "$REPO/python/Makefile" | cut -d= -f2- | xargs)"
LIBS="$(grep -m1 '^LIBS = ' "$REPO/python/Makefile" | cut -d= -f2- | sed 's/\${PYLIB}//')"
LDFLAGS="$(grep -m1 '^LDFLAGS = ' "$REPO/python/Makefile" | cut -d= -f2-)"
CFLAGS_REPO="$(grep -m1 '^CFLAGS = ' "$REPO/faclib/Makefile" | cut -d= -f2-)"
# same optimisation level as the library (taken from its Makefile), no extra workaround flags
$CC_BIN $HARNESS_DEFS -I"$REPO/faclib" -I"$REPO" $CFLAGS_REPO "$HERE/densitytosz_harness.c" -o "$OUT" -L"$REPO" $LDFLAGS $LIBS
