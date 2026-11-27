#!/bin/sh
set -eu
cd "$(dirname "$0")"
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT HUP INT TERM
export SOURCE_DATE_EPOCH=1795737600 FORCE_SOURCE_DATE=1 TZ=UTC
for pass in 1 2 3; do
  if ! pdflatex -jobname=main -output-directory="$work" -interaction=nonstopmode -halt-on-error -file-line-error main.tex >"$work/console.log" 2>&1; then
     tail -100 "$work/console.log" >&2; exit 1
  fi
done
if grep -E 'undefined|multiply defined|Overfull|Rerun to' "$work/main.log"; then exit 1; fi
cp "$work/main.pdf" manuscript.pdf
