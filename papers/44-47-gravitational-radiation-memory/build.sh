#!/bin/sh
set -eu
cd "$(dirname "$0")"
out="$(mktemp -d)"
trap 'rm -rf "$out"' EXIT HUP INT TERM
cp manuscript.tex material.tex recovery.tex rank.tex inverse-square.tex "$out/"
export SOURCE_DATE_EPOCH=1795737600 FORCE_SOURCE_DATE=1 TZ=UTC
cd "$out"
if [ -n "${AREA_PDF_MAP:-}" ]; then
  entry="\\pdfoutput=1\\pdfmapfile{$AREA_PDF_MAP}\\input{manuscript.tex}"
else
  entry='\pdfoutput=1\input{manuscript.tex}'
fi
for pass in 1 2 3; do
  if ! pdflatex -jobname=manuscript -halt-on-error -interaction=nonstopmode "$entry" > "pass-$pass.log"; then
    cat "pass-$pass.log" >&2
    exit 1
  fi
done
if grep -E '^(!|Overfull|Underfull)|Warning:' manuscript.log; then
  echo 'Unexpected TeX diagnostics on final pass.' >&2
  exit 1
fi
cd "$OLDPWD"
cp "$out/manuscript.pdf" manuscript.pdf
