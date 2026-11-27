#!/usr/bin/env bash
set -euo pipefail
cd -- "$(dirname -- "${BASH_SOURCE[0]}")"
export SOURCE_DATE_EPOCH="${SOURCE_DATE_EPOCH:-1795737600}"
export FORCE_SOURCE_DATE=1
tmp="$(mktemp -d)"
trap 'rm -rf "$tmp"' EXIT
for pass in 1 2 3; do
  pdflatex -interaction=nonstopmode -halt-on-error -output-directory="$tmp" main.tex > "$tmp/pass.log" || { cat "$tmp/pass.log"; exit 1; }
done
if grep -Ei 'undefined references|Citation .+ undefined|Reference .+ undefined|There were multiply-defined' "$tmp/main.log" ; then exit 1; fi
cp "$tmp/main.pdf" manuscript.pdf
