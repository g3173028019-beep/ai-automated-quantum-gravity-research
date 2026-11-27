#!/usr/bin/env bash
set -euo pipefail
cd -- "$(dirname -- "${BASH_SOURCE[0]}")"
if [[ -n "${TEX_ENV:-}" ]]; then source "$TEX_ENV"; fi
export SOURCE_DATE_EPOCH="${SOURCE_DATE_EPOCH:-1795737600}"
export FORCE_SOURCE_DATE=1
tmp="$(mktemp -d)"
trap 'rm -rf "$tmp"' EXIT
for pass in 1 2 3; do
  pdflatex -interaction=nonstopmode -halt-on-error -output-directory="$tmp" main.tex
done
cp "$tmp/main.pdf" manuscript.pdf
