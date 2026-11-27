#!/usr/bin/env bash
set -euo pipefail
cd -- "$(dirname -- "$0")"
if [[ -n ${TEX_ENV:-} ]]; then
    source "$TEX_ENV"
fi
export SOURCE_DATE_EPOCH="${SOURCE_DATE_EPOCH:-1795737600}"
export FORCE_SOURCE_DATE=1
tmp="$(mktemp -d)"
trap 'rm -rf "$tmp"' EXIT
for pass in 1 2 3; do
    if ! pdflatex -output-format=pdf -halt-on-error -interaction=nonstopmode \
           -output-directory="$tmp" main.tex >"$tmp/pass-$pass.txt"; then
        tail -80 "$tmp/pass-$pass.txt" >&2
        exit 1
    fi
done
if grep -E '(^!|undefined|Overfull|Token not allowed)' "$tmp/main.log"; then
    echo "LaTeX diagnostics need review" >&2
    exit 1
fi
cp "$tmp/main.pdf" manuscript.pdf
pdfinfo manuscript.pdf | grep -E '^(Title|Author|Pages|File size):'
