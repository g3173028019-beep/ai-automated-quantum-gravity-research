#!/bin/sh
set -eu
cd "$(dirname "$0")"
job=main
build_dir="$(mktemp -d)"
trap 'rm -rf "$build_dir"' EXIT HUP INT TERM
export SOURCE_DATE_EPOCH=1795737600
export FORCE_SOURCE_DATE=1
for pass in 1 2; do
    if ! pdflatex -jobname="$job" -output-directory="$build_dir" -interaction=nonstopmode -halt-on-error '\pdfoutput=1 \input{main.tex}' >"$build_dir/run.log" 2>&1; then
        cat "$build_dir/run.log" >&2
        exit 1
    fi
done
if grep -E 'undefined|Overfull|Underfull|Citation.*Warning|Rerun to' "$build_dir/run.log"; then
    exit 1
fi
cp "$build_dir/$job.pdf" manuscript.pdf
