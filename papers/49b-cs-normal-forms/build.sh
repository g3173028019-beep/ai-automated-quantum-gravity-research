#!/usr/bin/env bash
set -euo pipefail
cd -- "$(dirname -- "$0")"
export SOURCE_DATE_EPOCH=1795737600
export FORCE_SOURCE_DATE=1
export TZ=UTC
export LC_ALL=C
builddir=$(mktemp -d)
trap 'rm -rf -- "$builddir"' EXIT
export TEXMFVAR="$builddir/texmf-var"
export VARTEXFONTS="$builddir/fonts"
if ! kpsewhich article.cls >/dev/null; then
  # TeX Live data may be present without registered kpathsea defaults.
  data=$(printenv TEXLIVE_DIST || printf /usr/share/texlive/texmf-dist)
  [[ -r "$data/tex/latex/base/article.cls" ]] || {
    echo "TeX Live's article.cls was not found." >&2; exit 1;
  }
  export TEXINPUTS=".:$data/tex//:"
  export TFMFONTS="$data/fonts/tfm//:"
  export T1FONTS="$data/fonts/type1//:"
  export ENCFONTS="$data/fonts/enc//:"
  cat "$data"/fonts/map/dvips/amsfonts/*.map >"$builddir/pdftex.map"
  export TEXFONTMAPS="$builddir:$data/fonts/map//:"
  export TEXFORMATS="$builddir:"
  (cd -- "$builddir"; pdftex -ini -jobname=pdflatex -halt-on-error -interaction=nonstopmode '*latex.ini') >"$builddir/format.log" || {
    cat "$builddir/format.log" >&2; exit 1;
  }
fi
for n in 1 2; do
  pdflatex -file-line-error -halt-on-error -interaction=nonstopmode -jobname=main \
    -output-directory="$builddir" '\pdfoutput=1 \input{main.tex}' \
    >"$builddir/stdout" || { cat "$builddir/stdout" >&2; exit 1; }
done
if grep -E 'undefined references|undefined citations|multiply defined|Overfull' "$builddir/main.log"; then
  echo "Reference or layout error; see pdflatex." >&2
  exit 1
fi
cp -- "$builddir/main.pdf" manuscript.pdf
