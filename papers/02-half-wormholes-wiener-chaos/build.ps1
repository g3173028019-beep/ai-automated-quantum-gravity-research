$ErrorActionPreference = "Stop"
Push-Location $PSScriptRoot
try {
    & latexmk -pdf -bibtex -interaction=nonstopmode -halt-on-error -outdir=build main.tex
    if ($LASTEXITCODE -ne 0) { throw "LaTeX compilation failed" }
    Copy-Item -LiteralPath (Join-Path $PSScriptRoot "build/main.pdf") -Destination (Join-Path $PSScriptRoot "manuscript.pdf") -Force
} finally { Pop-Location }
