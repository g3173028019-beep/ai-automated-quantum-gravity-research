$ErrorActionPreference = 'Stop'
$repoRoot = Split-Path -Parent $PSScriptRoot
$papers = @('01-projection-sign-rigidity', '02-half-wormholes-wiener-chaos')
foreach ($paper in $papers) {
    $paperDir = Join-Path $repoRoot ('papers/' + $paper)
    Push-Location $paperDir
    try {
        & latexmk -pdf -interaction=nonstopmode -halt-on-error -outdir=build main.tex
        if ($LASTEXITCODE -ne 0) { throw "LaTeX build failed: $paper" }
        Copy-Item -LiteralPath 'build/main.pdf' -Destination 'manuscript.pdf' -Force
    } finally {
        Pop-Location
    }
}
