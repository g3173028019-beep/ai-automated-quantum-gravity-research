param([string[]]$Paper)
$ErrorActionPreference = 'Stop'
$repoRoot = Split-Path -Parent $PSScriptRoot
$papers = if ($Paper) { $Paper } else {
    Get-ChildItem -LiteralPath (Join-Path $repoRoot 'papers') -Directory |
        Where-Object { Test-Path -LiteralPath (Join-Path $_.FullName 'main.tex') } |
        Sort-Object Name | Select-Object -ExpandProperty Name
}
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
