$ErrorActionPreference = "Stop"

$RepoRoot = Resolve-Path (Join-Path $PSScriptRoot "..\..")
$PortablePy = Join-Path $RepoRoot ".tools\python\python\python.exe"

if (Test-Path $PortablePy) {
    $PyExe = $PortablePy
} else {
    $PyExe = "python"
}

Write-Host "[INFO] Running V1 demo pipeline + visualization"

& $PyExe (Join-Path $RepoRoot "V1\bin\riceg2g.py") `
    --trait HZ_Awn_length `
    --keyword Awn_length `
    --tissue young_panicle `
    --step2-trait HZ_Awn_length `
    --demo-viz

Write-Host "[DONE] Open report:"
Write-Host (Join-Path $RepoRoot "V1\demo\results\report.html")
