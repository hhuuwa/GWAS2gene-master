$ErrorActionPreference = "Stop"

$RepoRoot = Resolve-Path (Join-Path $PSScriptRoot "..\..")
$PortablePy = Join-Path $RepoRoot ".tools\python\python\python.exe"

if (Test-Path $PortablePy) {
    $PyExe = $PortablePy
} else {
    $PyExe = "python"
}

Write-Host "[INFO] Using Python: $PyExe"
Write-Host "[INFO] Running RiceG2G V1 demo (skip SelectFarmCPUPeakByCor.pl)"

& $PyExe (Join-Path $RepoRoot "V1\bin\riceg2g.py") `
    --trait HZ_Awn_length `
    --keyword Awn_length `
    --tissue young_panicle `
    --step2-trait HZ_Awn_length

Write-Host "[DONE] Demo output:"
Write-Host (Join-Path $RepoRoot "RAP_Step2_HZ_Awn_length.FarmCPUpeak_info")
