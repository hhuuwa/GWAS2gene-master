$ErrorActionPreference = "Stop"

param(
    [Parameter(Mandatory = $true)][string]$Trait,
    [Parameter(Mandatory = $true)][string]$Keyword,
    [Parameter(Mandatory = $true)][string]$Tissue,
    [string]$Step2Trait = ""
)

$RepoRoot = Resolve-Path (Join-Path $PSScriptRoot "..\..")
$PortablePy = Join-Path $RepoRoot ".tools\python\python\python.exe"

if (Test-Path $PortablePy) {
    $PyExe = $PortablePy
} else {
    $PyExe = "python"
}

if ([string]::IsNullOrWhiteSpace($Step2Trait)) {
    $Step2Trait = $Trait
}

& $PyExe (Join-Path $RepoRoot "V1\bin\riceg2g.py") `
    --trait $Trait `
    --keyword $Keyword `
    --tissue $Tissue `
    --step2-trait $Step2Trait
