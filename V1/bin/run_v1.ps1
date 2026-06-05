<#
.SYNOPSIS
Runs the RiceG2G V1 candidate-gene workflow on Windows.

.DESCRIPTION
This launcher wraps V1\bin\riceg2g.py with PowerShell-friendly parameters.
It can extract FarmCPU significant sites and peaks, optionally run MLM and
WinQTLcart preprocessing, collect candidate genes using WinQTLcart intervals
or a 200kb window fallback, run Step2 scoring, annotate candidates, and write
Top10 candidate-gene recommendations.

Use -CandidateMode window for FarmCPU-only datasets such as data_TGW when no
matching WinQTLcart files are available. Use -CandidateMode both when both
WinQTLcart-backed and 200kb-window outputs are needed.

.PARAMETER Trait
Trait/file prefix, for example HNHZ or all_TGW_quality.

.PARAMETER Keyword
Step2 trait keyword, for example Heading_date, Awn_length, or TGW.

.PARAMETER Tissue
Main tissue name, for example young_panicle, leaves_flowering_stage, or developing_seeds.

.PARAMETER Step2Trait
Optional Step2 prefix when existing Step2 input/output files use a different name.

.PARAMETER WinQtlTrait
Optional WinQTLcart trait key, for example HZ_Heading_date when -Trait is HNHZ.

.PARAMETER CandidateMode
Candidate interval mode: auto, winqtl, window, or both.

.PARAMETER FarmCpuPThreshold
P-value threshold for full FarmCPU CSV filtering. The default is 1e-5.

.PARAMETER CandidateWindowBp
Window size around each FarmCPU peak for FarmCPU-only candidate genes. The default is 200000.

.PARAMETER FarmCpuWindowBp
Window size used to merge nearby significant FarmCPU sites into peaks. The default is 2000000.

.PARAMETER ForceFarmCpu
Regenerate <Trait>.FarmCPU.significant.tsv and <Trait>.FarmCPU.peak.

.PARAMETER ForceRap
Regenerate RAP_<Trait>.FarmCPUpeak_info candidate tables.

.PARAMETER SkipFarmCpu
Reuse an existing <Trait>.FarmCPU.peak.

.PARAMETER SkipMlm
Skip MLM peak generation.

.PARAMETER SkipWinQtl
Skip WinQTLcart preprocessing.

.PARAMETER SkipAnnotation
Skip annotation and Top10 recommendation outputs.

.PARAMETER AnnotationOutputPrefix
Prefix for <prefix>_top10_candidate_genes.tsv.

.EXAMPLE
powershell -ExecutionPolicy Bypass -File .\V1\bin\run_v1.ps1 -Trait HNHZ -Keyword Heading_date -Tissue leaves_flowering_stage -WinQtlTrait HZ_Heading_date -CandidateMode both -ForceFarmCpu -ForceRap

Runs HNHZ with both WinQTLcart-backed and 200kb-window candidate outputs.

.EXAMPLE
powershell -ExecutionPolicy Bypass -File .\V1\bin\run_v1.ps1 -Trait all_TGW_quality -Keyword TGW -Tissue developing_seeds -CandidateMode window -SkipMlm -SkipWinQtl -ForceFarmCpu -ForceRap -FarmCpuPThreshold 1e-5

Runs TGW from data_TGW/all_TGW_quality.FarmCPU.csv using the 200kb-window fallback.

.OUTPUTS
<Trait>.FarmCPU.significant.tsv
<Trait>.FarmCPU.peak
RAP_<Trait>.FarmCPUpeak_info
RAP_Step2_<Trait>.FarmCPUpeak_info
RAP_Step2_<Trait>.annotated.tsv
<Trait>_top10_candidate_genes.tsv
#>
param(
    [Parameter(Mandatory = $true)][string]$Trait,
    [Parameter(Mandatory = $true)][string]$Keyword,
    [Parameter(Mandatory = $true)][string]$Tissue,
    [string]$Step2Trait = "",
    [string]$WinQtlTrait = "",
    [ValidateSet("auto", "winqtl", "window", "both")][string]$CandidateMode = "auto",
    [switch]$SkipMlm,
    [switch]$SkipWinQtl,
    [switch]$SkipFarmCpu,
    [switch]$SkipAnnotation,
    [switch]$ForceFarmCpu,
    [switch]$ForceRap,
    [string]$AnnotationOutputPrefix = "",
    [int]$FarmCpuWindowBp = 2000000,
    [double]$FarmCpuPThreshold = 1e-5,
    [int]$CandidateWindowBp = 200000
)

$ErrorActionPreference = "Stop"

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
if ([string]::IsNullOrWhiteSpace($WinQtlTrait)) {
    $WinQtlTrait = $Trait
}

$ArgsList = @(
    (Join-Path $RepoRoot "V1\bin\riceg2g.py"),
    "--trait", $Trait,
    "--keyword", $Keyword,
    "--tissue", $Tissue,
    "--step2-trait", $Step2Trait,
    "--winqtl-trait", $WinQtlTrait,
    "--candidate-mode", $CandidateMode,
    "--farmcpu-window-bp", $FarmCpuWindowBp,
    "--farmcpu-p-threshold", $FarmCpuPThreshold,
    "--candidate-window-bp", $CandidateWindowBp
)

if ($SkipMlm) {
    $ArgsList += "--skip-mlm"
}
if ($SkipWinQtl) {
    $ArgsList += "--skip-winqtl"
}
if ($SkipFarmCpu) {
    $ArgsList += "--skip-farmcpu"
}
if ($ForceFarmCpu) {
    $ArgsList += "--force-farmcpu"
}
if ($ForceRap) {
    $ArgsList += "--force-rap"
}
if ($SkipAnnotation) {
    $ArgsList += "--skip-annotation"
}
if (-not [string]::IsNullOrWhiteSpace($AnnotationOutputPrefix)) {
    $ArgsList += @("--annotation-output-prefix", $AnnotationOutputPrefix)
}

& $PyExe @ArgsList
