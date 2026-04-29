# scripts/replace_github_repo.ps1
<#
.SYNOPSIS
  Safely replace the contents of the benchtopc/Dual-Barcode-Demux GitHub repository
  with this clean workflow source tree.

.EXAMPLE
  powershell -ExecutionPolicy Bypass -File .\scripts\replace_github_repo.ps1

.EXAMPLE
  powershell -ExecutionPolicy Bypass -File .\scripts\replace_github_repo.ps1 `
    -RepoUrl "https://github.com/benchtopc/Dual-Barcode-Demux.git" `
    -CommitMessage "Replace corrupted files with clean v6 workflow"

.NOTES
  Requires Git and GitHub authentication.
  This script preserves the cloned repository's .git directory, deletes old tracked/source
  files, copies this clean workflow tree, commits, and pushes.
#>

param(
    [string]$RepoUrl = "https://github.com/benchtopc/Dual-Barcode-Demux.git",
    [string]$Branch = "main",
    [string]$CommitMessage = "Replace corrupted files with clean Dual-Barcode-Demux v6 workflow",
    [string]$WorkDir = "",
    [switch]$DryRun
)

$ErrorActionPreference = "Stop"

function Write-Step {
    param([string]$Message)
    Write-Host ""
    Write-Host "==> $Message" -ForegroundColor Cyan
}

function Resolve-WorkflowRoot {
    $scriptPath = $PSCommandPath
    if (-not $scriptPath) {
        $scriptPath = $MyInvocation.MyCommand.Path
    }

    if ($scriptPath) {
        $scriptDir = Split-Path -Parent $scriptPath
        $candidate = Resolve-Path (Join-Path $scriptDir "..")
        if (Test-Path (Join-Path $candidate "main.nf")) {
            return $candidate.Path
        }
    }

    $current = Get-Location
    if (Test-Path (Join-Path $current "main.nf")) {
        return $current.Path
    }

    throw "Cannot find workflow root. Run this from the unzipped workflow folder or from scripts/replace_github_repo.ps1."
}

function Copy-WorkflowContents {
    param(
        [string]$SourceRoot,
        [string]$DestinationRoot
    )

    $excludeNames = @(
        ".git",
        ".nextflow",
        ".nextflow.log",
        "work",
        "output",
        "local_test_output",
        "test_output"
    )

    Get-ChildItem -LiteralPath $SourceRoot -Force | ForEach-Object {
        if ($excludeNames -contains $_.Name) {
            return
        }

        $target = Join-Path $DestinationRoot $_.Name
        Copy-Item -LiteralPath $_.FullName -Destination $target -Recurse -Force
    }
}

function Clear-RepoContents {
    param([string]$RepoPath)

    Get-ChildItem -LiteralPath $RepoPath -Force | ForEach-Object {
        if ($_.Name -eq ".git") {
            return
        }
        Remove-Item -LiteralPath $_.FullName -Recurse -Force
    }
}

Write-Step "Checking Git"
git --version | Out-Host

$sourceRoot = Resolve-WorkflowRoot
Write-Step "Using workflow source root: $sourceRoot"

if (-not (Test-Path (Join-Path $sourceRoot "README.md"))) {
    throw "README.md not found in source root: $sourceRoot"
}
if (-not (Test-Path (Join-Path $sourceRoot "main.nf"))) {
    throw "main.nf not found in source root: $sourceRoot"
}
if (-not (Test-Path (Join-Path $sourceRoot "nextflow.config"))) {
    throw "nextflow.config not found in source root: $sourceRoot"
}
if (-not (Test-Path (Join-Path $sourceRoot "nextflow_schema.json"))) {
    throw "nextflow_schema.json not found in source root: $sourceRoot"
}

if ([string]::IsNullOrWhiteSpace($WorkDir)) {
    $timestamp = Get-Date -Format "yyyyMMdd_HHmmss"
    $WorkDir = Join-Path ([System.IO.Path]::GetTempPath()) "Dual-Barcode-Demux_replace_$timestamp"
}

if (Test-Path $WorkDir) {
    Remove-Item -LiteralPath $WorkDir -Recurse -Force
}
New-Item -ItemType Directory -Path $WorkDir | Out-Null

Write-Step "Cloning repository"
Push-Location $WorkDir
try {
    git clone $RepoUrl repo
    Set-Location repo

    git checkout $Branch

    Write-Step "Removing old repository contents except .git"
    if (-not $DryRun) {
        Clear-RepoContents -RepoPath (Get-Location).Path
    }

    Write-Step "Copying clean workflow contents"
    if (-not $DryRun) {
        Copy-WorkflowContents -SourceRoot $sourceRoot -DestinationRoot (Get-Location).Path
    }

    Write-Step "Validating required root files"
    foreach ($required in @("README.md", "main.nf", "nextflow.config", "nextflow_schema.json")) {
        if (-not (Test-Path $required)) {
            throw "Required root file missing after copy: $required"
        }
    }

    Write-Step "Git status"
    git status --short | Out-Host

    if ($DryRun) {
        Write-Host "Dry run complete. No commit or push was performed." -ForegroundColor Yellow
        exit 0
    }

    $changes = git status --porcelain
    if ([string]::IsNullOrWhiteSpace($changes)) {
        Write-Host "No changes detected. Repository already matches this workflow." -ForegroundColor Green
        exit 0
    }

    Write-Step "Committing changes"
    git add -A
    git commit -m $CommitMessage

    Write-Step "Pushing to GitHub"
    git push origin $Branch

    Write-Host ""
    Write-Host "Done. Now import this URL in EPI2ME:" -ForegroundColor Green
    Write-Host "https://github.com/benchtopc/Dual-Barcode-Demux"
}
finally {
    Pop-Location
}
