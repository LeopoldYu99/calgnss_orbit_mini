$ErrorActionPreference = "Stop"
$projectDir = Split-Path -Parent $PSScriptRoot
$requirements = Join-Path $PSScriptRoot "requirements-windows.txt"
$tester = Join-Path $PSScriptRoot "windows_mq_tester.py"

py -3 -c "import paramiko" 2>$null
if ($LASTEXITCODE -ne 0) {
    Write-Host "Installing the Windows tester dependency (Paramiko)..."
    py -3 -m pip install --user -r $requirements
}

Set-Location $projectDir
py -3 $tester
