param(
    [string]$ComsolBin = 'D:\tools\comosol\COMSOL62\Multiphysics\bin\win64'
)
$ErrorActionPreference = 'Stop'
$taskSource = $PSScriptRoot
$taskOutput = Split-Path -Parent $taskSource
$taskWorkspace = Split-Path -Parent (Split-Path -Parent $taskOutput)
$taskWork = Join-Path $taskWorkspace 'work\q1_rebuild'
$taskPrefs = Join-Path $taskWork 'prefs'
New-Item -ItemType Directory -Force -Path $taskWork, $taskPrefs | Out-Null
@'
security.external.enable=on
security.external.filepermission=full
security.external.netpermission=off
security.external.runtimepermission=off
security.external.propertypermission=off
'@ | Set-Content -LiteralPath (Join-Path $taskPrefs 'comsol.prefs') -Encoding ascii
$taskDestination = ($taskOutput.Replace('\','/') + '/')
foreach ($taskName in @('Q1SolidWire','Q1Post')) {
    $taskJava = Get-Content -LiteralPath (Join-Path $taskSource ($taskName + '.java')) -Raw
    $taskOutLine = 'static final String OUT="' + $taskDestination + '";'
    $taskJava = [regex]::Replace($taskJava, 'static final String OUT="[^"]*";', [System.Text.RegularExpressions.MatchEvaluator]{param($match) $taskOutLine})
    $taskJavaPath = Join-Path $taskWork ($taskName + '.java')
    $taskClassPath = Join-Path $taskWork ($taskName + '.class')
    Set-Content -LiteralPath $taskJavaPath -Value $taskJava -Encoding ascii
    $taskStarted = Get-Date
    & (Join-Path $ComsolBin 'comsolcompile.exe') $taskJavaPath
    if (-not (Test-Path -LiteralPath $taskClassPath)) { throw "Compilation did not create $taskClassPath" }
    if ((Get-Item -LiteralPath $taskClassPath).LastWriteTime -lt $taskStarted.AddSeconds(-2)) { throw 'Compilation output is stale.' }
    $taskLog = Join-Path $taskWork ($taskName + '.log')
    & (Join-Path $ComsolBin 'comsolbatch.exe') -inputfile $taskClassPath -nosave -prefsdir $taskPrefs -batchlog $taskLog -np 4 -graphics -3drend sw
    $taskLogText = Get-Content -LiteralPath $taskLog -Raw
    if ($LASTEXITCODE -ne 0 -or $taskLogText.Contains('/*****')) { throw "COMSOL reported an error. See $taskLog" }
}
Write-Output "COMSOL models and raw data rebuilt in $taskOutput"
Write-Output 'To regenerate comparison figures, run analyze_q1.py with NumPy, SciPy, and Matplotlib installed.'
