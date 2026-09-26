param(
    [string]$ComsolBin = 'D:\tools\comosol\COMSOL62\Multiphysics\bin\win64',
    [string]$OutputDirectory = '',
    [string]$PythonExe = 'D:\tools\anaconda\python.exe',
    [string]$WorkDirectory = ''
)
$ErrorActionPreference = 'Stop'
if (-not $OutputDirectory) {
    $OutputDirectory = Join-Path (Split-Path $PSScriptRoot -Parent) ('rebuild_' + (Get-Date -Format 'yyyyMMdd_HHmmss'))
}
$q2Out = [System.IO.Path]::GetFullPath($OutputDirectory)
if (-not $WorkDirectory) { $WorkDirectory = Join-Path $env:TEMP ('Q2Circular_' + (Get-Date -Format 'yyyyMMdd_HHmmss')) }
$q2Work = [System.IO.Path]::GetFullPath($WorkDirectory)
$q2Logs = Join-Path $q2Out 'verification\logs'
$q2Prefs = Join-Path $q2Work 'prefs'
New-Item -ItemType Directory -Force -Path $q2Out,$q2Work,$q2Logs,$q2Prefs | Out-Null
$utf8NoBom = New-Object System.Text.UTF8Encoding($false)
# File access is needed by the CSV and MPH exports. Only this run's preference
# directory is changed; installed COMSOL and personal preferences are unchanged.
$prefsText = @'
security.comsol.allowbatch=on
security.comsol.allowmethods=on
security.external.enable=on
security.external.filepermission=full
security.external.netpermission=off
security.external.propertypermission=off
security.external.runtimepermission=off
'@
[System.IO.File]::WriteAllText((Join-Path $q2Prefs 'comsol.prefs'), $prefsText, $utf8NoBom)
$q2OutJava = $q2Out.Replace('\','/') + '/'
$q2Seed = [System.IO.Path]::GetFullPath((Split-Path $PSScriptRoot -Parent))
if ($q2Out -ne $q2Seed) {
    # These references have identical copper area and are reused, not rerun.
    foreach ($q2Ref in @('Q2_EqualAreaSolid_200kHz.mph','equal_area_solid_results.csv','layout_parameters.json','circular_layout.png')) {
        Copy-Item -LiteralPath (Join-Path $q2Seed $q2Ref) -Destination (Join-Path $q2Out $q2Ref)
    }
    Copy-Item -LiteralPath (Join-Path $q2Seed 'reference_hexagonal') -Destination (Join-Path $q2Out 'reference_hexagonal') -Recurse
}
foreach ($q2Name in @('Q2CircularBundle','Q2CircularGeometry','Q2CircularBoundaryCheck','Q2CircularGeometryImage')) {
    Write-Output "Building $q2Name"
    $q2Text = [System.IO.File]::ReadAllText((Join-Path $PSScriptRoot ($q2Name+'.java')))
    $q2Text = [regex]::Replace($q2Text, 'String OUT="[^"]*"', ('String OUT="'+$q2OutJava+'"'))
    $q2Java = Join-Path $q2Work ($q2Name+'.java')
    [System.IO.File]::WriteAllText($q2Java,$q2Text,$utf8NoBom)
    & (Join-Path $ComsolBin 'comsolcompile.exe') $q2Java
    if ($LASTEXITCODE -ne 0 -or -not (Test-Path -LiteralPath (Join-Path $q2Work ($q2Name+'.class')))) {
        throw "Compilation failed: $q2Name"
    }
    $q2Log = Join-Path $q2Logs ($q2Name+'.log')
    & (Join-Path $ComsolBin 'comsolbatch.exe') -inputfile (Join-Path $q2Work ($q2Name+'.class')) -nosave -prefsdir $q2Prefs -batchlog $q2Log -np 4 -graphics -3drend sw
    if ($LASTEXITCODE -ne 0) { throw "COMSOL process failed: $q2Name" }
    $q2LogText = [System.IO.File]::ReadAllText($q2Log)
    if ($q2LogText -match 'IMAGE_ERROR|Exception:|java\.lang\.[A-Za-z]*Exception|com\.comsol\.util\.exceptions') {
        throw "COMSOL reported an exception. See $q2Log"
    }
}
if ((Test-Path -LiteralPath $PythonExe) -and (Test-Path -LiteralPath (Join-Path $PSScriptRoot 'analyze_results.py'))) {
    & $PythonExe (Join-Path $PSScriptRoot 'analyze_results.py') $q2Out
    if ($LASTEXITCODE -ne 0) { throw 'Numerical validation or figure generation failed.' }
}
Write-Output "Q2 rebuild completed: $q2Out"
