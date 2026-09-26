param(
    [string]$ComsolBin = 'D:\tools\comosol\COMSOL62\Multiphysics\bin\win64',
    [string]$OutputDirectory = '',
    [string]$PythonExe = 'D:\tools\anaconda\python.exe'
)
$ErrorActionPreference = 'Stop'
if (-not $OutputDirectory) {
    $OutputDirectory = Join-Path (Split-Path $PSScriptRoot -Parent) ('rebuild_' + (Get-Date -Format 'yyyyMMdd_HHmmss'))
}
$twistOut = [System.IO.Path]::GetFullPath($OutputDirectory)
$twistWork = Join-Path $env:TEMP ('Q2Twist_' + (Get-Date -Format 'yyyyMMdd_HHmmss'))
$twistLogs = Join-Path $twistOut 'verification\logs'
$twistPrefs = Join-Path $twistWork 'prefs'
New-Item -ItemType Directory -Force -Path $twistOut,$twistWork,$twistLogs,$twistPrefs | Out-Null
$utf8NoBom = New-Object System.Text.UTF8Encoding($false)
$prefsText = @'
security.comsol.allowbatch=on
security.comsol.allowmethods=on
security.external.enable=on
security.external.filepermission=full
security.external.netpermission=off
security.external.propertypermission=off
security.external.runtimepermission=off
'@
[System.IO.File]::WriteAllText((Join-Path $twistPrefs 'comsol.prefs'), $prefsText, $utf8NoBom)
$twistOutJava = $twistOut.Replace('\','/') + '/'
foreach ($twistName in @('Q2RegularTwist','Q2RegularTwistGeometry')) {
    $twistText = [System.IO.File]::ReadAllText((Join-Path $PSScriptRoot ($twistName+'.java')))
    $twistText = [regex]::Replace($twistText, 'String OUT="[^"]*"', ('String OUT="'+$twistOutJava+'"'))
    $twistJava = Join-Path $twistWork ($twistName+'.java')
    [System.IO.File]::WriteAllText($twistJava,$twistText,$utf8NoBom)
    & (Join-Path $ComsolBin 'comsolcompile.exe') $twistJava
    if ($LASTEXITCODE -ne 0) { throw "Compilation failed: $twistName" }
    $twistLog = Join-Path $twistLogs ($twistName+'.log')
    $twistConsole = Join-Path $twistLogs ($twistName+'_console.txt')
    & (Join-Path $ComsolBin 'comsolbatch.exe') -inputfile (Join-Path $twistWork ($twistName+'.class')) -nosave -prefsdir $twistPrefs -batchlog $twistLog -np 4 -graphics -3drend sw 2>&1 | Tee-Object -FilePath $twistConsole
    if ($LASTEXITCODE -ne 0) { throw "COMSOL failed: $twistName" }
    $expected = if ($twistName -eq 'Q2RegularTwist') { 'REGULAR_TWIST_COMPLETE' } else { 'REGULAR_TWIST_GEOMETRY_COMPLETE' }
    if (-not (Select-String -LiteralPath $twistConsole -SimpleMatch $expected -Quiet)) {
        throw "COMSOL did not report completion. See $twistLog"
    }
}
if (Test-Path -LiteralPath $PythonExe) {
    & $PythonExe (Join-Path $PSScriptRoot 'check_geometry.py') $twistOut
    if ($LASTEXITCODE -ne 0) { throw 'Geometry check failed.' }
    $twistReferenceSource = Join-Path (Split-Path $PSScriptRoot -Parent) 'references'
    $twistReferenceDestination = Join-Path $twistOut 'references'
    if ([System.IO.Path]::GetFullPath($twistReferenceSource) -ne [System.IO.Path]::GetFullPath($twistReferenceDestination)) {
        New-Item -ItemType Directory -Force -Path $twistReferenceDestination | Out-Null
        Get-ChildItem -LiteralPath $twistReferenceSource -File | Copy-Item -Destination $twistReferenceDestination
    }
    & $PythonExe (Join-Path $PSScriptRoot 'analyze_results.py') $twistOut
    if ($LASTEXITCODE -ne 0) { throw 'Result validation failed.' }
}
Write-Output "Rebuilt solved models and verification data: $twistOut"
