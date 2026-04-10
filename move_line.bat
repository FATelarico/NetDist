@echo off
setlocal

if "%~3"=="" (
    echo Usage: %~nx0 "source.txt" "target string" "anchor string"
    exit /b 1
)

set "SRC=%~1"
set "TARGET=%~2"
set "ANCHOR=%~3"

powershell -NoProfile -ExecutionPolicy Bypass -Command ^
  "$ErrorActionPreference='Stop';" ^
  "$src=$env:SRC; $target=$env:TARGET; $anchor=$env:ANCHOR;" ^
  "if (-not (Test-Path -LiteralPath $src)) { throw 'Source file not found.' }" ^
  "$lines=[System.Collections.Generic.List[string]]::new();" ^
  "$lines.AddRange([string[]](Get-Content -LiteralPath $src -Encoding UTF8));" ^
  "$targetIndex=-1;" ^
  "for($i=0; $i -lt $lines.Count; $i++) { if($lines[$i].IndexOf($target,[System.StringComparison]::OrdinalIgnoreCase) -ge 0) { $targetIndex=$i; break } }" ^
  "if ($targetIndex -lt 0) { throw 'Target string not found in any line.' }" ^
  "$lines.RemoveAt($targetIndex);" ^
  "$anchorIndex=-1;" ^
  "for($i=0; $i -lt $lines.Count; $i++) { if($lines[$i].IndexOf($anchor,[System.StringComparison]::OrdinalIgnoreCase) -ge 0) { $anchorIndex=$i; break } }" ^
  "if ($anchorIndex -lt 0) { throw 'Anchor string not found in any line.' }" ^
  "$insertIndex=$anchorIndex-3; if ($insertIndex -lt 0) { $insertIndex=0 }" ^
  "$lines.Insert($insertIndex,$target);" ^
  "Set-Content -LiteralPath $src -Value $lines -Encoding UTF8;"

if errorlevel 1 (
    echo Failed.
    exit /b 1
)

echo Done.
exit /b 0