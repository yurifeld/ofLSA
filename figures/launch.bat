@echo off
rem Common launcher for fig*.bat:  launch.bat <figure 1-4> <256 | 50 | case dir>
setlocal
set "HERE=%~dp0"
set "FIG=%~1"
set "GRID=%~2"
for %%I in ("%HERE%..") do set "ROOT=%%~fI"

set "CASEDIR="
if "%GRID%"=="256" if %FIG% LEQ 2 set "CASEDIR=%ROOT%\tutorials\paper\straightCavity256"
if "%GRID%"=="50"  if %FIG% LEQ 2 set "CASEDIR=%ROOT%\tutorials\lidDrivenCavity3D"
if "%GRID%"=="256" if %FIG% GEQ 3 set "CASEDIR=%ROOT%\tutorials\paper\diagonalCavity256"
if "%GRID%"=="50"  if %FIG% GEQ 3 set "CASEDIR=%ROOT%\tutorials\diagonalCavity3D"
if not defined CASEDIR for %%I in ("%GRID%") do set "CASEDIR=%%~fI"

if exist "%CASEDIR%\1\leadingEigenModeRealU" goto haveData
if exist "%CASEDIR%\processor0\1\leadingEigenModeRealU" goto haveData
echo No results in "%CASEDIR%\1".
echo Run the case first, or copy the case folder with its 1\ directory from the cluster.
pause
exit /b 1
:haveData

rem ParaView: PATH first, then Program Files (the newest ParaView* folder wins)
set "PV="
for %%P in (paraview.exe) do if not "%%~$PATH:P"=="" set "PV=%%~$PATH:P"
if defined PV goto havePV
for /d %%D in ("%ProgramFiles%\ParaView*") do if exist "%%~fD\bin\paraview.exe" set "PV=%%~fD\bin\paraview.exe"
if defined PV goto havePV
echo ParaView not found. Install it from paraview.org or add its bin folder to PATH.
pause
exit /b 1
:havePV

set "OFLSA_CASE=%CASEDIR%"
set "OFLSA_FIGDIR=%HERE:~0,-1%"
echo Fig. %FIG% from %CASEDIR%
start "" "%PV%" --script="%HERE%fig%FIG%.py"
