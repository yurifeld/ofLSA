@echo off
rem Fig. 1 (straight cavity, mode amplitudes) on the 256^3 grid. Double-click, or: fig1.bat <case dir>
set "OFLSA_GRID=%~1"
if "%OFLSA_GRID%"=="" set "OFLSA_GRID=256"
call "%~dp0launch.bat" 1 "%OFLSA_GRID%"
