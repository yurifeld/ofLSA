@echo off
rem Fig. 2 (straight cavity, real/imaginary parts) on the 50^3 grid. Double-click, or: fig2.bat <case dir>
set "OFLSA_GRID=%~1"
if "%OFLSA_GRID%"=="" set "OFLSA_GRID=50"
call "%~dp0launch.bat" 2 "%OFLSA_GRID%"
