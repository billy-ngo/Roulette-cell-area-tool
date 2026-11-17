@echo off
rem Double-click this file to launch the TEM Cell Area Tool on Windows.
rem Delayed expansion would strip "!" from folder and file names.
setlocal DisableDelayedExpansion
rem Remember this folder before changing directory: %~dp0 can resolve
rem differently afterwards.
set "HERE=%~dp0"
cd /d "%HERE%" 2>nul

if not exist "%HERE%cell_area_tool.py" (
    echo.
    echo cell_area_tool.py was not found next to this launcher.
    echo Extract the whole folder first, then run the launcher from there.
    echo.
    pause
    exit /b 1
)

rem Use the first Python 3 that has Tkinter, NumPy and Pillow.  "call" is
rem needed because "python" can be a batch file (pyenv-win).
set "PY="
for %%C in ("py -3" "python" "python3") do (
    if not defined PY (
        call %%~C -c "import tkinter, numpy, PIL" >nul 2>nul && set "PY=%%~C"
    )
)

if not defined PY (
    echo.
    echo Could not find a Python 3 with Tkinter, NumPy and Pillow.
    echo.
    echo Install Python 3 from https://www.python.org/downloads/ then run:
    echo     py -m pip install -r "%HERE%requirements.txt"
    echo.
    pause
    exit /b 1
)

rem Prefer pythonw.exe so that no console window stays open.  The path is
rem printed in the console's code page so that accented letters survive.
set "PYW="
for /f "usebackq delims=" %%P in (`%PY% -c "import os, sys, ctypes; cp = ctypes.windll.kernel32.GetConsoleOutputCP(); cp and sys.stdout.reconfigure(encoding='cp' + str(cp), errors='replace'); print(os.path.join(os.path.dirname(sys.executable), 'pythonw.exe'))"`) do set "PYW=%%P"

rem No parentheses from here on: a dropped file name may contain them.
if exist "%PYW%" start "" "%PYW%" "%HERE%cell_area_tool.py" %*
if exist "%PYW%" exit /b 0
%PY% "%HERE%cell_area_tool.py" %*
if errorlevel 1 pause
