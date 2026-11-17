#!/bin/bash
# ============================================================
#  Double-click this file (in Finder) to launch the TEM Cell
#  Area Tool on macOS.
#
#  If double-clicking is blocked ("cannot be opened because it
#  is from an unidentified developer" or "permission denied"),
#  run this once in Terminal from this folder:
#      chmod +x "Launch Cell Area Tool.command"
#  and, if Gatekeeper complains, right-click the file -> Open
#  (macOS 15 and later: System Settings -> Privacy & Security
#  -> Open Anyway).
# ============================================================

# Work from the folder this script lives in, so relative paths resolve.
cd "$(dirname "$0")" || exit 1

SCRIPT="cell_area_tool.py"

if [ ! -f "$SCRIPT" ]; then
    echo
    echo "$SCRIPT was not found next to this launcher."
    echo "Keep the launcher in the same folder as the rest of the tool."
    echo
    read -n 1 -s -r -p "Press any key to close..."
    exit 1
fi

# A usable interpreter must (a) have Tk 8.6 or newer, (b) create a Tk window
# without aborting and (c) import numpy + Pillow. Apple's system python3 ships
# Tcl/Tk 8.5, which is buggy and aborts on recent macOS ("macOS ... or later
# required"), so we probe each candidate for real rather than trusting the
# first python3 on PATH. The subshell + redirect hides the "Abort trap: 6"
# notice for candidates that crash, so a failed probe is silent.
probe() {
    ( "$1" - <<'PY'
import sys, tkinter, numpy, PIL
if tkinter.TkVersion < 8.6:
    sys.exit(1)
tkinter.Tk().destroy()
PY
    ) >/dev/null 2>&1
}

# Candidate interpreters, best first: a virtual environment in this folder,
# python.org frameworks (bundle Tk 8.6+), then versioned Homebrew / local
# builds, then whatever is on PATH. Homebrew often installs only a versioned
# name (python3.11), with no bare "python3".
CANDIDATES=("./.venv/bin/python3" "./venv/bin/python3")
for v in 3.14 3.13 3.12 3.11 3.10; do
    CANDIDATES+=("/Library/Frameworks/Python.framework/Versions/$v/bin/python3")
    CANDIDATES+=("/opt/homebrew/bin/python$v")
    CANDIDATES+=("/usr/local/bin/python$v")
done
CANDIDATES+=("/opt/homebrew/bin/python3" "/usr/local/bin/python3")
CANDIDATES+=("$(command -v python3 2>/dev/null)" "$(command -v python 2>/dev/null)")

PY=""
for c in "${CANDIDATES[@]}"; do
    [ -n "$c" ] && [ -x "$c" ] || continue
    if probe "$c"; then PY="$c"; break; fi
done

if [ -z "$PY" ]; then
    echo
    echo "No Python 3 with a working Tk GUI and NumPy + Pillow was found."
    echo
    echo "Fix (recommended): install the official build from"
    echo "    https://www.python.org/downloads/macos/"
    echo "then, in Terminal, from this folder:"
    echo "    python3 -m pip install -r requirements.txt"
    echo
    echo "Or with Homebrew, which needs a virtual environment in this folder:"
    echo "    brew install python-tk"
    echo "    python3 -m venv .venv"
    echo "    .venv/bin/python3 -m pip install -r requirements.txt"
    echo
    echo "(Apple's built-in python3 uses an old Tk that crashes on this macOS,"
    echo " so it cannot run this tool.)"
    echo
    read -n 1 -s -r -p "Press any key to close..."
    exit 1
fi

# Launch, passing through any files given as arguments.
"$PY" "$SCRIPT" "$@"
STATUS=$?

# If the app crashed, keep the window open so the error stays visible.
if [ $STATUS -ne 0 ]; then
    echo
    echo "The tool exited with an error (status $STATUS)."
    echo "See cell_area_tool_error.log in this folder for details."
    echo
    read -n 1 -s -r -p "Press any key to close..."
fi
