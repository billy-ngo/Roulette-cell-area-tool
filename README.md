# Roulette cell area tool

A desktop tool for outlining cells in TEM (or any other) images and measuring
their area and perimeter in physical units. Several images can be open at once,
each with its own scale. Results are exported as one CSV table plus an annotated
JPEG per image.

It is a standalone Python/Tkinter program and does not require ImageJ.

## Install

Roulette needs Python 3.8 or newer with Tkinter, which the
[python.org](https://www.python.org/downloads/) installers for Windows and
macOS include. Installing it also installs NumPy, Pillow, SciPy and tifffile.

Windows:

```bash
py -m pip install roulette-cell-tool
```

macOS and Linux:

```bash
python3 -m pip install roulette-cell-tool
```

Platform notes:

- **macOS:** use Python from python.org. Apple's built-in `/usr/bin/python3`
  ships a Tk that crashes on recent macOS. Homebrew Python also works after
  `brew install python-tk`, but it refuses the `pip install` above; install
  into a virtual environment instead and start the tool with
  `~/roulette-env/bin/roulette`:

  ```bash
  python3 -m venv ~/roulette-env
  ~/roulette-env/bin/python3 -m pip install roulette-cell-tool
  ```

- **Linux:** install Tkinter from the distribution, e.g.
  `sudo apt install python3-tk python3-venv`. Recent Debian and Ubuntu also
  refuse the `pip install` above; use a virtual environment as for Homebrew.

## Running

```bash
roulette
```

`roulette image1.tif image2.tif` opens those images at start. If the command is
not found, the folder that pip installs commands into is not on `PATH`; start
the tool with `py -m cell_area_tool` on Windows or `python3 -m cell_area_tool`
on macOS and Linux instead.

### Running from a download

The tool is a single file, `cell_area_tool.py`, and also runs from a copy of
this repository without being installed. It then needs NumPy and Pillow (9.0 or
newer); SciPy (for *Detect scale bar*) and tifffile (for reading pixel size from
TIFF metadata, and for opening TIFF types that Pillow cannot read) are optional.

Install the packages from a terminal opened in the folder. On Windows, click
the address bar of the folder window, type `cmd` and press Enter, then run
`py -m pip install -r requirements.txt`. On macOS and Linux run
`python3 -m pip install -r requirements.txt`; where pip refuses (Homebrew,
Debian, Ubuntu), create the virtual environment in the folder
(`python3 -m venv .venv`, then
`.venv/bin/python3 -m pip install -r requirements.txt`), where the macOS
launcher looks for it.

Then start the tool:

- **Windows:** double-click `Launch Cell Area Tool.bat`.
- **macOS:** double-click `Launch Cell Area Tool.command`.
  - If macOS reports that the file cannot be executed, run
    `chmod +x "Launch Cell Area Tool.command"` once in Terminal from the folder.
  - If macOS blocks it as being from an unidentified developer, right-click it
    and choose *Open*. On macOS 15 and later, allow it under
    *System Settings → Privacy & Security → Open Anyway*.
- **From a terminal:** `py cell_area_tool.py [image ...]` on Windows,
  `python3 cell_area_tool.py [image ...]` on macOS and Linux (with a virtual
  environment, `.venv/bin/python3 cell_area_tool.py`).

## Workflow

1. **Open images** with `Open image(s)…`. TIFF, PNG, JPG and BMP are supported.
   16-bit and float TIFFs are contrast-stretched for display only. Only the
   first frame of a multi-page TIFF is used. Switch between open images with
   the `Image:` drop-down.
2. **Set the scale of each image**, using one of:
   - *File metadata.* Pixel size found in a TIFF (FEI/Thermo, Zeiss, ImageJ,
     OME-XML, description text, resolution tags) is applied when the image is
     opened. `Scale from file` re-reads it. Values that look like a print
     resolution (dpi), or that describe non-square pixels, are only applied
     after confirmation.
   - `Detect scale bar` finds the bar, highlights it in magenta and asks for
     its length. A bar thinner than 3 px is only picked when nothing thicker
     and bar-like is present; if something else is highlighted, cancel and use
     `Set scale (bar)`. The detected length can be a pixel or two short for
     bars with rounded or blurred ends; trace such bars yourself.
   - `Set scale (bar)`: click the two ends of the bar and enter its length.
   - `Pixel size…`: enter the pixel size directly.

   A warning appears if the new scale differs at least 100-fold from the other
   open images, which usually means the wrong unit was selected.
3. **Outline cells.** The `Draw:` button switches between two modes:
   - *Freehand:* drag around the cell and release to close the outline.
   - *Points:* click each vertex. Close the outline by clicking the first
     point, pressing `Enter`, right-clicking or double-clicking.
4. **Review.** Each outline is numbered and listed in the table with its area and
   perimeter. Click a row to highlight that cell; `Delete cell` removes it.
5. **Export** with `Export table + JPEGs…`.

Areas are reported in nm² and perimeters in nm. Images without a scale are
reported in pixels.

## Controls

| Action | Control |
| --- | --- |
| Zoom | Mouse wheel, `+` / `−` buttons, `+` / `-` keys |
| Pan | Middle-button drag, or `Pan tool` and left-drag |
| Fit image to window | `Fit` button or `f` |
| Undo last point (Points mode) | `Backspace` |
| Cancel current outline | `Esc` |
| Delete a cell | Click its row in the table, then `Delete` (the delete key on a Mac) or `Delete cell` |

## Output

- `cell_measurements.csv`: one row per cell with its number, image, pixel size,
  vertex count, area and perimeter (in pixels and in calibrated units), and a
  `self_intersecting` flag. Encoded as UTF-8 with a byte-order mark so that
  Excel displays the units correctly.
- `<image>_outlines.jpg`: the image with its numbered outlines.

The `image` column holds the file name only. If two open images share a file
name, the message shown after the export lists which cell numbers belong to
which file.

## Sessions

`Save session` writes the image paths, scales and outlines to a JSON file;
`Load session` restores them. A session file contains the full path of every
image. It still loads after the project folder is moved or copied to another
computer, as long as the images keep their position relative to the session
file or sit in the same folder as it. Images that cannot be found or opened are
listed when the session is loaded, and their outlines are not loaded.

## Measurement notes

- Area is computed from the outline vertices with the shoelace formula and
  converted with the square of that image's pixel size.
- Pixels are assumed to be square.
- Freehand outlines are simplified with the Ramer–Douglas–Peucker algorithm
  (tolerance 1.5 image pixels). This cuts the corners of a curved outline, so
  the measured area is slightly too small: by less than 1% for outlines more
  than about 300 px across, 1–2% at 150 px, 2–3% at 80 px and 4–7% at 40 px.
  Outlines drawn in Points mode are stored exactly as clicked.
- Outline points placed beyond the edge of the image are moved onto the edge,
  and the part of an outline that runs along the edge is stored as a straight
  line.
- An outline that crosses itself is marked `⚠` in the table and flagged in the
  CSV. The usual cause is a freehand trace that runs past its starting point
  before the button is released, which typically changes the area by less than
  1%. An outline that crosses itself anywhere else should be retraced.

## Tests

From a copy of this repository:

```bash
python3 test_geometry.py
```

(`py test_geometry.py` on Windows.) Tests that need SciPy or tifffile are
reported as `SKIP` when the package is not installed.

`make_test_image.py` generates `test_tem_image.png`, a synthetic image with a
500 nm scale bar and shapes of known area (listed in
`test_tem_image_truth.txt`), for trying the tool without real data.
