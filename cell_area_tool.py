#!/usr/bin/env python
"""
TEM Cell Area Tool
==================

An interactive program for outlining cells in electron-microscopy (or any)
images and measuring their area and perimeter in real physical units.

Workflow
--------
1. Open one or more images (TIFF / PNG / JPG / BMP ...).
2. Calibrate each image: from file metadata, from its scale bar, or by entering
   the pixel size directly.
3. Trace cells freehand or point by point.
4. Export a CSV of measurements and an annotated JPEG per image; save the
   session as JSON to resume later.

The geometry and calibration math live in pure functions near the top of this
file so they can be unit-tested without a display (see test_geometry.py).

Run
---
    python3 cell_area_tool.py               # start empty
    python3 cell_area_tool.py image.tif     # open image(s) on launch

(``py`` instead of ``python3`` on Windows.)

Requires Python 3.8+ with Tkinter, NumPy and Pillow 9.0+.  SciPy (scale-bar
detection) and tifffile (TIFF pixel-size metadata, and TIFF types that Pillow
cannot read) are optional.
"""

from __future__ import annotations

import csv
import json
import math
import os
import re
import sys
from dataclasses import dataclass, field
from typing import List, Optional, Sequence, Tuple


def _pip_command(package: str) -> str:
    """The command that installs ``package`` for the Python running this tool.

    A bare ``pip`` may belong to another Python installation, or be missing.
    """
    exe = sys.executable or "python3"
    folder, name = os.path.split(exe)
    if name.lower() == "pythonw.exe":       # no console: pip's output would be lost
        exe = os.path.join(folder, "python.exe")
    if re.search(r"[^\w.:\\/+~-]", exe):        # a space, &, ( ...: must be quoted
        only_spaces = not re.search(r"[^\w.:\\/+~ -]", exe)
        exe = f'"{exe}"'
        # PowerShell does not run a quoted path that is followed by arguments.
        # Through "cmd /c" the line works in PowerShell and in cmd alike, as
        # long as the path holds none of cmd's own special characters.
        if os.name == "nt" and only_spaces:
            exe = "cmd /c " + exe
    return f"{exe} -m pip install {package}"


# ---------------------------------------------------------------------------
# Dependency checks with friendly messages
# ---------------------------------------------------------------------------
try:
    import numpy as np
except ImportError:  # pragma: no cover
    sys.stderr.write(f"This tool needs NumPy.  Install it with:  {_pip_command('numpy')}\n")
    raise

try:
    from PIL import Image, ImageDraw, ImageFont
except ImportError:  # pragma: no cover
    sys.stderr.write(f"This tool needs Pillow.  Install it with:  {_pip_command('Pillow')}\n")
    raise


# Candidate bold fonts for the numbered labels burned into exported JPEGs,
# tried in order.  The bare names work on Windows (Pillow searches the Fonts
# folder); the absolute paths cover macOS and common Linux locations, where the
# Windows-style names are not found.  Falls back to Pillow's bitmap default.
_LABEL_FONT_CANDIDATES = (
    "arialbd.ttf",                                            # Windows
    "arial.ttf",
    "/System/Library/Fonts/Supplemental/Arial Bold.ttf",     # macOS
    "/System/Library/Fonts/Supplemental/Arial.ttf",
    "/System/Library/Fonts/Helvetica.ttc",
    "/Library/Fonts/Arial.ttf",
    "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf",  # Linux
    "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
    "DejaVuSans-Bold.ttf",
    "DejaVuSans.ttf",
)


def load_label_font(size: int):
    """Return a TrueType label font at ``size`` px, trying platform candidates.

    Falls back to Pillow's small built-in font if no candidate is found.
    """
    for name in _LABEL_FONT_CANDIDATES:
        try:
            return ImageFont.truetype(name, size)
        except Exception:
            continue
    return ImageFont.load_default()


Point = Tuple[float, float]

# Areas are always reported in nm^2, so every calibration is normalised to
# nanometres per pixel.  This maps the unit a scale bar may be labelled in to
# its length in nanometres.
LENGTH_UNITS_NM = {"nm": 1.0, "µm": 1000.0, "um": 1000.0, "mm": 1_000_000.0}


def to_nm(value: float, unit: str) -> float:
    """Convert a length in nm/µm/mm to nanometres."""
    return value * LENGTH_UNITS_NM.get(unit, 1.0)


# ===========================================================================
# Pure geometry / calibration helpers (unit-tested in test_geometry.py)
# ===========================================================================
def polygon_area(points: Sequence[Point]) -> float:
    """Absolute area of a simple polygon using the shoelace formula.

    ``points`` are (x, y) vertices in any consistent unit; the returned area is
    in that unit squared.  The polygon is treated as closed (the last vertex is
    joined back to the first).  Fewer than three vertices give an area of 0.

    Note: the shoelace formula is only correct for *simple* (non
    self-intersecting) polygons; see :func:`polygon_is_simple`.
    """
    n = len(points)
    if n < 3:
        return 0.0
    s = 0.0
    for i in range(n):
        x1, y1 = points[i]
        x2, y2 = points[(i + 1) % n]
        s += x1 * y2 - x2 * y1
    return abs(s) / 2.0


def polygon_perimeter(points: Sequence[Point], closed: bool = True) -> float:
    """Total edge length of a polyline/polygon.

    With ``closed=True`` the outline is treated as a closed polygon (the last
    vertex connects back to the first).
    """
    n = len(points)
    if n < 2:
        return 0.0
    total = 0.0
    last = n if closed else n - 1
    for i in range(last):
        x1, y1 = points[i]
        x2, y2 = points[(i + 1) % n]
        total += math.hypot(x2 - x1, y2 - y1)
    return total


def _seg_orientation(a: Point, b: Point, c: Point) -> float:
    """Cross product (b-a) x (c-a); sign gives turn direction."""
    return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])


def _on_segment(a: Point, b: Point, c: Point) -> bool:
    """True if collinear point c lies on segment ab."""
    return (
        min(a[0], b[0]) <= c[0] <= max(a[0], b[0])
        and min(a[1], b[1]) <= c[1] <= max(a[1], b[1])
    )


def segments_intersect(p1: Point, p2: Point, p3: Point, p4: Point) -> bool:
    """Return True if segment p1p2 crosses segment p3p4 (including touching)."""
    d1 = _seg_orientation(p3, p4, p1)
    d2 = _seg_orientation(p3, p4, p2)
    d3 = _seg_orientation(p1, p2, p3)
    d4 = _seg_orientation(p1, p2, p4)
    if ((d1 > 0) != (d2 > 0)) and ((d3 > 0) != (d4 > 0)):
        return True
    if d1 == 0 and _on_segment(p3, p4, p1):
        return True
    if d2 == 0 and _on_segment(p3, p4, p2):
        return True
    if d3 == 0 and _on_segment(p1, p2, p3):
        return True
    if d4 == 0 and _on_segment(p1, p2, p4):
        return True
    return False


def polygon_is_simple(points: Sequence[Point]) -> bool:
    """True if no two non-adjacent edges of the closed polygon intersect.

    Self-intersecting ("bow-tie") polygons make the shoelace area meaningless,
    so the GUI warns when this returns False.  O(n^2), fine for hand-traced
    outlines.
    """
    n = len(points)
    if n < 4:
        return True
    for i in range(n):
        a1 = points[i]
        a2 = points[(i + 1) % n]
        for j in range(i + 1, n):
            # Skip adjacent edges (they legitimately share a vertex).
            if (j + 1) % n == i or (i + 1) % n == j:
                continue
            b1 = points[j]
            b2 = points[(j + 1) % n]
            if segments_intersect(a1, a2, b1, b2):
                return False
    return True


def point_line_distance(p: Point, a: Point, b: Point) -> float:
    """Perpendicular distance from point ``p`` to the infinite line through a, b."""
    ax, ay = a
    bx, by = b
    px, py = p
    dx, dy = bx - ax, by - ay
    if dx == 0 and dy == 0:
        return math.hypot(px - ax, py - ay)
    return abs(dy * px - dx * py + bx * ay - by * ax) / math.hypot(dx, dy)


def rdp_simplify(points: Sequence[Point], epsilon: float) -> List[Point]:
    """Ramer–Douglas–Peucker polyline simplification (iterative, stack-based).

    Removes vertices that lie within ``epsilon`` of the line between kept
    neighbours.  Used to tidy dense freehand traces without changing their shape
    (and therefore barely changing the measured area).  Fewer than 3 points are
    returned unchanged.
    """
    pts = list(points)
    n = len(pts)
    if n < 3 or epsilon <= 0:
        return pts
    keep = [False] * n
    keep[0] = keep[n - 1] = True
    stack = [(0, n - 1)]
    while stack:
        s, e = stack.pop()
        dmax, idx = 0.0, -1
        for i in range(s + 1, e):
            d = point_line_distance(pts[i], pts[s], pts[e])
            if d > dmax:
                dmax, idx = d, i
        if idx != -1 and dmax > epsilon:
            keep[idx] = True
            stack.append((s, idx))
            stack.append((idx, e))
    return [pts[i] for i in range(n) if keep[i]]


def dedupe_polygon(points: Sequence[Point], tol: float = 0.5) -> List[Point]:
    """Drop consecutive (and wrap-around) vertices closer than ``tol``.

    Removes the coincident vertices that a click-then-double-click can produce,
    so a degenerate zero-area outline is never created.
    """
    out: List[Point] = []
    for p in points:
        p = (float(p[0]), float(p[1]))
        if not out or math.hypot(p[0] - out[-1][0], p[1] - out[-1][1]) > tol:
            out.append(p)
    # The polygon is closed, so a final vertex on top of the first is redundant.
    if len(out) >= 2 and math.hypot(out[-1][0] - out[0][0], out[-1][1] - out[0][1]) <= tol:
        out.pop()
    return out


def collapse_edge_runs(points: Sequence[Point], width: float, height: float) -> List[Point]:
    """Drop every vertex whose two neighbours lie with it on one edge of the image.

    Outline points beyond the image are moved onto its edge.  Where the pointer
    went further along the edge outside the image than the places it left and
    came back, that leaves a run of vertices on the edge that doubles back on
    itself: no area, but extra perimeter and an outline that touches itself.
    Only the two ends of each run are kept.  The vertices dropped lie on a
    straight line with their neighbours, so the area is unchanged.
    """
    def edges(p: Point) -> Tuple[bool, bool, bool, bool]:
        return (p[0] == 0.0, p[0] == width, p[1] == 0.0, p[1] == height)

    pts = list(points)
    changed = True
    while changed:      # repeat: removing a run can bring two other runs together
        changed = False
        i = 0
        while i < len(pts) and len(pts) > 2:
            before, here, after = edges(pts[i - 1]), edges(pts[i]), edges(pts[(i + 1) % len(pts)])
            if any(a and b and c for a, b, c in zip(before, here, after)):
                del pts[i]
                changed = True
            else:
                i += 1
    return pts


def area_in_units(area_px: float, unit_per_px: float) -> float:
    """Convert an area in px^2 to physical units^2.

    Area scales with the *square* of the linear pixel size.
    """
    return area_px * (unit_per_px ** 2)


def length_in_units(length_px: float, unit_per_px: float) -> float:
    """Convert a length in px to physical units."""
    return length_px * unit_per_px


def fit_transform(
    img_w: int, img_h: int, canvas_w: int, canvas_h: int, margin: float = 0.98
) -> Tuple[float, Tuple[float, float]]:
    """Compute a (scale, offset) that fits the image centred in the canvas.

    ``offset`` is the canvas coordinate of image pixel (0, 0).  Returns a scale
    and an (ox, oy) such that the whole image is visible and centred.
    """
    if img_w <= 0 or img_h <= 0 or canvas_w <= 0 or canvas_h <= 0:
        return 1.0, (0.0, 0.0)
    scale = min(canvas_w / img_w, canvas_h / img_h) * margin
    ox = (canvas_w - img_w * scale) / 2.0
    oy = (canvas_h - img_h * scale) / 2.0
    return scale, (ox, oy)


def canvas_to_image(cx: float, cy: float, scale: float, offset: Tuple[float, float]) -> Point:
    """Map a canvas coordinate to an image coordinate."""
    return (cx - offset[0]) / scale, (cy - offset[1]) / scale


def image_to_canvas(ix: float, iy: float, scale: float, offset: Tuple[float, float]) -> Point:
    """Map an image coordinate to a canvas coordinate."""
    return ix * scale + offset[0], iy * scale + offset[1]


def zoom_about(
    cx: float,
    cy: float,
    old_scale: float,
    offset: Tuple[float, float],
    factor: float,
    min_scale: float = 0.02,
    max_scale: float = 64.0,
) -> Tuple[float, Tuple[float, float]]:
    """Zoom by ``factor`` keeping the image point under (cx, cy) fixed.

    Returns the new (scale, offset).
    """
    new_scale = max(min_scale, min(max_scale, old_scale * factor))
    # Image point currently under the cursor.
    ix, iy = canvas_to_image(cx, cy, old_scale, offset)
    # Choose offset so that image point maps back to the same canvas point.
    ox = cx - ix * new_scale
    oy = cy - iy * new_scale
    return new_scale, (ox, oy)


def normalize_for_display(arr: "np.ndarray", low_pct: float = 0.5, high_pct: float = 99.5) -> "np.ndarray":
    """Convert an image array of any dtype to a display-ready uint8 array.

    16-bit and float TEM images are contrast-stretched between the given
    percentiles so faint structure is visible.  8-bit images pass through
    unchanged.  Colour images keep their channels; single-channel images stay
    2-D.  Pixel *values* never affect measurements, so this is purely cosmetic.
    """
    a = np.asarray(arr)
    if a.dtype == np.uint8:
        return a
    a = a.astype(np.float64)
    finite = a[np.isfinite(a)]
    if finite.size == 0:
        return np.zeros(a.shape, dtype=np.uint8)
    lo = np.percentile(finite, low_pct)
    hi = np.percentile(finite, high_pct)
    if hi <= lo:
        lo, hi = float(finite.min()), float(finite.max())
    if hi <= lo:
        return np.zeros(a.shape, dtype=np.uint8)
    a = np.clip((a - lo) / (hi - lo), 0.0, 1.0)
    return (a * 255.0 + 0.5).astype(np.uint8)


def format_measurement(value: float, unit: str) -> str:
    """Format a measurement with ~4 significant figures and thousands grouping.

    Very large (>= 1e6) or very small (< 1e-3) magnitudes fall back to
    scientific notation so the string stays short.
    """
    if value == 0:
        return "0 " + unit
    mag = abs(value)
    # nm^2 areas are naturally large; keep them as grouped integers up to 1e7
    # and only fall back to scientific notation beyond that (or for tiny values).
    if mag >= 1e7 or mag < 1e-3:
        return f"{value:.3e} {unit}"
    # Decimal places needed for ~4 significant figures at this magnitude.
    decimals = max(0, 3 - int(math.floor(math.log10(mag))))
    return f"{value:,.{decimals}f} {unit}"


# ===========================================================================
# Data model
# ===========================================================================
@dataclass
class Calibration:
    """Maps pixels to physical units.  ``unit_per_px`` == None means uncalibrated."""

    unit_per_px: Optional[float] = None
    unit: str = "px"

    @property
    def calibrated(self) -> bool:
        return self.unit_per_px is not None and self.unit_per_px > 0

    def area_unit(self) -> str:
        return (self.unit + "²") if self.calibrated else "px²"

    def length_unit(self) -> str:
        return self.unit if self.calibrated else "px"


@dataclass
class Cell:
    """One traced outline, stored in *image* pixel coordinates.

    Cells are numbered by their position in the session, not by a stored label,
    so numbers stay contiguous as cells are added and removed.
    """

    points: List[Point] = field(default_factory=list)
    # Cached result of the O(n^2) self-intersection test; a cell's points are
    # never edited after it is created.
    _simple: Optional[bool] = field(default=None, init=False, repr=False, compare=False)

    def area_px(self) -> float:
        return polygon_area(self.points)

    def perimeter_px(self) -> float:
        return polygon_perimeter(self.points, closed=True)

    def is_simple(self) -> bool:
        if self._simple is None:
            self._simple = polygon_is_simple(self.points)
        return self._simple


@dataclass
class ImageDoc:
    """One open image plus its own calibration, cells and view state.

    Each image carries its *own* calibration so areas are computed with the
    correct pixel size for that image, even when several images are open at
    once.
    """

    path: str
    pil_image: "Image.Image"
    width: int
    height: int
    calibration: Calibration = field(default_factory=Calibration)
    cells: List[Cell] = field(default_factory=list)
    scale_source: str = ""   # where the calibration came from (metadata / manual)
    # Per-image view so switching images returns you to where you were.
    scale: float = 1.0
    offset: Tuple[float, float] = (0.0, 0.0)
    fitted: bool = False

    @property
    def name(self) -> str:
        return os.path.basename(self.path)


def _is_tiff_file(path: str) -> bool:
    """True if ``path`` starts with a TIFF / BigTIFF byte-order mark."""
    try:
        with open(path, "rb") as f:
            return f.read(4) in (b"II*\x00", b"MM\x00*", b"II+\x00", b"MM\x00+")
    except OSError:
        return False


def _display_image(img: "Image.Image") -> "Image.Image":
    """Return a detached copy of ``img`` for display: 8-bit grayscale or RGB.

    Grayscale stays single-channel, which takes a quarter of the memory of RGB
    (Pillow stores RGB at 4 bytes per pixel) and matters when a session holds
    a hundred micrographs.  Pixel values are identical either way.
    """
    return img.copy() if img.mode == "L" else img.convert("RGB")


def load_image(path: str) -> "Image.Image":
    """Load ``path`` as an 8-bit grayscale or RGB PIL image ready for display.

    Handles 8/16-bit and float TIFFs (contrast-stretched), multi-frame TIFFs
    (first frame), palette and grayscale images.  Falls back to ``tifffile`` for
    exotic TIFFs that Pillow cannot open.
    """
    try:
        img = Image.open(path)
        img.load()
    except Exception as exc:
        # Fall back to tifffile for TEM-specific TIFF flavours.  For anything
        # that is not a TIFF, Pillow's own error is the one worth reporting.
        if not _is_tiff_file(path):
            raise
        try:
            import tifffile
        except ImportError:
            raise OSError(f"{exc} (Pillow cannot read this kind of TIFF; installing the "
                          "optional 'tifffile' package may allow it to be opened)") from exc
        with tifffile.TiffFile(path) as tif:
            arr = tif.asarray()
            samples = tif.pages[0].samplesperpixel
        if arr.ndim > 2 and arr.shape[0] not in (3, 4) and arr.shape[-1] not in (3, 4):
            arr = arr[0]  # multi-frame stack -> first frame
        elif arr.ndim == 3 and arr.shape[0] in (3, 4) and arr.shape[-1] not in (3, 4):
            if samples == arr.shape[0]:
                arr = np.moveaxis(arr, 0, -1)   # planar RGB(A) -> rows x columns x samples
            else:
                arr = arr[0]                    # a stack of exactly 3 or 4 frames
        disp = normalize_for_display(arr)
        return _display_image(Image.fromarray(disp))

    # Multi-frame (e.g. TIFF stack): use the first frame.
    if getattr(img, "n_frames", 1) > 1:
        img.seek(0)

    if img.mode in ("I", "I;16", "I;16B", "I;16L", "F"):
        arr = np.asarray(img)
        disp = normalize_for_display(arr)
        return _display_image(Image.fromarray(disp))

    # Even 8-bit L / RGB images may occasionally be low-contrast; leave them as
    # captured (no stretch) so the displayed image matches the file.
    return _display_image(img)


# ===========================================================================
# Reading pixel-size / scale calibration from image metadata
# ===========================================================================
_UNIT_TO_NM = {
    "nm": 1.0, "nanometer": 1.0, "nanometers": 1.0, "nanometre": 1.0, "nanometres": 1.0,
    "um": 1e3, "µm": 1e3, "μm": 1e3, "micron": 1e3, "microns": 1e3,
    "micrometer": 1e3, "micrometers": 1e3, "micrometre": 1e3, "micrometres": 1e3,
    "mm": 1e6, "millimeter": 1e6, "millimetre": 1e6,
    "cm": 1e7, "centimeter": 1e7, "centimetre": 1e7,
    "m": 1e9, "meter": 1e9, "metre": 1e9,
    "a": 0.1, "å": 0.1, "angstrom": 0.1, "angstroms": 0.1, "ang": 0.1,
    "pm": 1e-3, "picometer": 1e-3,
    "in": 2.54e7, "inch": 2.54e7, "inches": 2.54e7,
    "\\u00b5m": 1e3,   # ImageJ writes µ as a literal escape sequence
}

# A decimal number, optionally in scientific notation.
_NUM = r'[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?'

# Plain TIFF resolution tags coarser than this (10 µm/px, i.e. below ~2540 dpi)
# are assumed to be a print/display resolution, not a microscope calibration.
MAX_TRUSTED_RESOLUTION_NM_PER_PX = 1e4


def _unit_to_nm(unit: Optional[str]) -> Optional[float]:
    """Return the length of one ``unit`` in nanometres, or None if unrecognised."""
    if not unit:
        return None
    return _UNIT_TO_NM.get(str(unit).strip().lower())


def _rational_to_float(v) -> Optional[float]:
    """Coerce a TIFF RATIONAL (tuple, Fraction, or number) to a float."""
    try:
        if isinstance(v, (tuple, list)) and len(v) == 2:
            num, den = v
            return (num / den) if den else None
        return float(v)
    except Exception:
        return None


def _ome_size_nm(ome: str, axis: str) -> Optional[float]:
    m = re.search(r'PhysicalSize' + axis + r'="(' + _NUM + r')"', ome)
    if not m:
        return None
    mu = re.search(r'PhysicalSize' + axis + r'Unit="([^"]+)"', ome)
    nm = _unit_to_nm(mu.group(1)) if mu else _unit_to_nm("µm")  # OME default is µm
    return float(m.group(1)) * nm if nm else None


def _pixsize_from_ome(ome: str) -> Optional[Tuple[float, str]]:
    x_nm = _ome_size_nm(ome, "X")
    return (x_nm, "OME-XML metadata") if x_nm else None


def _pixsize_from_zeiss(sem: dict) -> Optional[Tuple[float, str]]:
    for key in ("ap_image_pixel_size", "ap_pixel_size", "dp_pixel_size", "ap_pixelsize"):
        if key not in sem:
            continue
        v = sem[key]
        parts = list(v) if isinstance(v, (tuple, list)) else [v]
        number = next((float(x) for x in parts
                       if isinstance(x, (int, float)) and not isinstance(x, bool)), None)
        unit_nm = None
        for s in parts:
            if isinstance(s, str):
                u = _unit_to_nm(s.strip())
                if u:
                    unit_nm = u
                    break
                m = re.search(r'([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)\s*([a-zAµμÅ]+)', s)
                if m:
                    n = _unit_to_nm(m.group(2))
                    if n:
                        return (float(m.group(1)) * n, "Zeiss SEM metadata")
        if number is not None:
            return (number * (unit_nm or 1.0), "Zeiss SEM metadata")
    return None


_SENSOR_WORDS = re.compile(r'(?:detector|sensor|camera|ccd|cmos)\W*$', re.IGNORECASE)


def _pixsize_from_description(text: str) -> Optional[Tuple[float, str]]:
    """Find an explicit pixel size ("pixel size = 2.5 nm", "2.5 nm/px") in text.

    A bare "scale = 500 nm" is deliberately not matched: it is usually the
    length of the scale bar, not the size of a pixel.
    """
    patterns = [
        r'pixel\s*(?:size|width|spacing)\s*[=:]?\s*(' + _NUM + r')\s*(nm|µm|μm|um|mm|cm|angstroms?|å|pm)\b',
        # The look-behinds keep "0,25 nm/pixel" from being read as 25 nm/pixel.
        r'(?<![\d.])(?<!\d,)(' + _NUM + r')\s*(nm|µm|μm|um|mm)\s*(?:/|per)\s*(?:pixel|px)\b',
    ]
    for pat in patterns:
        for m in re.finditer(pat, text, re.IGNORECASE):
            nm = _unit_to_nm(m.group(2))
            if not nm:
                continue
            value = float(m.group(1)) * nm
            # "Detector pixel size = 14 um" is the pitch of the camera sensor,
            # not the size of an image pixel.  A sensor pitch is several µm, so
            # a smaller value is an image pixel even when the text just before
            # it names the camera ("Camera = CCD; Pixel size = 0.5 nm").
            if value >= 1000.0 and _SENSOR_WORDS.search(text[max(0, m.start() - 24):m.start()]):
                continue
            return (value, "image description")
    return None


def _not_square(x: Optional[float], y: Optional[float], tol: float = 0.01) -> bool:
    """True if both pixel dimensions are known and differ by more than ``tol``."""
    return bool(x and y and x > 0 and y > 0 and abs(x / y - 1.0) > tol)


def read_pixel_size_nm(path: str) -> Optional[Tuple[float, str, bool]]:
    """Best-effort read of the pixel size from a TIFF's metadata.

    Returns ``(nm_per_px, source, trusted)`` or None if nothing was found.
    ``trusted`` is False for print/display resolutions (e.g. 72 dpi), which are
    usually NOT a real microscope calibration, and for pixels that are not
    square, which this tool cannot measure correctly.  Handles FEI/Thermo,
    Zeiss SEM, ImageJ, OME-XML, description text, and TIFF resolution tags.
    """
    ext = os.path.splitext(path)[1].lower()
    if ext not in (".tif", ".tiff"):
        return None
    try:
        import tifffile
    except ImportError:
        return None
    try:
        with tifffile.TiffFile(path) as tif:
            page = tif.pages[0]

            # 1) FEI / Thermo Fisher (PixelWidth is in metres).
            try:
                fei = tif.fei_metadata
                if fei:
                    for section in ("Scan", "EScan", "IonScan"):
                        sec = fei.get(section)
                        if isinstance(sec, dict) and sec.get("PixelWidth"):
                            pw = float(sec["PixelWidth"])
                            if pw > 0:
                                return (pw * 1e9, f"FEI/Thermo metadata ({section})", True)
            except Exception:
                pass

            # 2) Zeiss SEM.
            try:
                sem = tif.sem_metadata
                if sem:
                    res = _pixsize_from_zeiss(sem)
                    if res and res[0] > 0:
                        return (res[0], res[1], True)
            except Exception:
                pass

            def resolution(tag_name):
                tag = page.tags.get(tag_name)
                return _rational_to_float(tag.value) if tag is not None else None

            # 3) ImageJ (unit in metadata, spacing from XResolution).
            try:
                ij = tif.imagej_metadata
                if ij and ij.get("unit"):
                    nm_unit = _unit_to_nm(ij.get("unit"))
                    ppu = resolution("XResolution")
                    if nm_unit and ppu:
                        nm_per_px = (1.0 / ppu) * nm_unit
                        # ImageJ labels a plain print resolution "inch" (or
                        # "cm"), so those get the same doubt as in step 6.
                        unit = str(ij.get("unit")).strip().lower()
                        trusted = not (
                            unit in ("in", "inch", "inches")
                            or (unit in ("cm", "centimeter", "centimetre")
                                and nm_per_px > MAX_TRUSTED_RESOLUTION_NM_PER_PX)
                            or _not_square(ppu, resolution("YResolution")))
                        return (nm_per_px, "ImageJ metadata", trusted)
            except Exception:
                pass

            # 4) OME-XML.
            try:
                ome = tif.ome_metadata
                if ome:
                    res = _pixsize_from_ome(ome)
                    if res and res[0] > 0:
                        trusted = not _not_square(res[0], _ome_size_nm(ome, "Y"))
                        return (res[0], res[1], trusted)
            except Exception:
                pass

            # 5) Description text (e.g. "Pixel size = 2.5 nm").
            try:
                desc = page.tags.get("ImageDescription")
                if desc is not None:
                    res = _pixsize_from_description(str(desc.value))
                    if res and res[0] > 0:
                        return (res[0], res[1], True)
            except Exception:
                pass

            # 6) Plain TIFF resolution tags (least reliable: often print DPI).
            try:
                ru = page.tags.get("ResolutionUnit")
                ppu = resolution("XResolution")
                if ppu and ppu > 0:
                    unit_code = int(ru.value) if ru is not None else 2
                    nm_per_unit = {3: 1e7, 2: 2.54e7}.get(unit_code)  # cm / inch
                    if nm_per_unit:
                        # cm-based resolution is usually a real calibration; inch
                        # (dpi) usually is not, and neither is a cm-based value
                        # in the print-resolution range, so flag those untrusted.
                        nm_per_px = (1.0 / ppu) * nm_per_unit
                        trusted = (unit_code == 3
                                   and nm_per_px <= MAX_TRUSTED_RESOLUTION_NM_PER_PX
                                   and not _not_square(ppu, resolution("YResolution")))
                        return (nm_per_px, "TIFF resolution tag", trusted)
            except Exception:
                pass
    except Exception:
        return None
    return None


def detect_scale_bar(gray: "np.ndarray") -> Optional[Tuple[int, int, int, int, float]]:
    """Locate a burned-in scale bar in a grayscale image.

    A scale bar is a solid, high-contrast (near-white or near-black), long and
    thin horizontal rectangle — usually in a corner or an info band.  Returns
    ``(x0, y0, x1, y1, length_px)`` for the best candidate's bounding box, or
    None if nothing bar-like is found.  Heuristic: the caller highlights the
    result so the user can confirm it before calibrating.
    """
    try:
        from scipy import ndimage
    except ImportError:
        return None
    a = np.asarray(gray)
    if a.ndim != 2 or a.size == 0:
        return None
    h, w = a.shape
    # Best candidate at least 3 px high, and best thinner one.  A 1-2 px line
    # (a rule in the data band, a saturated streak) has a huge aspect ratio and
    # would otherwise outrank a real bar, so it is only used as a last resort.
    best = {True: None, False: None}
    best_score = {True: -1.0, False: -1.0}
    for mask in (a >= 240, a <= 15):        # bright bar or dark bar
        if not mask.any():
            continue
        labeled, n = ndimage.label(mask)
        if n == 0:
            continue
        slices = ndimage.find_objects(labeled)
        areas = np.bincount(labeled.ravel())
        for idx, sl in enumerate(slices, start=1):
            if sl is None:
                continue
            ys, xs = sl
            x0, x1 = int(xs.start), int(xs.stop)
            y0, y1 = int(ys.start), int(ys.stop)
            bw, bh = x1 - x0, y1 - y0
            if bw <= 0 or bh <= 0:
                continue
            aspect = bw / bh
            solidity = areas[idx] / (bw * bh)
            if (0.03 * w <= bw <= 0.95 * w        # spans a plausible fraction of width
                    and bh <= max(0.02 * h, 25)   # thin
                    and aspect >= 6.0             # long and thin
                    and solidity >= 0.7):         # a solid filled rectangle
                # Prefer long, wide bars lower in the image (where bars usually sit).
                score = aspect + (bw / w) * 20.0 + (y0 / h) * 8.0
                thick = bh >= 3
                if score > best_score[thick]:
                    best_score[thick] = score
                    # Stray near-saturated pixels touching the ends of the bar
                    # widen the bounding box, so the length is the span of the
                    # columns in which the object is more than a quarter as
                    # thick as it typically (median) is, and more than 1 px
                    # (rounded or JPEG-softened ends still belong to the bar).
                    cols = (labeled[sl] == idx).sum(axis=0)
                    med = float(np.median(cols))
                    keep = np.flatnonzero(cols > min(med / 2.0, max(1.0, med / 4.0)))
                    length = float(keep[-1] - keep[0] + 1)
                    best[thick] = (x0, y0, x1, y1, length)
    return best[True] or best[False]


def resolve_session_image(saved_path: str, session_dir: str, max_up: int = 4) -> Optional[str]:
    """Locate an image referenced by a session file.

    Sessions store absolute paths, which break when the project folder is moved
    or copied to another machine.  If ``saved_path`` no longer exists, the
    trailing part of it is looked up relative to the session file's folder and
    up to ``max_up`` of its parent folders, longest match first, so a folder
    tree that moved as a whole resolves to the same relative location.  A bare
    file name is only looked up next to the session file.  The saved path may
    use Windows or POSIX separators regardless of the current platform.
    """
    if saved_path and os.path.isfile(saved_path):
        return saved_path
    parts = [p for p in re.split(r"[\\/]", saved_path or "") if p]
    if parts and re.fullmatch(r"[A-Za-z]:", parts[0]):
        parts = parts[1:]                      # drop a Windows drive letter
    if not parts or any(p in (".", "..") for p in parts):
        return None
    bases = [os.path.abspath(session_dir)]
    for _ in range(max_up):
        parent = os.path.dirname(bases[-1])
        if not parent or parent == bases[-1]:
            break
        bases.append(parent)
    for n in range(len(parts), 0, -1):
        for base in (bases if n > 1 else bases[:1]):
            candidate = os.path.join(base, *parts[-n:])
            if os.path.isfile(candidate):
                return candidate
    return None


def scale_is_outlier(nm_per_px: float, other_nm_per_px: Sequence[float],
                     factor: float = 100.0) -> bool:
    """True if a pixel size differs from the median of the others by ``factor``.

    A scale bar entered with the wrong unit (0.5 nm instead of 0.5 µm) is off by
    1000x, which makes areas wrong by 10^6; the GUI uses this to ask the user to
    double-check such a value.
    """
    others = sorted(v for v in other_nm_per_px if v and v > 0)
    if not others or not nm_per_px or nm_per_px <= 0:
        return False
    mid = len(others) // 2
    median = others[mid] if len(others) % 2 else (others[mid - 1] + others[mid]) / 2.0
    ratio = nm_per_px / median
    return ratio >= factor or ratio <= 1.0 / factor


def parse_positive_float(text: str) -> Optional[float]:
    """Parse a user-typed positive number; accepts a decimal comma.

    Text that looks like thousands grouping ("2,000") is rejected rather than
    read as 2.0: guessing wrong would make the scale 1000x off.
    """
    s = str(text).strip()
    if re.fullmatch(r"[1-9]\d{0,2}(?:,\d{3})+(?:\.\d*)?", s):
        return None
    try:
        value = float(s.replace(",", "."))
    except ValueError:
        return None
    return value if value > 0 and math.isfinite(value) else None


NUMBER_HINT = "Enter a positive number without a thousands separator, e.g. 2000 or 0.5."


def _listing(items: Sequence[str], limit: int = 15) -> str:
    """Join lines for a message box, cutting a long list short."""
    lines = list(items[:limit])
    if len(items) > limit:
        lines.append(f"… and {len(items) - limit} more")
    return "\n".join(lines)


def _module_available(name: str) -> bool:
    """True if an optional package can be imported.

    Importing it is the only reliable test: a package can be present and still
    fail to load (for example SciPy built for another NumPy version).
    """
    import importlib
    try:
        importlib.import_module(name)
    except Exception:
        return False
    return True


def display_names(paths: Sequence[str]) -> List[str]:
    """Short labels for image paths: the file name, made unique where needed.

    Images that share a file name (``control/001.tif`` and ``treated/001.tif``)
    get as many parent folders as it takes to tell them apart.  Entries with
    one and the same path cannot be told apart and keep the same label.
    """
    split = [[p for p in re.split(r"[\\/]", path) if p] or [path] for path in paths]
    depth = [1] * len(split)
    while True:
        labels = ["/".join(parts[-d:]) for parts, d in zip(split, depth)]
        groups = {}
        for i, label in enumerate(labels):
            groups.setdefault(label, []).append(i)
        grew = False
        for members in groups.values():
            # More folders only help where the paths differ; for one path
            # listed twice they would just grow into the whole path.
            if len({tuple(split[i]) for i in members}) > 1:
                for i in members:
                    if depth[i] < len(split[i]):
                        depth[i] += 1
                        grew = True
        if not grew:
            return labels


# ===========================================================================
# GUI
# ===========================================================================
def _build_gui():
    """Import tkinter lazily and build the application class.

    Kept inside a function so the module can be imported (to test the pure
    functions above) on systems without a display.
    """
    import tkinter as tk
    from tkinter import ttk, filedialog, messagebox
    from PIL import ImageTk

    # Distinct outline colours, cycled by cell number.
    PALETTE = [
        "#ff3b30", "#34c759", "#0a84ff", "#ffcc00", "#ff2d55",
        "#5ac8fa", "#af52de", "#ff9500", "#30d158", "#64d2ff",
        "#ff6482", "#a2845e",
    ]

    def hex_to_rgb(h):
        return tuple(int(h[i:i + 2], 16) for i in (1, 3, 5))

    def _make_modal(dialog):
        """Grab input for ``dialog`` once it is on screen.

        X11 refuses a grab on a window that is not yet mapped, so wait for it
        first.  If the grab still cannot be taken the dialog simply stays
        non-modal rather than failing to open.
        """
        try:
            dialog.wait_visibility()
            dialog.grab_set()
        except tk.TclError:
            pass

    class ScaleDialog(tk.Toplevel):
        """Modal dialog: enter a physical length + unit for scale-bar calibration."""

        def __init__(self, parent, pixel_length, prompt):
            super().__init__(parent)
            self.title("Set scale")
            self.resizable(False, False)
            self.result = None
            self.transient(parent)

            frm = ttk.Frame(self, padding=12)
            frm.grid(sticky="nsew")
            ttk.Label(frm, text=prompt, wraplength=340, justify="left").grid(
                row=0, column=0, columnspan=2, sticky="w", pady=(0, 8))
            ttk.Label(frm, text=f"Traced length: {pixel_length:.1f} px").grid(
                row=1, column=0, columnspan=2, sticky="w", pady=(0, 8))
            ttk.Label(frm, text="Known length:").grid(row=2, column=0, sticky="w")
            self.len_var = tk.StringVar()
            len_entry = ttk.Entry(frm, textvariable=self.len_var, width=14)
            len_entry.grid(row=2, column=1, sticky="w", padx=(6, 0))
            ttk.Label(frm, text="Unit:").grid(row=3, column=0, sticky="w", pady=(6, 0))
            self.unit_var = tk.StringVar(value="nm")
            ttk.Combobox(frm, textvariable=self.unit_var, width=11, state="readonly",
                         values=["nm", "µm", "mm"]).grid(
                row=3, column=1, sticky="w", padx=(6, 0), pady=(6, 0))
            ttk.Label(frm, text="Areas are reported in nm².", foreground="#666").grid(
                row=4, column=0, columnspan=2, sticky="w", pady=(8, 0))
            btns = ttk.Frame(frm)
            btns.grid(row=5, column=0, columnspan=2, pady=(12, 0), sticky="e")
            ttk.Button(btns, text="Cancel", command=self._cancel).grid(row=0, column=0, padx=4)
            ttk.Button(btns, text="OK", command=self._ok).grid(row=0, column=1)

            self.pixel_length = pixel_length
            len_entry.focus_set()
            self.bind("<Return>", lambda e: self._ok())
            self.bind("<Escape>", lambda e: self._cancel())
            self.protocol("WM_DELETE_WINDOW", self._cancel)
            _make_modal(self)

        def _ok(self):
            length = parse_positive_float(self.len_var.get())
            if length is None:
                messagebox.showerror("Invalid length", NUMBER_HINT, parent=self)
                return
            # Normalise to nm/px so areas come out in nm^2 regardless of the
            # unit the scale bar was labelled in.
            length_nm = to_nm(length, self.unit_var.get().strip())
            self.result = (length_nm / self.pixel_length, "nm")
            self.destroy()

        def _cancel(self):
            self.result = None
            self.destroy()

    class PixelSizeDialog(tk.Toplevel):
        """Modal dialog: enter pixel size directly (physical units per pixel)."""

        def __init__(self, parent, current):
            super().__init__(parent)
            self.title("Enter pixel size")
            self.resizable(False, False)
            self.result = None
            self.transient(parent)

            frm = ttk.Frame(self, padding=12)
            frm.grid(sticky="nsew")
            ttk.Label(frm, text="Physical size of one pixel.\n"
                      "e.g. 2.5 nm means each pixel is 2.5 nm across.",
                      justify="left").grid(row=0, column=0, columnspan=2, sticky="w", pady=(0, 8))
            ttk.Label(frm, text="Pixel size:").grid(row=1, column=0, sticky="w")
            # Show the current pixel size back in nm (calibration is stored in nm/px).
            self._prefill = "" if not current.calibrated else f"{current.unit_per_px:g}"
            self.val_var = tk.StringVar(value=self._prefill)
            entry = ttk.Entry(frm, textvariable=self.val_var, width=14)
            entry.grid(row=1, column=1, sticky="w", padx=(6, 0))
            ttk.Label(frm, text="Unit:").grid(row=2, column=0, sticky="w", pady=(6, 0))
            self.unit_var = tk.StringVar(value="nm")
            ttk.Combobox(frm, textvariable=self.unit_var, width=11, state="readonly",
                         values=["nm", "µm", "mm"]).grid(
                row=2, column=1, sticky="w", padx=(6, 0), pady=(6, 0))
            ttk.Label(frm, text="Areas are reported in nm².", foreground="#666").grid(
                row=3, column=0, columnspan=2, sticky="w", pady=(8, 0))
            btns = ttk.Frame(frm)
            btns.grid(row=4, column=0, columnspan=2, pady=(12, 0), sticky="e")
            ttk.Button(btns, text="Cancel", command=self._cancel).grid(row=0, column=0, padx=4)
            ttk.Button(btns, text="OK", command=self._ok).grid(row=0, column=1)

            entry.focus_set()
            self.bind("<Return>", lambda e: self._ok())
            self.bind("<Escape>", lambda e: self._cancel())
            self.protocol("WM_DELETE_WINDOW", self._cancel)
            _make_modal(self)

        def _ok(self):
            unit = self.unit_var.get().strip()
            if self._prefill and unit == "nm" and self.val_var.get().strip() == self._prefill:
                # Nothing was edited.  The prefill is rounded for display, so
                # applying it would replace the stored scale with a rounded one.
                self._cancel()
                return
            val = parse_positive_float(self.val_var.get())
            if val is None:
                messagebox.showerror("Invalid value", NUMBER_HINT, parent=self)
                return
            self.result = (to_nm(val, unit), "nm")
            self.destroy()

        def _cancel(self):
            self.result = None
            self.destroy()

    class App(tk.Tk):
        MODE_IDLE = "idle"
        MODE_CALIBRATE = "calibrate"
        STYLE_FREEHAND = "freehand"
        STYLE_POLYGON = "polygon"

        CLOSE_THRESHOLD_PX = 12   # snap-close a polygon near its first vertex
        FREEHAND_SAMPLE_PX = 2.5  # min canvas travel between freehand samples
        RDP_EPSILON_PX = 1.5      # freehand simplification tolerance (image px)
        MIN_TRACE_PX = 4.0        # ignore freehand traces smaller than this (image px)

        def __init__(self, initial_paths=None):
            super().__init__()
            self.title("TEM Cell Area Tool")
            self.geometry("1280x820")
            self.minsize(900, 600)

            # Tk on macOS ("aqua") differs from Windows/X11 in ways handled below.
            self._aqua = self.tk.call("tk", "windowingsystem") == "aqua"

            # --- session state ----------------------------------------------
            self.docs: List[ImageDoc] = []
            self.active: Optional[int] = None

            self.draw_style = self.STYLE_FREEHAND
            self.mode = self.MODE_IDLE
            self.current_points: List[Point] = []   # in-progress trace (image coords)
            self._freehanding = False
            self.calib_points: List[Point] = []

            self.pan_tool = False
            self._panning = False
            self._pan_start = (0, 0)
            self._pan_offset_start = (0.0, 0.0)

            self._cursor_img: Optional[Point] = None
            # Highlighted scale bar: (x0, y0, x1, y1, length_px).
            self._detected_bar: Optional[Tuple[int, int, int, int, float]] = None
            self._selected: Optional[Tuple[int, int]] = None  # (doc_index, cell_index)
            self._row_map: List[Tuple[int, int]] = []
            self._suppress_tree_event = False   # guard against re-entrant tree selection
            self._tk_photo = None
            self._photo_src = None       # the image self._tk_photo was scaled from ...
            self._photo_key = None       # ... and the crop box, size and filter used
            self._labels: List[str] = []       # display name of each open image
            self._label_paths: List[str] = []  # the paths self._labels was built from
            self._unsaved = False        # outlines / scales changed since the last session save
            self._not_loaded: List[str] = []   # images the loaded session could not open

            self._build_widgets()
            self._bind_events()
            self.protocol("WM_DELETE_WINDOW", self._on_close)
            if self._aqua:
                try:    # Cmd+Q does not go through WM_DELETE_WINDOW on macOS
                    self.createcommand("::tk::mac::Quit", self._on_close)
                except tk.TclError:
                    pass
            # Open wide enough for the toolbar, whatever the platform's fonts,
            # and no taller than a small screen (less the task bar and title).
            self.update_idletasks()
            width = min(max(1280, self.winfo_reqwidth()), self.winfo_screenwidth())
            height = max(600, min(820, self.winfo_screenheight() - 80))
            self.geometry(f"{width}x{height}")
            if height < 820 and sys.platform == "win32":
                # Windows does not keep a new window clear of the task bar, which
                # would cover the status bar on a screen this low: fill the screen.
                self.state("zoomed")
            self._refresh_tree()
            self._update_status()

            if initial_paths:
                self.after(120, lambda: [self.add_image(p) for p in initial_paths])

        def report_callback_exception(self, exc, val, tb):
            """Show errors raised inside event handlers instead of losing them.

            Tk's default prints to stderr, which does not exist when the tool is
            launched without a console (pythonw / double-click).
            """
            import traceback
            text = "".join(traceback.format_exception(exc, val, tb))
            try:
                sys.stderr.write(text)
            except Exception:
                pass
            try:
                messagebox.showerror("Unexpected error", text[-1500:], parent=self)
            except Exception:
                pass

        def _on_close(self):
            """Close the window, asking first if outlines would be lost."""
            try:
                n = self.total_cells()
                if self._unsaved and n and not messagebox.askyesno(
                        "Quit",
                        f"{n} outlined cell(s) have not been saved to a session file "
                        "since they were last changed.\n\nQuit without saving a session?",
                        parent=self):
                    return
            except Exception:
                pass    # never let a failed question make the window impossible to close
            self.destroy()

        REPORT_BOX_LINES = 20     # longer messages go into a scrolling window

        def _show_report(self, title: str, text: str):
            """Show ``text`` in a message box, or in a scrolling window if it is long.

            A message box cannot scroll, so a long list in one is either cut
            short or runs off the screen.
            """
            lines = text.split("\n")
            if len(lines) <= self.REPORT_BOX_LINES:
                messagebox.showinfo(title, text, parent=self)
                return
            win = tk.Toplevel(self)
            win.title(title)
            win.transient(self)
            frm = ttk.Frame(win, padding=10)
            frm.pack(fill="both", expand=True)
            close = ttk.Button(frm, text="Close", command=win.destroy)
            close.pack(side="bottom", anchor="e", pady=(8, 0))
            box = tk.Text(frm, wrap="word", font="TkDefaultFont", relief="flat",
                          width=min(110, max(60, max(len(s) for s in lines) + 2)),
                          height=min(28, len(lines) + 1))
            scroll = ttk.Scrollbar(frm, orient="vertical", command=box.yview)
            box.configure(yscrollcommand=scroll.set)
            scroll.pack(side="right", fill="y")
            box.pack(side="left", fill="both", expand=True)
            box.insert("1.0", text)
            box.configure(state="disabled")
            # A read-only Text does not take the focus on every platform, and
            # without the focus its text cannot be selected and copied.
            box.bind("<Button-1>", lambda e: box.focus_set())
            for key in ("<Return>", "<KP_Enter>", "<Escape>"):
                win.bind(key, lambda e: win.destroy())
            win.update_idletasks()
            x = self.winfo_rootx() + max(0, (self.winfo_width() - win.winfo_reqwidth()) // 2)
            y = self.winfo_rooty() + max(0, (self.winfo_height() - win.winfo_reqheight()) // 3)
            win.geometry(f"+{x}+{y}")
            close.focus_set()
            _make_modal(win)
            self.wait_window(win)

        # ------------------------------------------------------------------
        # Active-document helpers
        # ------------------------------------------------------------------
        @property
        def doc(self) -> Optional[ImageDoc]:
            if self.active is None or not (0 <= self.active < len(self.docs)):
                return None
            return self.docs[self.active]

        def global_number(self, doc_index: int, cell_index: int) -> int:
            n = sum(len(self.docs[i].cells) for i in range(doc_index))
            return n + cell_index + 1

        def total_cells(self) -> int:
            return sum(len(d.cells) for d in self.docs)

        # ------------------------------------------------------------------
        # Widgets
        # ------------------------------------------------------------------
        def _build_widgets(self):
            bar1 = ttk.Frame(self, padding=(6, 4))
            bar1.pack(side="top", fill="x")
            ttk.Button(bar1, text="Open image(s)…", command=self.open_dialog).pack(side="left", padx=2)
            ttk.Label(bar1, text="Image:").pack(side="left", padx=(10, 2))
            self.image_selector = ttk.Combobox(bar1, state="readonly", width=42, values=[])
            self.image_selector.pack(side="left", padx=2)
            self.image_selector.bind("<<ComboboxSelected>>", self._on_select_image)
            ttk.Button(bar1, text="Close image", command=self.close_image).pack(side="left", padx=6)
            # These three sit in the first row so that they stay on screen with
            # the wider fonts and buttons of Linux and macOS.
            ttk.Separator(bar1, orient="vertical").pack(side="left", fill="y", padx=6)
            ttk.Button(bar1, text="Export table + JPEGs…", command=self.export_all).pack(side="left", padx=2)
            ttk.Button(bar1, text="Save session", command=self.save_session).pack(side="left", padx=2)
            ttk.Button(bar1, text="Load session", command=self.load_session).pack(side="left", padx=2)

            bar2 = ttk.Frame(self, padding=(6, 0))
            bar2.pack(side="top", fill="x")

            self.style_btn = ttk.Button(bar2, text="Draw: Freehand", command=self.toggle_draw_style)
            self.style_btn.pack(side="left", padx=2)
            ttk.Separator(bar2, orient="vertical").pack(side="left", fill="y", padx=6)
            ttk.Button(bar2, text="Scale from file", command=self.detect_scale_from_file).pack(side="left", padx=2)
            ttk.Button(bar2, text="Detect scale bar", command=self.autodetect_scale_bar).pack(side="left", padx=2)
            ttk.Button(bar2, text="Set scale (bar)", command=self.start_calibration).pack(side="left", padx=2)
            ttk.Button(bar2, text="Pixel size…", command=self.enter_pixel_size).pack(side="left", padx=2)
            ttk.Separator(bar2, orient="vertical").pack(side="left", fill="y", padx=6)
            ttk.Button(bar2, text="Finish (Enter)", command=self.finish_polygon).pack(side="left", padx=2)
            ttk.Button(bar2, text="Undo point", command=self.undo_point).pack(side="left", padx=2)
            ttk.Button(bar2, text="Delete cell", command=self.delete_selected_cell).pack(side="left", padx=2)
            ttk.Separator(bar2, orient="vertical").pack(side="left", fill="y", padx=6)
            self.pan_btn = ttk.Button(bar2, text="Pan tool: off", command=self.toggle_pan_tool)
            self.pan_btn.pack(side="left", padx=2)
            ttk.Button(bar2, text="Fit", command=self.fit_to_window).pack(side="left", padx=2)
            ttk.Button(bar2, text="+", width=3, command=lambda: self.zoom_center(1.25)).pack(side="left", padx=1)
            ttk.Button(bar2, text="−", width=3, command=lambda: self.zoom_center(0.8)).pack(side="left", padx=1)

            body = ttk.Frame(self)
            body.pack(side="top", fill="both", expand=True)

            self.canvas = tk.Canvas(body, bg="#141414", highlightthickness=0, cursor="crosshair")
            self.canvas.pack(side="left", fill="both", expand=True)

            # Wide enough for the four table columns below (320 px in all).
            side = ttk.Frame(body, padding=(8, 6), width=340)
            side.pack(side="right", fill="y")
            side.pack_propagate(False)

            ttk.Label(side, text="Cells (all images)", font=("", 10, "bold")).pack(anchor="w")
            cols = ("no", "image", "area", "perim")
            # The table asks for a few rows only and takes whatever height is
            # left, so that the summary and the help text below it stay visible.
            self.tree = ttk.Treeview(side, columns=cols, show="headings", height=4, selectmode="browse")
            for c, txt, w, anchor in (
                ("no", "No.", 40, "center"),
                ("image", "Image", 110, "w"),
                ("area", "Area", 96, "e"),
                ("perim", "Perim.", 74, "e"),
            ):
                self.tree.heading(c, text=txt)
                self.tree.column(c, width=w, anchor=anchor)
            self.tree.pack(fill="both", expand=True, pady=(4, 6))
            self.tree.bind("<<TreeviewSelect>>", self._on_tree_select)
            self.tree.bind("<Delete>", lambda e: self.delete_selected_cell())
            if self._aqua:
                self.tree.bind("<BackSpace>", self._on_tree_backspace)

            self.summary_var = tk.StringVar()
            ttk.Label(side, textvariable=self.summary_var, justify="left").pack(anchor="w")

            help_text = (
                "How to use\n"
                "1. Open images. TIFF pixel-size metadata is\n"
                "   read and applied automatically.\n"
                "2. Otherwise set the scale for EACH image:\n"
                "   • 'Detect scale bar' finds the bar,\n"
                "     highlights it (magenta) to verify, then\n"
                "     asks its length.\n"
                "   • 'Set scale (bar)' to trace it yourself.\n"
                "   • 'Pixel size…' to type it in.\n"
                "3. Outline cells:\n"
                "   • Freehand — press and drag around\n"
                "     the cell, release to close.\n"
                "   • Points — click vertices; click the\n"
                "     first point / Enter to close.\n"
                "4. Each cell is numbered and listed with\n"
                "   its area in nm² (from the scale bar).\n"
                "   Switch images from the dropdown;\n"
                "   outlines are kept per image.\n"
                "5. Export writes the table (CSV) and one\n"
                "   annotated JPEG per image.\n\n"
                "Wheel = zoom · middle-drag = pan\n"
                "Backspace = undo point · Esc = cancel"
            )
            ttk.Label(side, text=help_text, justify="left", foreground="#666").pack(anchor="w", pady=(10, 0))

            self.status_var = tk.StringVar()
            ttk.Label(self, textvariable=self.status_var, relief="sunken", anchor="w",
                      padding=(6, 2)).pack(side="bottom", fill="x")

        def _bind_events(self):
            self.canvas.bind("<Configure>", self._on_canvas_configure)
            self.canvas.bind("<Button-1>", self.on_left_down)
            self.canvas.bind("<B1-Motion>", self.on_left_drag)
            self.canvas.bind("<ButtonRelease-1>", self.on_left_up)
            self.canvas.bind("<Double-Button-1>", self.on_double_click)
            # Tk 8.6 on macOS numbers the right button 2 and the middle button 3,
            # the reverse of Windows/X11 (and of Tk 8.7+ on macOS).
            right, middle = (2, 3) if (self._aqua and tk.TkVersion < 8.7) else (3, 2)
            self.canvas.bind(f"<Button-{right}>", self.on_right_click)
            if self._aqua:
                self.canvas.bind("<Control-Button-1>", self.on_right_click)
            self.canvas.bind("<Motion>", self.on_motion)
            self.canvas.bind(f"<Button-{middle}>", self.on_pan_start)
            self.canvas.bind(f"<B{middle}-Motion>", self.on_pan_move)
            self.canvas.bind(f"<ButtonRelease-{middle}>", self.on_pan_end)
            self.canvas.bind("<MouseWheel>", self.on_wheel)
            if self.tk.call("tk", "windowingsystem") == "x11" and tk.TkVersion < 8.7:
                # Only X11 with Tk 8.6 reports the wheel as buttons 4 and 5;
                # elsewhere those are the mouse's back / forward buttons.
                self.canvas.bind("<Button-4>", lambda e: self.on_wheel(e, delta=120))
                self.canvas.bind("<Button-5>", lambda e: self.on_wheel(e, delta=-120))

            self.bind("<Return>", lambda e: self.finish_polygon())
            self.bind("<KP_Enter>", lambda e: self.finish_polygon())
            self.bind("<BackSpace>", lambda e: self.undo_point())
            self.bind("<Escape>", lambda e: self.cancel_current())
            self.bind("<f>", lambda e: self.fit_to_window())
            self.bind("<plus>", lambda e: self.zoom_center(1.25))
            self.bind("<equal>", lambda e: self.zoom_center(1.25))
            self.bind("<minus>", lambda e: self.zoom_center(0.8))

        def _on_tree_backspace(self, event):
            """macOS: the "delete" key is BackSpace, which is also "undo point".

            With the table focused it deletes the selected cell, unless an
            outline or a scale bar is being drawn, when it undoes a point.
            """
            if self.current_points or self.mode == self.MODE_CALIBRATE:
                self.undo_point()
            else:
                self.delete_selected_cell()
            return "break"   # keep the window-level undo binding from also firing

        # ------------------------------------------------------------------
        # Opening / switching / closing images
        # ------------------------------------------------------------------
        def open_dialog(self):
            images = "*.tif *.tiff *.png *.jpg *.jpeg *.bmp *.gif"
            tiff = "*.tif *.tiff"
            if self.tk.call("tk", "windowingsystem") == "x11":
                # Tk's own file dialog on X11 matches patterns case-sensitively.
                images += " " + images.upper()
                tiff += " " + tiff.upper()
            paths = filedialog.askopenfilenames(
                title="Open image(s)",
                filetypes=[("Images", images), ("TIFF", tiff), ("All files", "*.*")],
            )
            for p in paths:
                self.add_image(p)

        def add_image(self, path: str):
            norm = os.path.normcase(os.path.abspath(path))
            for i, d in enumerate(self.docs):
                if os.path.normcase(os.path.abspath(d.path)) == norm:
                    self.set_active(i)   # already open -> just switch to it
                    return
            try:
                img = load_image(path)
            except Exception as exc:
                messagebox.showerror("Cannot open image", f"{path}\n\n{exc}", parent=self)
                return
            doc = ImageDoc(path=path, pil_image=img, width=img.size[0], height=img.size[1])
            # Auto-read the scale from the file's metadata, if present and trusted.
            auto_msg = None
            try:
                found = read_pixel_size_nm(path)
            except Exception:
                found = None
            if found and found[2]:
                nm, source, _ = found
                doc.calibration = Calibration(unit_per_px=nm, unit="nm")
                doc.scale_source = source
                auto_msg = f"Auto-detected scale for {doc.name}: {nm:g} nm/px ({source})"
            self.docs.append(doc)
            self.cancel_current()
            self._refresh_selector()          # populate combobox before selecting
            self.set_active(len(self.docs) - 1)
            self._refresh_tree()
            if auto_msg:
                self._flash_status(auto_msg)

        def set_active(self, index: int):
            if not (0 <= index < len(self.docs)):
                return
            self.active = index
            self.cancel_current()
            self._clear_selection()
            self._detected_bar = None   # highlight is specific to the previous image
            self.mode = self.MODE_IDLE
            if 0 <= index < len(self.image_selector["values"]):
                self.image_selector.current(index)
            if not self.doc.fitted:
                self.fit_to_window()
            else:
                self.redraw()
            self._refresh_summary()
            self._update_status()

        def _on_select_image(self, event):
            idx = self.image_selector.current()
            if idx >= 0:
                self.set_active(idx)

        def _image_labels(self) -> List[str]:
            """Name shown for each open image; folders are added to duplicates."""
            paths = [d.path for d in self.docs]
            if paths != self._label_paths:
                self._label_paths = paths
                self._labels = display_names(paths)
            return self._labels

        def _refresh_selector(self):
            labels = self._image_labels()
            vals = [f"{i + 1}: {labels[i]}  ({len(d.cells)} cells)" for i, d in enumerate(self.docs)]
            self.image_selector["values"] = vals
            if self.active is not None:
                self.image_selector.current(self.active)
            elif vals:
                self.image_selector.set(vals[0])
            else:
                self.image_selector.set("")

        def close_image(self):
            if self.doc is None:
                return
            if self.doc.cells and not messagebox.askyesno(
                "Close image",
                f"'{self.doc.name}' has {len(self.doc.cells)} outlined cell(s).\n"
                "Close it and discard those outlines?",
                parent=self):
                return
            del self.docs[self.active]
            if not self.docs:
                self.active = None
                self.canvas.delete("all")
                self._tk_photo = self._photo_src = self._photo_key = None   # release the image
            else:
                self.active = min(self.active, len(self.docs) - 1)
            self.cancel_current()
            self._selected = None
            self._detected_bar = None
            self._refresh_selector()
            self._refresh_tree()
            if self.doc is not None:
                if not self.doc.fitted:
                    self.fit_to_window()
                else:
                    self.redraw()
            self._update_status()

        # ------------------------------------------------------------------
        # View transforms (operate on the active doc)
        # ------------------------------------------------------------------
        def _on_canvas_configure(self, event):
            d = self.doc
            if d is not None and not d.fitted and event.width > 1:
                self.fit_to_window()
            else:
                self.redraw()

        def fit_to_window(self):
            d = self.doc
            if d is None:
                return
            cw = max(self.canvas.winfo_width(), 1)
            ch = max(self.canvas.winfo_height(), 1)
            d.scale, d.offset = fit_transform(d.width, d.height, cw, ch)
            if cw > 1 and ch > 1:
                d.fitted = True
            self.redraw()

        def zoom_center(self, factor: float):
            d = self.doc
            if d is None:
                return
            cx = self.canvas.winfo_width() / 2.0
            cy = self.canvas.winfo_height() / 2.0
            d.scale, d.offset = zoom_about(cx, cy, d.scale, d.offset, factor)
            self.redraw()
            self._update_status()

        def on_wheel(self, event, delta=None):
            d = self.doc
            if d is None:
                return
            step = delta if delta is not None else event.delta
            factor = 1.25 if step > 0 else 0.8
            d.scale, d.offset = zoom_about(event.x, event.y, d.scale, d.offset, factor)
            self.redraw()
            self._update_status()

        # ------------------------------------------------------------------
        # Panning
        # ------------------------------------------------------------------
        def toggle_pan_tool(self):
            self.pan_tool = not self.pan_tool
            self.pan_btn.config(text=f"Pan tool: {'on' if self.pan_tool else 'off'}")
            self.canvas.config(cursor="fleur" if self.pan_tool else "crosshair")

        def on_pan_start(self, event):
            if self.doc is None:
                return
            self._panning = True
            self._pan_start = (event.x, event.y)
            self._pan_offset_start = self.doc.offset

        def on_pan_move(self, event):
            if not self._panning or self.doc is None:
                return
            dx = event.x - self._pan_start[0]
            dy = event.y - self._pan_start[1]
            self.doc.offset = (self._pan_offset_start[0] + dx, self._pan_offset_start[1] + dy)
            self.redraw()

        def on_pan_end(self, event):
            self._panning = False

        # ------------------------------------------------------------------
        # Drawing style
        # ------------------------------------------------------------------
        def toggle_draw_style(self):
            self.cancel_current()
            self.draw_style = (self.STYLE_POLYGON if self.draw_style == self.STYLE_FREEHAND
                               else self.STYLE_FREEHAND)
            self.style_btn.config(
                text="Draw: Freehand" if self.draw_style == self.STYLE_FREEHAND else "Draw: Points")
            self._update_status()

        # ------------------------------------------------------------------
        # Mouse interaction
        # ------------------------------------------------------------------
        @staticmethod
        def _clamp_to_image(d: ImageDoc, ix: float, iy: float) -> Point:
            """Keep an outline vertex on the image: area beyond its edge is not real."""
            return (min(max(ix, 0.0), float(d.width)), min(max(iy, 0.0), float(d.height)))

        def on_left_down(self, event):
            d = self.doc
            if d is None:
                return
            if self._aqua:
                # Take the keyboard focus back from the table, so that the Mac
                # "delete" key undoes a point instead of deleting a cell.
                self.canvas.focus_set()
            if self.pan_tool:
                self.on_pan_start(event)
                return
            ix, iy = canvas_to_image(event.x, event.y, d.scale, d.offset)

            if self.mode == self.MODE_CALIBRATE:
                # Lock the scale bar horizontal: the second point shares the
                # first point's y, so only the horizontal distance is measured.
                if self.calib_points:
                    iy = self.calib_points[0][1]
                self.calib_points.append((ix, iy))
                if len(self.calib_points) == 2:
                    self._finish_calibration()
                self.redraw()
                return

            if self.draw_style == self.STYLE_FREEHAND:
                self._freehanding = True
                self.current_points = [self._clamp_to_image(d, ix, iy)]
                self.redraw()
            else:  # polygon
                if len(self.current_points) >= 3:
                    fx, fy = image_to_canvas(*self.current_points[0], d.scale, d.offset)
                    if math.hypot(event.x - fx, event.y - fy) <= self.CLOSE_THRESHOLD_PX:
                        self.finish_polygon()
                        return
                self.current_points.append(self._clamp_to_image(d, ix, iy))
                self.redraw()
                self._update_status()

        def on_left_drag(self, event):
            d = self.doc
            if d is None:
                return
            if self.pan_tool and self._panning:
                self.on_pan_move(event)
                return
            if self.draw_style == self.STYLE_FREEHAND and self._freehanding:
                if self.current_points:
                    lx, ly = image_to_canvas(*self.current_points[-1], d.scale, d.offset)
                    if math.hypot(event.x - lx, event.y - ly) < self.FREEHAND_SAMPLE_PX:
                        return
                ix, iy = canvas_to_image(event.x, event.y, d.scale, d.offset)
                self.current_points.append(self._clamp_to_image(d, ix, iy))
                self.redraw()

        def on_left_up(self, event):
            if self.pan_tool and self._panning:
                self.on_pan_end(event)
                return
            if self.draw_style == self.STYLE_FREEHAND and self._freehanding:
                self._finish_freehand()

        def on_double_click(self, event):
            # Tk delivers the second press of a double-click here instead of to
            # on_left_down.  In Points mode that closes the outline on the
            # points already placed.  When panning or drawing freehand it is
            # just a quick second press, which must start a drag as usual.
            if self.pan_tool or (self.mode != self.MODE_CALIBRATE
                                 and self.draw_style == self.STYLE_FREEHAND):
                self.on_left_down(event)
            elif self.draw_style == self.STYLE_POLYGON and self.current_points:
                self.finish_polygon()

        def on_right_click(self, event):
            if self.draw_style == self.STYLE_POLYGON and self.current_points:
                self.finish_polygon()

        def on_motion(self, event):
            d = self.doc
            if d is None:
                return
            self._cursor_img = canvas_to_image(event.x, event.y, d.scale, d.offset)
            if self.mode == self.MODE_CALIBRATE or (
                    self.draw_style == self.STYLE_POLYGON and self.current_points):
                self.redraw()
            self._update_status()

        # ------------------------------------------------------------------
        # Calibration (per active image)
        # ------------------------------------------------------------------
        def detect_scale_from_file(self):
            if self.doc is None:
                messagebox.showinfo("No image", "Open an image first.", parent=self)
                return
            is_tiff = os.path.splitext(self.doc.path)[1].lower() in (".tif", ".tiff")
            if is_tiff and not _module_available("tifffile"):
                messagebox.showinfo(
                    "tifffile is not available",
                    "Reading the scale from a TIFF needs the 'tifffile' package, which "
                    "is not installed or cannot be loaded. To install it:\n\n"
                    f"    {_pip_command('tifffile')}\n\n"
                    "Until then, use 'Set scale (bar)' or 'Pixel size…'.", parent=self)
                return
            try:
                found = read_pixel_size_nm(self.doc.path)
            except Exception as exc:
                found = None
                messagebox.showerror("Read failed", str(exc), parent=self)
                return
            if not found:
                messagebox.showinfo(
                    "No scale metadata",
                    "No pixel-size information was found in this file's metadata.\n\n"
                    "Use 'Set scale (bar)' to trace the scale bar, or 'Pixel size…'.",
                    parent=self)
                return
            nm, source, trusted = found
            if not trusted and not messagebox.askyesno(
                    "Low-confidence scale",
                    f"Found {nm:g} nm/px from the {source}, but it may not be a usable "
                    "microscope calibration: it looks like a print/display resolution, "
                    "or the pixels are not square.\n\n"
                    "Use it anyway?", parent=self):
                return
            self.doc.calibration = Calibration(unit_per_px=nm, unit="nm")
            self.doc.scale_source = source
            self._unsaved = True
            self._refresh_tree()
            self._update_status()
            self.redraw()
            messagebox.showinfo("Scale set", f"{nm:g} nm/px\nSource: {source}", parent=self)

        def autodetect_scale_bar(self):
            if self.doc is None:
                messagebox.showinfo("No image", "Open an image first.", parent=self)
                return
            if not _module_available("scipy.ndimage"):
                messagebox.showinfo(
                    "SciPy is not available",
                    "Detecting the scale bar needs the 'scipy' package, which is not "
                    "installed or cannot be loaded. To install it:\n\n"
                    f"    {_pip_command('scipy')}\n\n"
                    "Until then, trace the bar with 'Set scale (bar)'.", parent=self)
                return
            try:
                gray = np.asarray(self.doc.pil_image.convert("L"))
                bar = detect_scale_bar(gray)
            except Exception as exc:
                messagebox.showerror("Detection failed", str(exc), parent=self)
                return
            if bar is None:
                self._detected_bar = None
                self.redraw()
                messagebox.showinfo(
                    "No scale bar found",
                    "Could not automatically find a scale bar in this image.\n\n"
                    "Trace it manually with 'Set scale (bar)'.", parent=self)
                return
            x0, y0, x1, y1, length = bar
            self._detected_bar = bar
            # Bring the bar into view and highlight it before asking for a length.
            self._ensure_visible(x0, y0, x1, y1)
            self.redraw()
            dlg = ScaleDialog(
                self, length,
                "A scale bar was auto-detected and highlighted in magenta. If it is the "
                "correct bar, enter its real length; otherwise Cancel and trace it manually.")
            self.wait_window(dlg)
            self._detected_bar = None
            if dlg.result is not None:
                upp, unit = dlg.result
                self._apply_scale(upp, unit, "auto-detected scale bar")
            self._update_status()
            self.redraw()

        def _ensure_visible(self, x0, y0, x1, y1):
            """Pan so the given image-space box is within the canvas (no zoom change)."""
            d = self.doc
            if d is None:
                return
            cw = max(self.canvas.winfo_width(), 1)
            ch = max(self.canvas.winfo_height(), 1)
            s = d.scale
            cx0, cy0 = image_to_canvas(x0, y0, s, d.offset)
            cx1, cy1 = image_to_canvas(x1, y1, s, d.offset)
            ox, oy = d.offset
            margin = 40
            if cx0 < margin:
                ox += margin - cx0
            elif cx1 > cw - margin:
                ox -= cx1 - (cw - margin)
            if cy0 < margin:
                oy += margin - cy0
            elif cy1 > ch - margin:
                oy -= cy1 - (ch - margin)
            d.offset = (ox, oy)

        def start_calibration(self):
            if self.doc is None:
                messagebox.showinfo("No image", "Open an image first.", parent=self)
                return
            if self.pan_tool:            # otherwise the pan tool would eat the clicks
                self.toggle_pan_tool()
            self.cancel_current()
            self.mode = self.MODE_CALIBRATE
            self.calib_points = []
            self.canvas.config(cursor="tcross")
            self._update_status()
            self.redraw()

        def _finish_calibration(self):
            (x1, y1), (x2, y2) = self.calib_points
            pixel_length = math.hypot(x2 - x1, y2 - y1)
            self.mode = self.MODE_IDLE
            self.canvas.config(cursor="crosshair")
            if pixel_length < 1e-6:
                messagebox.showwarning("Too short", "The two points are identical.", parent=self)
                self.calib_points = []
                self.redraw()
                return
            dlg = ScaleDialog(self, pixel_length,
                              "Enter the real length of the scale bar you just traced "
                              "(measured horizontally).")
            self.wait_window(dlg)
            self.calib_points = []
            if dlg.result is not None:
                upp, unit = dlg.result
                self._apply_scale(upp, unit, "traced scale bar")
            self._refresh_tree()
            self._update_status()
            self.redraw()

        def enter_pixel_size(self):
            if self.doc is None:
                messagebox.showinfo("No image", "Open an image first.", parent=self)
                return
            dlg = PixelSizeDialog(self, self.doc.calibration)
            self.wait_window(dlg)
            if dlg.result is not None:
                val, unit = dlg.result
                self._apply_scale(val, unit, "manual pixel size")

        def _apply_scale(self, nm_per_px: float, unit: str, source: str):
            """Set the active image's scale; warn if it is wildly unlike the others."""
            d = self.doc
            if d is None:
                return
            d.calibration = Calibration(unit_per_px=nm_per_px, unit=unit)
            d.scale_source = source
            self._unsaved = True
            self._refresh_tree()
            self._update_status()
            self.redraw()
            others = [o.calibration.unit_per_px for o in self.docs
                      if o is not d and o.calibration.calibrated]
            if scale_is_outlier(nm_per_px, others):
                messagebox.showwarning(
                    "Check the scale",
                    f"{d.name} is now {nm_per_px:g} nm/px, which is at least 100x "
                    "different from the other open images.\n\n"
                    "If that is not intended, set the scale again and check the unit "
                    "(nm / µm / mm).", parent=self)

        # ------------------------------------------------------------------
        # Cell tracing
        # ------------------------------------------------------------------
        def _add_cell(self, points: List[Point]):
            self.doc.cells.append(Cell(points=list(points)))
            self._unsaved = True
            new_idx = len(self.doc.cells) - 1
            self._selected = (self.active, new_idx)   # restored by _refresh_tree
            self._refresh_selector()
            self._refresh_tree()
            num = self.global_number(self.active, new_idx)
            area_str, _ = self._cell_metrics(self.doc, self.doc.cells[new_idx])
            self._flash_status(f"Added cell {num}: area {area_str}")

        def _finish_freehand(self):
            self._freehanding = False
            pts = self.current_points
            self.current_points = []
            if self.doc is not None:
                pts = collapse_edge_runs(pts, self.doc.width, self.doc.height)
            if len(pts) >= 3:
                xs = [p[0] for p in pts]
                ys = [p[1] for p in pts]
                diag = math.hypot(max(xs) - min(xs), max(ys) - min(ys))
                if diag >= self.MIN_TRACE_PX:
                    simplified = dedupe_polygon(rdp_simplify(pts, self.RDP_EPSILON_PX))
                    if len(simplified) >= 3:
                        self._add_cell(simplified)
            self.redraw()
            self._update_status()

        def finish_polygon(self):
            if self.draw_style != self.STYLE_POLYGON or not self.current_points:
                return
            pts = dedupe_polygon(self.current_points)
            if len(pts) < 3:
                messagebox.showinfo("Need more points",
                                    "A cell outline needs at least 3 distinct points.", parent=self)
                return
            if self.doc is not None:
                pts = collapse_edge_runs(pts, self.doc.width, self.doc.height)
            if len(pts) < 3 or polygon_area(pts) < 1e-6:
                messagebox.showinfo(
                    "No area",
                    "These points lie on one line, so the outline encloses no area and "
                    "was not added.\n\nPoints placed beside the image are moved onto its "
                    "edge. Add a point inside the image, or press Esc to start again.",
                    parent=self)
                return
            if not polygon_is_simple(pts):
                if not messagebox.askyesno(
                        "Self-intersecting outline",
                        "This outline crosses itself, so its area may be wrong.\n\nKeep it anyway?",
                        parent=self):
                    return
            self.current_points = []
            self._add_cell(pts)
            self.redraw()
            self._update_status()

        def undo_point(self):
            if self.mode == self.MODE_CALIBRATE and self.calib_points:
                self.calib_points.pop()
                self.redraw()
            elif self.draw_style == self.STYLE_POLYGON and self.current_points:
                self.current_points.pop()
                self.redraw()
                self._update_status()

        def cancel_current(self):
            self.current_points = []
            self._freehanding = False
            self.calib_points = []
            if self.mode == self.MODE_CALIBRATE:
                self.mode = self.MODE_IDLE
                self.canvas.config(cursor="fleur" if self.pan_tool else "crosshair")
            self.redraw()
            self._update_status()

        def delete_selected_cell(self):
            # Prefer the highlighted cell; fall back to the table row selection.
            target = self._selected
            if target is None:
                sel = self.tree.selection()
                if not sel:
                    return
                idx = self.tree.index(sel[0])
                if not (0 <= idx < len(self._row_map)):
                    return
                target = self._row_map[idx]
            di, ci = target
            if 0 <= di < len(self.docs) and 0 <= ci < len(self.docs[di].cells):
                del self.docs[di].cells[ci]
                self._unsaved = True
            self._clear_selection()
            self._refresh_selector()
            self._refresh_tree()
            self.redraw()
            self._update_status()

        # ------------------------------------------------------------------
        # Rendering
        # ------------------------------------------------------------------
        def redraw(self):
            self.canvas.delete("all")
            d = self.doc
            if d is None:
                return
            cw = self.canvas.winfo_width()
            ch = self.canvas.winfo_height()
            s = d.scale
            ox, oy = d.offset

            il = max(0, int(math.floor((0 - ox) / s)))
            it = max(0, int(math.floor((0 - oy) / s)))
            ir = min(d.width, int(math.ceil((cw - ox) / s)))
            ib = min(d.height, int(math.ceil((ch - oy) / s)))
            if ir > il and ib > it:
                new_w = max(1, int(round((ir - il) * s)))
                new_h = max(1, int(round((ib - it) * s)))
                resample = Image.NEAREST if s >= 1.0 else Image.LANCZOS
                # Resampling a large image takes far longer than a mouse move,
                # and this runs for every point of a trace, so the scaled
                # picture is reused for as long as the view stays the same.
                key = (il, it, ir, ib, new_w, new_h, resample)
                if (self._tk_photo is None or self._photo_src is not d.pil_image
                        or self._photo_key != key):
                    crop = d.pil_image.crop((il, it, ir, ib)).resize((new_w, new_h), resample)
                    self._tk_photo = ImageTk.PhotoImage(crop)
                    self._photo_src = d.pil_image
                    self._photo_key = key
                self.canvas.create_image(ox + il * s, oy + it * s, anchor="nw", image=self._tk_photo)

            self._draw_overlays(d)

        def _draw_overlays(self, d: ImageDoc):
            s = d.scale
            off = d.offset

            for ci, cell in enumerate(d.cells):
                num = self.global_number(self.active, ci)
                color = PALETTE[(num - 1) % len(PALETTE)]
                selected = self._selected == (self.active, ci)
                pts = [image_to_canvas(x, y, s, off) for (x, y) in cell.points]
                if len(pts) >= 2 and self._aqua:
                    # Tk on macOS ignores stipple and would paint the cell
                    # solid, hiding the image; draw the outline only.
                    self.canvas.create_polygon(
                        [v for p in pts for v in p],
                        outline=color, fill="", width=4 if selected else 2)
                elif len(pts) >= 2:
                    self.canvas.create_polygon(
                        [v for p in pts for v in p],
                        outline=color, fill=color,
                        stipple="gray25" if selected else "gray12",
                        width=3 if selected else 2)
                if pts:
                    cx = sum(p[0] for p in pts) / len(pts)
                    cy = sum(p[1] for p in pts) / len(pts)
                    r = 10
                    self.canvas.create_oval(cx - r, cy - r, cx + r, cy + r,
                                            fill="#000000", outline=color, width=2)
                    self.canvas.create_text(cx, cy, text=str(num), fill="white",
                                            font=("", 10, "bold"))

            # In-progress freehand trace.
            if self._freehanding and self.current_points:
                pts = [image_to_canvas(x, y, s, off) for (x, y) in self.current_points]
                if len(pts) >= 2:
                    self.canvas.create_line([v for p in pts for v in p], fill="#00e5ff", width=2)
                    fx, fy = pts[0]
                    lx, ly = pts[-1]
                    self.canvas.create_line(lx, ly, fx, fy, fill="#00e5ff", width=1, dash=(4, 3))

            # In-progress polygon.
            if self.draw_style == self.STYLE_POLYGON and self.current_points and not self._freehanding:
                pts = [image_to_canvas(x, y, s, off) for (x, y) in self.current_points]
                for i in range(len(pts) - 1):
                    self.canvas.create_line(*pts[i], *pts[i + 1], fill="#00e5ff", width=2)
                if self._cursor_img is not None:
                    lx, ly = pts[-1]
                    cxp, cyp = image_to_canvas(*self._cursor_img, s, off)
                    self.canvas.create_line(lx, ly, cxp, cyp, fill="#00e5ff", width=1, dash=(4, 3))
                for i, (vx, vy) in enumerate(pts):
                    r = 4 if i == 0 else 3
                    fill = "#ffd400" if i == 0 else "#00e5ff"
                    self.canvas.create_oval(vx - r, vy - r, vx + r, vy + r, fill=fill, outline="white")

            # Auto-detected scale bar highlight (for the user to verify).
            if self._detected_bar is not None:
                x0, y0, x1, y1, length = self._detected_bar
                cx0, cy0 = image_to_canvas(x0, y0, s, off)
                cx1, cy1 = image_to_canvas(x1, y1, s, off)
                pad = 4
                self.canvas.create_rectangle(cx0 - pad, cy0 - pad, cx1 + pad, cy1 + pad,
                                             outline="#ff00ff", width=2, dash=(5, 3))
                self.canvas.create_text((cx0 + cx1) / 2, cy0 - pad - 8,
                                        text=f"detected bar: {length:g} px",
                                        fill="#ff00ff", font=("", 10, "bold"))

            # Calibration line.
            if self.mode == self.MODE_CALIBRATE:
                pts = [image_to_canvas(x, y, s, off) for (x, y) in self.calib_points]
                for (vx, vy) in pts:
                    self.canvas.create_oval(vx - 4, vy - 4, vx + 4, vy + 4, fill="#ff3b30", outline="white")
                if len(pts) == 1 and self._cursor_img is not None:
                    # Horizontal preview: hold the y of the first point.
                    cxp, _ = image_to_canvas(*self._cursor_img, s, off)
                    self.canvas.create_line(pts[0][0], pts[0][1], cxp, pts[0][1],
                                            fill="#ff3b30", width=2)
                elif len(pts) == 2:
                    self.canvas.create_line(*pts[0], *pts[1], fill="#ff3b30", width=2)

        # ------------------------------------------------------------------
        # Table + status
        # ------------------------------------------------------------------
        def _cell_metrics(self, doc: ImageDoc, cell: Cell) -> Tuple[str, str]:
            area_px = cell.area_px()
            perim_px = cell.perimeter_px()
            cal = doc.calibration
            if cal.calibrated:
                return (format_measurement(area_in_units(area_px, cal.unit_per_px), cal.area_unit()),
                        format_measurement(length_in_units(perim_px, cal.unit_per_px), cal.length_unit()))
            return (format_measurement(area_px, "px²"), format_measurement(perim_px, "px"))

        def _select_cell(self, di: int, ci: int):
            """Highlight cell (di, ci) on the canvas and select its table row."""
            self._selected = (di, ci)
            self._suppress_tree_event = True
            try:
                for item, pair in zip(self.tree.get_children(), self._row_map):
                    if pair == (di, ci):
                        self.tree.selection_set(item)
                        self.tree.see(item)
                        break
            finally:
                self._suppress_tree_event = False

        def _clear_selection(self):
            """Clear both the canvas highlight and the table row selection."""
            self._selected = None
            self._suppress_tree_event = True
            try:
                sel = self.tree.selection()
                if sel:
                    self.tree.selection_remove(*sel)
            finally:
                self._suppress_tree_event = False

        def _on_tree_select(self, event):
            if self._suppress_tree_event:
                return
            sel = self.tree.selection()
            if not sel:
                return
            idx = self.tree.index(sel[0])
            if not (0 <= idx < len(self._row_map)):
                return
            di, ci = self._row_map[idx]
            # Ignore the virtual event our own selection_set() echoes back;
            # only a genuine user click changes the selection away from
            # self._selected. Without this guard the programmatic re-select
            # would re-fire <<TreeviewSelect>> forever.
            if self._selected == (di, ci):
                return
            if di != self.active:
                self.set_active(di)   # clears selection; re-apply below
            self._select_cell(di, ci)
            self.redraw()

        def _refresh_tree(self):
            self._suppress_tree_event = True
            try:
                self.tree.delete(*self.tree.get_children())
                self._row_map = []
                labels = self._image_labels()
                for di, doc in enumerate(self.docs):
                    # The column is narrow and cut on the right, so the file
                    # name goes first and the folder that tells it apart after.
                    label = labels[di]
                    if label != doc.name and label.endswith("/" + doc.name):
                        label = f"{doc.name}  ({label[:-len(doc.name) - 1]})"
                    for ci, cell in enumerate(doc.cells):
                        num = self.global_number(di, ci)
                        area_str, perim_str = self._cell_metrics(doc, cell)
                        flag = "" if cell.is_simple() else " ⚠"
                        self.tree.insert("", "end", values=(num, label, area_str + flag, perim_str))
                        self._row_map.append((di, ci))
            finally:
                self._suppress_tree_event = False
            # Restore the row selection for the still-highlighted cell, if any.
            if self._selected is not None:
                di, ci = self._selected
                if 0 <= di < len(self.docs) and 0 <= ci < len(self.docs[di].cells):
                    self._select_cell(di, ci)
                else:
                    self._selected = None
            self._refresh_summary()

        def _refresh_summary(self):
            total = self.total_cells()
            lines = [f"Images: {len(self.docs)}    Cells: {total}"]
            d = self.doc
            if d is not None:
                cal = ("uncalibrated (px²)" if not d.calibration.calibrated
                       else f"{d.calibration.unit_per_px:g} {d.calibration.unit}/px")
                lines.append(f"Active: {d.name}")
                lines.append(f"  scale: {cal}")
                if d.calibration.calibrated and d.scale_source:
                    lines.append(f"  source: {d.scale_source}")
                lines.append(f"  cells here: {len(d.cells)}")
                if d.cells and d.calibration.calibrated:
                    upp = d.calibration.unit_per_px
                    areas = [area_in_units(c.area_px(), upp) for c in d.cells]
                    mean = sum(areas) / len(areas)
                    lines.append(f"  mean area: {format_measurement(mean, d.calibration.area_unit())}")
            self.summary_var.set("\n".join(lines))

        def _flash_status(self, msg: str):
            self.status_var.set(msg)
            self.after(2500, self._update_status)

        def _update_status(self):
            parts = []
            d = self.doc
            if d is None:
                self.status_var.set("Open one or more images to begin.")
                return
            if d.calibration.calibrated:
                parts.append(f"Scale: {d.calibration.unit_per_px:g} {d.calibration.unit}/px")
            else:
                parts.append("Scale: NOT calibrated — set it for this image (areas shown in px²)")

            if self.mode == self.MODE_CALIBRATE:
                parts.append("Calibrating — click the two ends of the scale bar "
                             "(locked horizontal; only the width is measured).")
            elif self.draw_style == self.STYLE_FREEHAND:
                parts.append("Freehand — press and drag around a cell, release to close.")
            else:
                parts.append(f"Points — click around the cell ({len(self.current_points)} placed); "
                             "click first point / Enter to close.")

            if self._cursor_img is not None:
                ix, iy = self._cursor_img
                if 0 <= ix < d.width and 0 <= iy < d.height:
                    parts.append(f"({ix:.0f}, {iy:.0f}) px  |  zoom {d.scale * 100:.0f}%")
            self.status_var.set("     ".join(parts))

        # ------------------------------------------------------------------
        # Export: CSV table + one annotated JPEG per image
        # ------------------------------------------------------------------
        def export_all(self):
            if self.total_cells() == 0:
                messagebox.showinfo("Nothing to export", "Outline at least one cell first.", parent=self)
                return
            labels = self._image_labels()
            uncal = [labels[i] for i, d in enumerate(self.docs)
                     if d.cells and not d.calibration.calibrated]
            if uncal:
                if not messagebox.askyesno(
                        "Some images are uncalibrated",
                        "These images have cells but no scale set, so their areas will be "
                        "exported in px² only:\n\n" + _listing(uncal) +
                        "\n\nExport anyway?", parent=self):
                    return
            out_dir = filedialog.askdirectory(title="Choose a folder for the export")
            if not out_dir:
                return
            csv_name = "cell_measurements.csv"
            jpg_names = self._overlay_names()
            existing = [n for n in [csv_name] + list(jpg_names.values())
                        if os.path.exists(os.path.join(out_dir, n))]
            if existing and not messagebox.askyesno(
                    "Replace existing files?",
                    f"{len(existing)} file(s) in this folder will be replaced:\n\n"
                    + _listing(existing) + "\n\nContinue?", parent=self):
                return
            written = []
            try:
                self._write_csv(os.path.join(out_dir, csv_name))
                written.append(csv_name)
                for di, fname in jpg_names.items():
                    jpg_path = os.path.join(out_dir, fname)
                    self._render_overlay_image(self.docs[di], di).save(jpg_path, "JPEG", quality=92)
                    written.append(fname)
            except Exception as exc:
                messagebox.showerror("Export failed", str(exc), parent=self)
                return
            msg = f"Saved to:\n{out_dir}\n\n" + _listing(written)
            # The CSV lists images by file name only, so say which rows belong
            # to which file when two images share a name.
            names = [d.name for d in self.docs]
            clashes = [f"cells {self.global_number(di, 0)}–"
                       f"{self.global_number(di, len(self.docs[di].cells) - 1)}: "
                       f"{labels[di]}  →  {fname}"
                       for di, fname in jpg_names.items() if names.count(names[di]) > 1]
            if clashes:
                msg += ("\n\nSome images have the same file name, and the image column of "
                        "the CSV holds the file name only. Their cells are:\n"
                        + "\n".join(clashes))
            self._show_report("Export complete", msg)

        def _overlay_names(self) -> dict:
            """File name of the annotated JPEG for each image that has cells."""
            names = {}
            used = set()
            for di, doc in enumerate(self.docs):
                if not doc.cells:
                    continue
                stem = os.path.splitext(doc.name)[0] or f"image{di + 1}"
                fname = stem + "_outlines.jpg"
                k = 1
                while fname.lower() in used:
                    fname = f"{stem}_{k}_outlines.jpg"
                    k += 1
                used.add(fname.lower())
                names[di] = fname
            return names

        def _write_csv(self, path: str):
            with open(path, "w", newline="", encoding="utf-8-sig") as f:
                w = csv.writer(f)
                w.writerow(["cell_no", "image", "calibrated", "unit_per_px", "length_unit",
                            "n_vertices", "area_px2", "area_value", "area_unit",
                            "perimeter_px", "perimeter_value", "perimeter_unit",
                            "self_intersecting"])
                for di, doc in enumerate(self.docs):
                    cal = doc.calibration
                    for ci, cell in enumerate(doc.cells):
                        num = self.global_number(di, ci)
                        area_px = cell.area_px()
                        perim_px = cell.perimeter_px()
                        if cal.calibrated:
                            area_val = f"{area_in_units(area_px, cal.unit_per_px):.6g}"
                            perim_val = f"{length_in_units(perim_px, cal.unit_per_px):.6g}"
                            area_unit = cal.area_unit()
                            perim_unit = cal.length_unit()
                            upp = f"{cal.unit_per_px:g}"
                        else:
                            area_val = perim_val = ""
                            area_unit = "px²"
                            perim_unit = "px"
                            upp = ""
                        w.writerow([num, doc.name, cal.calibrated, upp, cal.unit if cal.calibrated else "px",
                                    len(cell.points), f"{area_px:.4f}", area_val, area_unit,
                                    f"{perim_px:.4f}", perim_val, perim_unit,
                                    (not cell.is_simple())])

        def _render_overlay_image(self, doc: ImageDoc, doc_index: int) -> "Image.Image":
            overlay = doc.pil_image.convert("RGB").copy()
            draw = ImageDraw.Draw(overlay, "RGBA")
            # Small label: a thin number with a dark halo (no big disc), so it
            # sits near the centroid without covering the outline.
            fsize = max(11, doc.width // 130)
            font = load_label_font(fsize)
            lw = max(2, doc.width // 500)
            stroke = max(1, fsize // 6)
            for ci, cell in enumerate(doc.cells):
                if len(cell.points) < 2:
                    continue
                num = self.global_number(doc_index, ci)
                rgb = hex_to_rgb(PALETTE[(num - 1) % len(PALETTE)])
                draw.polygon([(x, y) for x, y in cell.points],
                             outline=rgb + (255,), fill=rgb + (55,), width=lw)
                cx = sum(p[0] for p in cell.points) / len(cell.points)
                cy = sum(p[1] for p in cell.points) / len(cell.points)
                draw.text((cx, cy), str(num), fill=(255, 255, 255, 255), font=font,
                          anchor="mm", stroke_width=stroke, stroke_fill=(0, 0, 0, 255))
            return overlay

        # ------------------------------------------------------------------
        # Session save / load (JSON) — resume work later
        # ------------------------------------------------------------------
        def save_session(self):
            if not self.docs:
                messagebox.showinfo("Nothing to save", "Open an image first.", parent=self)
                return
            if self._not_loaded:
                if not messagebox.askyesno(
                        "Save session",
                        f"{len(self._not_loaded)} image(s) of the session you loaded could not "
                        "be opened, so their outlines are not in the workspace and will not "
                        "be in the file saved now:\n\n" + _listing(self._not_loaded, 10)
                        + "\n\nSaving over the original session file would discard them. "
                        "Save anyway?", parent=self):
                    return
            path = filedialog.asksaveasfilename(
                title="Save session", defaultextension=".json",
                initialfile="cell_session.json", filetypes=[("JSON", "*.json")])
            if not path:
                return
            data = {
                "version": 1,
                "images": [
                    {
                        "path": os.path.abspath(d.path),
                        "size": [d.width, d.height],
                        "calibration": {"unit_per_px": d.calibration.unit_per_px, "unit": d.calibration.unit},
                        "scale_source": d.scale_source,
                        "cells": [{"points": [[float(x), float(y)] for x, y in c.points]} for c in d.cells],
                    }
                    for d in self.docs
                ],
            }
            try:
                with open(path, "w", encoding="utf-8") as f:
                    json.dump(data, f, indent=2)
            except Exception as exc:
                messagebox.showerror("Save failed", f"Could not write:\n{path}\n\n{exc}", parent=self)
                return
            self._unsaved = False
            self._not_loaded = []     # asked, and a file was written: do not ask on every save
            messagebox.showinfo("Saved", f"Session written to:\n{path}", parent=self)

        def load_session(self):
            if self.total_cells() > 0 and not messagebox.askyesno(
                    "Load session",
                    f"Loading a session replaces the {len(self.docs)} open image(s) and "
                    f"{self.total_cells()} outlined cell(s) currently in the workspace.\n\n"
                    "Any unsaved outlines will be lost. Continue?",
                    parent=self):
                return
            path = filedialog.askopenfilename(
                title="Load session", filetypes=[("JSON", "*.json"), ("All files", "*.*")])
            if not path:
                return
            try:
                with open(path, "r", encoding="utf-8") as f:
                    data = json.load(f)
            except Exception as exc:
                messagebox.showerror("Cannot load", str(exc), parent=self)
                return
            missing = []      # saved paths with no file to be found
            failed = []       # files that exist but could not be opened, with the reason
            relocated = []    # (saved path, path used) for images found elsewhere
            repeated = []     # entries that resolved to a file another entry uses
            resized = []
            new_docs = []
            used_files = set()
            session_dir = os.path.dirname(os.path.abspath(path))
            try:
                for entry in data.get("images", []):
                    saved_path = entry.get("path", "")
                    img_path = resolve_session_image(saved_path, session_dir)
                    if img_path is None:
                        missing.append(saved_path)
                        continue
                    try:
                        img = load_image(img_path)
                    except Exception as exc:
                        failed.append(f"{saved_path}  ({exc})")
                        continue
                    cal = entry.get("calibration") or {}
                    upp = cal.get("unit_per_px")
                    unit = cal.get("unit", "px")
                    if not isinstance(unit, str):
                        raise ValueError("a calibration unit is not text")
                    if upp is not None:
                        upp = float(upp)
                        if not math.isfinite(upp):
                            raise ValueError("a pixel size is not a finite number")
                    if upp and unit in LENGTH_UNITS_NM:   # normalise to nm/px
                        upp = to_nm(upp, unit)
                        unit = "nm"
                    doc = ImageDoc(path=img_path, pil_image=img, width=img.size[0], height=img.size[1],
                                   calibration=Calibration(unit_per_px=upp, unit=unit),
                                   scale_source=entry.get("scale_source", ""))
                    for c in entry.get("cells", []):
                        pts = [(float(x), float(y)) for x, y in c.get("points", [])]
                        if not all(math.isfinite(v) for p in pts for v in p):
                            raise ValueError("an outline has a coordinate that is not a number")
                        doc.cells.append(Cell(points=pts))
                    saved_size = entry.get("size")
                    if saved_size and list(saved_size) != [doc.width, doc.height]:
                        resized.append(len(new_docs))     # the index this image will have
                    if img_path != saved_path:
                        relocated.append((saved_path, img_path))
                    file_id = os.path.normcase(os.path.abspath(img_path))
                    if file_id in used_files:
                        repeated.append(saved_path)
                    used_files.add(file_id)
                    new_docs.append(doc)
            except (AttributeError, TypeError, ValueError, OverflowError) as exc:
                messagebox.showerror("Cannot load",
                                     f"This does not look like a valid session file.\n\n{exc}",
                                     parent=self)
                return
            if not new_docs:
                if failed:
                    detail = ("None of the images in the session could be opened:\n\n"
                              + _listing(failed))
                else:
                    detail = ("None of the images in the session could be found.\n\n"
                              "Put the images in the same folder as the session file "
                              "and load it again.")
                messagebox.showwarning("Nothing loaded", detail, parent=self)
                return
            self.docs = new_docs
            self.active = 0
            self.cancel_current()
            self._clear_selection()
            self._detected_bar = None
            self._unsaved = False
            self._not_loaded = missing + failed
            self.docs[0].fitted = False
            self._refresh_selector()
            self._refresh_tree()
            self.fit_to_window()
            self._update_status()
            msg = f"Loaded {len(new_docs)} image(s), {self.total_cells()} cell(s)."
            if relocated:
                msg += (f"\n\n{len(relocated)} image(s) were not at the path saved in the "
                        "session and were loaded from the matching path near the session "
                        "file instead. The first one:\n"
                        f"  saved: {relocated[0][0]}\n  used: {relocated[0][1]}")
            if repeated:
                msg += ("\n\nWarning: these entries point to a file that another entry "
                        "also uses, so the same image is open more than once:\n"
                        + "\n".join(repeated))
            if missing:
                msg += "\n\nCould not find:\n" + "\n".join(missing)
            if failed:
                msg += "\n\nFound but could not open:\n" + "\n".join(failed)
            if missing or failed:
                msg += ("\n\nThe outlines of those images are not loaded. Saving this "
                        "session under the same name would discard them.")
            if resized:
                labels = self._image_labels()
                msg += ("\n\nThese images are a different size than when the session was "
                        "saved, so their outlines may not line up:\n"
                        + "\n".join(labels[i] for i in resized))
            self._show_report("Loaded", msg)

    return App


def main(argv: Optional[Sequence[str]] = None):
    argv = list(sys.argv[1:] if argv is None else argv)
    try:
        App = _build_gui()
        app = App(initial_paths=argv or None)
        app.mainloop()
    except Exception:
        # When launched without a console (pythonw / double-click), a crash would
        # otherwise vanish silently. Log it and show a dialog if we can.
        import traceback
        tb = traceback.format_exc()
        try:
            log_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                    "cell_area_tool_error.log")
            with open(log_path, "w", encoding="utf-8") as f:
                f.write(tb)
        except Exception:
            log_path = "(could not write log)"
        try:
            import tkinter as tk
            from tkinter import messagebox
            root = tk.Tk()
            root.withdraw()
            messagebox.showerror(
                "TEM Cell Area Tool crashed",
                f"{tb}\n\nA copy was saved to:\n{log_path}")
            root.destroy()
        except Exception:
            sys.stderr.write(tb)
        raise


if __name__ == "__main__":
    main()
