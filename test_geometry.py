#!/usr/bin/env python
"""Headless tests for the pure geometry / calibration functions.

These cover the parts of cell_area_tool that produce numbers (areas, lengths,
unit conversion, coordinate transforms) without needing a display.  Run with:

    python3 test_geometry.py        (py test_geometry.py on Windows)

Tests that need the optional SciPy / tifffile packages are reported as SKIP
when those are not installed.
"""
import math
import os
import tempfile

import numpy as np

from cell_area_tool import (
    Cell,
    resolve_session_image,
    parse_positive_float,
    display_names,
    load_image,
    collapse_edge_runs,
    _pip_command,
    _pixsize_from_description,
    polygon_area,
    polygon_perimeter,
    polygon_is_simple,
    segments_intersect,
    area_in_units,
    length_in_units,
    fit_transform,
    canvas_to_image,
    image_to_canvas,
    zoom_about,
    normalize_for_display,
    format_measurement,
    rdp_simplify,
    point_line_distance,
    dedupe_polygon,
    to_nm,
    read_pixel_size_nm,
    _unit_to_nm,
    detect_scale_bar,
)

_checks = 0


def have(module):
    try:
        __import__(module)
        return True
    except ImportError:
        return False


class Skip(Exception):
    """Raised by a test that cannot run without an optional package."""


def need(module):
    if not have(module):
        raise Skip(f"{module} is not installed")


def approx(a, b, tol=1e-6):
    return abs(a - b) <= tol * max(1.0, abs(a), abs(b))


def check(cond, msg):
    global _checks
    _checks += 1
    if not cond:
        raise AssertionError(msg)


def test_square_area():
    sq = [(0, 0), (100, 0), (100, 100), (0, 100)]
    check(approx(polygon_area(sq), 10000.0), "square area should be 10000 px^2")
    # Order independence (clockwise gives same absolute area).
    check(approx(polygon_area(list(reversed(sq))), 10000.0), "reversed square area")


def test_triangle_area():
    tri = [(0, 0), (4, 0), (0, 3)]
    check(approx(polygon_area(tri), 6.0), "triangle area should be 6")


def test_degenerate_area():
    check(polygon_area([(0, 0), (1, 1)]) == 0.0, "2 points -> 0 area")
    check(polygon_area([]) == 0.0, "empty -> 0 area")


def test_perimeter():
    sq = [(0, 0), (100, 0), (100, 100), (0, 100)]
    check(approx(polygon_perimeter(sq, closed=True), 400.0), "closed square perimeter 400")
    check(approx(polygon_perimeter(sq, closed=False), 300.0), "open square polyline 300")
    tri = [(0, 0), (3, 0), (0, 4)]
    check(approx(polygon_perimeter(tri, closed=True), 12.0), "3-4-5 triangle perimeter 12")


def test_unit_conversion():
    # 1 px = 2 nm -> area scales by 4.
    area_px = 10000.0
    check(approx(area_in_units(area_px, 2.0), 40000.0), "area scales with pixel size squared")
    check(approx(length_in_units(400.0, 2.0), 800.0), "length scales linearly")
    # Realistic TEM: 2.5 nm/px, a 200x200 px cell -> 250000 nm^2 = 0.25 um^2.
    a = area_in_units(200 * 200, 2.5)
    check(approx(a, 250000.0), "250000 nm^2 area")


def test_self_intersection():
    square = [(0, 0), (10, 0), (10, 10), (0, 10)]
    check(polygon_is_simple(square), "square is simple")
    bowtie = [(0, 0), (10, 10), (10, 0), (0, 10)]
    check(not polygon_is_simple(bowtie), "bow-tie is not simple")
    tri = [(0, 0), (5, 0), (2, 5)]
    check(polygon_is_simple(tri), "triangle is simple")


def test_segments_intersect():
    check(segments_intersect((0, 0), (10, 10), (0, 10), (10, 0)), "crossing X")
    check(not segments_intersect((0, 0), (1, 0), (0, 1), (1, 1)), "parallel no cross")
    check(not segments_intersect((0, 0), (1, 1), (2, 2), (3, 3)), "collinear disjoint")


def test_transform_roundtrip():
    scale = 2.5
    offset = (37.0, -12.5)
    for (ix, iy) in [(0, 0), (100, 50), (-3, 999), (12.3, 45.6)]:
        cx, cy = image_to_canvas(ix, iy, scale, offset)
        rx, ry = canvas_to_image(cx, cy, scale, offset)
        check(approx(rx, ix) and approx(ry, iy), f"roundtrip failed for {(ix, iy)}")


def test_fit_transform():
    scale, (ox, oy) = fit_transform(1000, 500, 800, 800, margin=1.0)
    # Width is the binding dimension: 800/1000 = 0.8.
    check(approx(scale, 0.8), "fit scale should be 0.8")
    # Horizontally the image exactly fills the canvas -> ox == 0.
    check(approx(ox, 0.0), "fit ox centred")
    # Vertically: image height 500*0.8 = 400, centred in 800 -> oy = 200.
    check(approx(oy, 200.0), "fit oy centred")


def test_zoom_about_keeps_point_fixed():
    scale = 1.0
    offset = (0.0, 0.0)
    cursor = (400, 300)
    # Image point currently under the cursor.
    before = canvas_to_image(*cursor, scale, offset)
    new_scale, new_offset = zoom_about(*cursor, scale, offset, 2.0)
    check(approx(new_scale, 2.0), "zoom doubled scale")
    # The same image point must still project to the cursor.
    proj = image_to_canvas(*before, new_scale, new_offset)
    check(approx(proj[0], cursor[0]) and approx(proj[1], cursor[1]), "zoom keeps cursor point fixed")


def test_zoom_clamped():
    s, _ = zoom_about(0, 0, 0.02, (0, 0), 0.5, min_scale=0.02, max_scale=64)
    check(approx(s, 0.02), "zoom clamps at min")
    s, _ = zoom_about(0, 0, 64, (0, 0), 2.0, min_scale=0.02, max_scale=64)
    check(approx(s, 64.0), "zoom clamps at max")


def test_normalize_display():
    # 16-bit ramp -> spans 0..255.
    arr = np.arange(0, 65536, dtype=np.uint16).reshape(256, 256)
    out = normalize_for_display(arr)
    check(out.dtype == np.uint8, "output is uint8")
    check(out.min() == 0 and out.max() == 255, "output spans full 8-bit range")
    # uint8 passthrough.
    a8 = np.array([[0, 128, 255]], dtype=np.uint8)
    check(np.array_equal(normalize_for_display(a8), a8), "uint8 passthrough")
    # All-constant image should not crash / divide by zero.
    const = np.full((10, 10), 500, dtype=np.uint16)
    out2 = normalize_for_display(const)
    check(out2.dtype == np.uint8, "constant image handled")


def test_format_measurement():
    check(format_measurement(0, "px²") == "0 px²", "zero format")
    check("e" in format_measurement(12_345_678, "nm²"), "very large uses scientific")
    check("e" in format_measurement(0.0001, "µm²"), "small uses scientific")
    check(format_measurement(12.3456, "µm²") == "12.35 µm²", "mid format 2dp")
    # Typical nm^2 cell areas stay as readable grouped integers.
    check(format_measurement(337500, "nm²") == "337,500 nm²", "nm² grouped, no scientific")
    check(format_measurement(1_350_000, "nm²") == "1,350,000 nm²", "1.35e6 nm² grouped")


def test_to_nm():
    check(approx(to_nm(500, "nm"), 500.0), "nm passthrough")
    check(approx(to_nm(0.5, "µm"), 500.0), "µm -> nm")
    check(approx(to_nm(2, "mm"), 2_000_000.0), "mm -> nm")
    # A 0.5 µm scale bar over 200 px == a 500 nm bar over 200 px == 2.5 nm/px.
    check(approx(to_nm(0.5, "µm") / 200.0, 2.5), "scale-bar unit normalisation")


def test_point_line_distance():
    check(approx(point_line_distance((0, 5), (0, 0), (10, 0)), 5.0), "distance to x-axis")
    check(approx(point_line_distance((5, 0), (0, 0), (10, 0)), 0.0), "on the line")
    # Degenerate line (a == b) falls back to point distance.
    check(approx(point_line_distance((3, 4), (0, 0), (0, 0)), 5.0), "degenerate line")


def test_rdp_collinear_collapses():
    line = [(0, 0), (1, 0), (2, 0), (3, 0), (4, 0)]
    out = rdp_simplify(line, 0.5)
    check(out == [(0, 0), (4, 0)], "collinear points collapse to endpoints")


def test_rdp_keeps_corners():
    # A dense trace of an L-shape: the corner must survive.
    pts = [(0, 0), (2, 0), (4, 0), (6, 0), (6, 2), (6, 4), (6, 6)]
    out = rdp_simplify(pts, 0.5)
    check((6, 0) in out, "corner preserved")
    check(len(out) < len(pts), "simplification reduced vertex count")
    check(out[0] == (0, 0) and out[-1] == (6, 6), "endpoints preserved")


def test_rdp_preserves_area_of_dense_polygon():
    # Sample a circle densely, simplify, and confirm the area barely changes.
    r = 100.0
    dense = [(r * math.cos(t), r * math.sin(t))
             for t in [i * 2 * math.pi / 400 for i in range(400)]]
    a_dense = polygon_area(dense)
    simp = rdp_simplify(dense, 0.75)
    a_simp = polygon_area(simp)
    check(len(simp) < len(dense), "circle simplified")
    check(abs(a_simp - a_dense) / a_dense < 0.01, "area preserved within 1%")


def test_rdp_short_input():
    check(rdp_simplify([(0, 0), (1, 1)], 1.0) == [(0, 0), (1, 1)], "2 pts unchanged")


def test_dedupe_polygon_removes_coincident():
    # A click + double-click produces [p1, p2, p2]; must collapse to <3 pts.
    pts = [(10, 10), (50, 50), (50, 50)]
    out = dedupe_polygon(pts)
    check(len(out) == 2, "coincident duplicate removed -> 2 pts (rejected as polygon)")
    # Exactly duplicated last vertex closing the loop is dropped.
    sq = [(0, 0), (10, 0), (10, 10), (0, 10), (0, 0)]
    check(dedupe_polygon(sq) == [(0, 0), (10, 0), (10, 10), (0, 10)], "wrap duplicate dropped")


def test_dedupe_polygon_keeps_distinct():
    tri = [(0, 0), (100, 0), (50, 80)]
    out = dedupe_polygon(tri)
    check(len(out) == 3, "distinct vertices preserved")
    check(approx(polygon_area(out), polygon_area(tri)), "area unchanged after dedupe")


def test_unit_to_nm():
    check(approx(_unit_to_nm("nm"), 1.0), "nm")
    check(approx(_unit_to_nm("µm"), 1000.0), "µm")
    check(approx(_unit_to_nm("micron"), 1000.0), "micron word")
    check(approx(_unit_to_nm("cm"), 1e7), "cm")
    check(_unit_to_nm("furlong") is None, "unknown unit -> None")


def _write_tiffs(tmp):
    """Create TIFFs with various calibration metadata; return their paths."""
    import tifffile
    arr = np.zeros((64, 64), dtype=np.uint8)
    paths = {}

    # ImageJ: unit=nm, 0.4 px/nm  -> pixel size 2.5 nm.
    p = os.path.join(tmp, "ij.tif")
    tifffile.imwrite(p, arr, imagej=True, resolution=(0.4, 0.4), metadata={"unit": "nm"})
    paths["imagej"] = p

    # OME-XML: PhysicalSizeX = 2.5 nm.
    p = os.path.join(tmp, "ome.ome.tif")
    tifffile.imwrite(p, arr, metadata={"axes": "YX", "PhysicalSizeX": 2.5,
                     "PhysicalSizeXUnit": "nm", "PhysicalSizeY": 2.5, "PhysicalSizeYUnit": "nm"})
    paths["ome"] = p

    # Description text: "Pixel size = 2.5 nm".
    p = os.path.join(tmp, "desc.tif")
    tifffile.imwrite(p, arr, description="Pixel size = 2.5 nm", metadata=None)
    paths["desc"] = p

    # Plain resolution tag in cm: 4e6 px/cm -> 2.5 nm/px (trusted).
    p = os.path.join(tmp, "rescm.tif")
    tifffile.imwrite(p, arr, resolution=(4e6, 4e6), resolutionunit="CENTIMETER")
    paths["cm"] = p

    # Print resolution written in cm (28.35 px/cm == 72 dpi) -> untrusted.
    p = os.path.join(tmp, "resprintcm.tif")
    tifffile.imwrite(p, arr, resolution=(28.35, 28.35), resolutionunit="CENTIMETER")
    paths["printcm"] = p

    # Plain resolution tag in inch (72 dpi) -> untrusted print resolution.
    p = os.path.join(tmp, "resin.tif")
    tifffile.imwrite(p, arr, resolution=(72, 72), resolutionunit="INCH")
    paths["inch"] = p

    # ImageJ file that only carries a print resolution (300 dpi) -> untrusted.
    p = os.path.join(tmp, "ij_inch.tif")
    tifffile.imwrite(p, arr, imagej=True, resolution=(300, 300), metadata={"unit": "inch"})
    paths["imagej_inch"] = p

    # ImageJ, non-square pixels (2 nm x 4 nm) -> untrusted.
    p = os.path.join(tmp, "ij_aniso.tif")
    tifffile.imwrite(p, arr, imagej=True, resolution=(0.5, 0.25), metadata={"unit": "nm"})
    paths["imagej_aniso"] = p

    # OME-XML, non-square pixels -> untrusted.
    p = os.path.join(tmp, "ome_aniso.ome.tif")
    tifffile.imwrite(p, arr, metadata={"axes": "YX", "PhysicalSizeX": 2.0,
                     "PhysicalSizeXUnit": "nm", "PhysicalSizeY": 4.0, "PhysicalSizeYUnit": "nm"})
    paths["ome_aniso"] = p

    # Resolution tags in cm, non-square pixels -> untrusted.
    p = os.path.join(tmp, "rescm_aniso.tif")
    tifffile.imwrite(p, arr, resolution=(5e6, 2.5e6), resolutionunit="CENTIMETER")
    paths["cm_aniso"] = p

    # No calibration at all.
    p = os.path.join(tmp, "plain.tif")
    tifffile.imwrite(p, arr)
    paths["plain"] = p
    return paths


def test_read_pixel_size_from_metadata():
    need("tifffile")
    with tempfile.TemporaryDirectory(prefix="cat_meta_") as tmp:
        paths = _write_tiffs(tmp)

        ij = read_pixel_size_nm(paths["imagej"])
        check(ij is not None and approx(ij[0], 2.5, 1e-3) and ij[2] is True,
              f"ImageJ -> 2.5 nm/px trusted, got {ij}")

        ome = read_pixel_size_nm(paths["ome"])
        check(ome is not None and approx(ome[0], 2.5, 1e-3) and ome[2] is True,
              f"OME -> 2.5 nm/px trusted, got {ome}")

        desc = read_pixel_size_nm(paths["desc"])
        check(desc is not None and approx(desc[0], 2.5, 1e-3) and desc[2] is True,
              f"description -> 2.5 nm/px trusted, got {desc}")

        cm = read_pixel_size_nm(paths["cm"])
        check(cm is not None and approx(cm[0], 2.5, 1e-3) and cm[2] is True,
              f"cm resolution -> 2.5 nm/px trusted, got {cm}")

        printcm = read_pixel_size_nm(paths["printcm"])
        check(printcm is not None and printcm[2] is False,
              f"print resolution in cm -> untrusted, got {printcm}")

        inch = read_pixel_size_nm(paths["inch"])
        check(inch is not None and inch[2] is False, f"inch (dpi) -> untrusted, got {inch}")

        ij_inch = read_pixel_size_nm(paths["imagej_inch"])
        check(ij_inch is not None and ij_inch[2] is False,
              f"ImageJ print resolution (inch) -> untrusted, got {ij_inch}")

        for key in ("imagej_aniso", "ome_aniso", "cm_aniso"):
            res = read_pixel_size_nm(paths[key])
            check(res is not None and approx(res[0], 2.0, 1e-3) and res[2] is False,
                  f"{key}: non-square pixels -> untrusted, got {res}")

        plain = read_pixel_size_nm(paths["plain"])
        check(plain is None, f"no metadata -> None, got {plain}")


def test_pixel_size_from_description_text():
    def nm(text):
        res = _pixsize_from_description(text)
        return res[0] if res else None

    check(approx(nm("Pixel size = 2.5 nm"), 2.5), "pixel size = N nm")
    check(approx(nm("PixelWidth: 1.2e-3 um"), 1.2), "scientific notation, um")
    check(approx(nm("resolution 0.5 µm/pixel"), 500.0), "N µm/pixel")
    check(approx(nm("3 nm per px"), 3.0), "N nm per px")
    check(approx(nm("Pixel size = 1.06 Angstroms"), 0.106), "plural angstroms")
    check(approx(nm("Pixel size = 1.06 Angstrom"), 0.106), "singular angstrom")
    # A decimal comma must give no match rather than the digits after the comma.
    check(nm("Pixelsize: 0,25 nm/pixel") is None, "decimal comma not read as 25 nm/pixel")
    check(nm("1,000 nm/px") is None, "thousands comma not read as 0 nm/px")
    check(approx(nm("scale: 2.5 nm/px"), 2.5), "N nm/px after other text")
    check(approx(nm("Mag 5000x,0.25 nm/px"), 0.25), "N nm/px after a list comma")
    # The pitch of the camera sensor is not the size of an image pixel.
    check(nm("Detector pixel size = 14 um; Mag = 20000x") is None, "detector pitch ignored")
    check(nm("CMOS-pixel size=6.5um") is None, "sensor pitch ignored")
    check(approx(nm("CCD pixel size: 14 um, calibrated pixel size: 0.7 nm"), 0.7),
          "calibrated value found after the sensor pitch")
    # ... but a nanometre value is an image pixel even right after a camera name.
    for text, value in (("Detector: CCD\nPixel size: 0.5 nm", 0.5),
                        ("Microscope=JEM-1400\nCamera=CCD\nPixelSize=0.42 nm\nMag=20000", 0.42),
                        ("kV=80;Mag=20000;Camera=AMT CMOS;PixelSize=0.52nm", 0.52),
                        ("TEM camera: 0.3 nm/px", 0.3)):
        check(approx(nm(text), value), f"image pixel size after a camera name: {text!r}")
    # A scale-bar length must never be mistaken for a pixel size.
    check(nm("Scale = 500 nm") is None, "bare 'scale = N nm' is not a pixel size")
    check(nm("scale bar: 2 um") is None, "scale bar text is not a pixel size")
    check(nm("acquired on a TEM") is None, "no numbers -> None")


def test_resolve_session_image():
    with tempfile.TemporaryDirectory(prefix="cat_sess_") as tmp:
        here = os.path.join(tmp, "cell1.tif")
        open(here, "wb").close()
        check(resolve_session_image(here, "/nowhere") == here, "existing path used as-is")
        # Paths saved on another machine fall back to the session's folder,
        # whichever separator they were written with.
        check(resolve_session_image(r"C:\Users\someone\data\cell1.tif", tmp) == here,
              "Windows path resolved by file name")
        check(resolve_session_image("/Users/someone/data/cell1.tif", tmp) == here,
              "POSIX path resolved by file name")
        check(resolve_session_image("/Users/someone/data/other.tif", tmp) is None,
              "missing image -> None")
        check(resolve_session_image("", tmp) is None, "empty path -> None")
        check(resolve_session_image("/Users/someone/../cell1.tif", tmp) is None,
              "relative components rejected")

        # A project tree that moved as a whole: the session sits in results/
        # and its images in sibling folders that reuse the same file names.
        for sub in ("groupA", "groupB", "results"):
            os.makedirs(os.path.join(tmp, "proj", "raw", sub) if sub != "results"
                        else os.path.join(tmp, "proj", "results"))
        img_a = os.path.join(tmp, "proj", "raw", "groupA", "img000.tif")
        img_b = os.path.join(tmp, "proj", "raw", "groupB", "img000.tif")
        for p in (img_a, img_b):
            open(p, "wb").close()
        sess_dir = os.path.join(tmp, "proj", "results")
        check(resolve_session_image("/Users/someone/proj/raw/groupA/img000.tif", sess_dir) == img_a,
              "sibling folder resolved (groupA)")
        check(resolve_session_image(r"D:\work\proj\raw\groupB\img000.tif", sess_dir) == img_b,
              "sibling folder resolved (groupB, Windows path)")
        # A bare file name is never taken from a parent folder.
        open(os.path.join(tmp, "proj", "stray.tif"), "wb").close()
        check(resolve_session_image("/elsewhere/stray.tif", sess_dir) is None,
              "bare name not matched in a parent folder")


def test_scale_is_outlier():
    from cell_area_tool import scale_is_outlier
    others = [18.0, 18.0, 16.0]
    check(scale_is_outlier(0.018, others), "1000x too small flagged")
    check(scale_is_outlier(18000.0, others), "1000x too large flagged")
    check(not scale_is_outlier(18.0, others), "same scale not flagged")
    check(not scale_is_outlier(0.75, others), "24x magnification change not flagged")
    check(not scale_is_outlier(2.5, []), "nothing to compare with -> not flagged")


def test_parse_positive_float():
    check(parse_positive_float("2.5") == 2.5, "plain number")
    check(parse_positive_float(" 0,5 ") == 0.5, "decimal comma accepted")
    check(parse_positive_float("1e3") == 1000.0, "scientific notation")
    for bad in ("", "abc", "0", "-1", "inf", "nan"):
        check(parse_positive_float(bad) is None, f"{bad!r} rejected")
    # Thousands grouping must not be read as a decimal comma (2,000 -> 2.0).
    for grouped in ("2,000", "1,000", "10,000", "100,000", "2,000,000", "1,000.5"):
        check(parse_positive_float(grouped) is None, f"{grouped!r} rejected as ambiguous")
    for text, value in (("2,5", 2.5), ("12,25", 12.25), ("0,125", 0.125), ("0,500", 0.5),
                        ("1,2345", 1.2345), ("2000", 2000.0)):
        check(parse_positive_float(text) == value, f"{text!r} -> {value}")


def test_display_names():
    check(display_names(["/a/x.tif", "/a/y.tif"]) == ["x.tif", "y.tif"], "unique names stay short")
    check(display_names(["/d/control/001.tif", "/d/treated/001.tif", "/d/treated/002.tif"])
          == ["control/001.tif", "treated/001.tif", "002.tif"], "duplicates get their folder")
    check(display_names([r"C:\d\run1\a\001.tif", r"C:\d\run2\a\001.tif"])
          == ["run1/a/001.tif", "run2/a/001.tif"], "as many folders as needed, Windows paths")
    # One file listed twice cannot be told apart: it keeps the short name.
    check(display_names(["/d/001.tif", "/d/001.tif"]) == ["001.tif", "001.tif"],
          "identical paths keep the file name")
    check(display_names(["/m/001.tif", "/m/001.tif", "/o/001.tif"])
          == ["m/001.tif", "m/001.tif", "o/001.tif"], "identical paths next to a different one")
    check(display_names([]) == [], "no images")


def test_pip_command():
    import sys
    old = sys.executable
    try:
        sys.executable = os.path.join("opt", "py", "python3")
        check(_pip_command("scipy") == os.path.join("opt", "py", "python3") + " -m pip install scipy",
              "command uses the running interpreter")
        sys.executable = os.path.join("C:", "Program Files", "Py", "pythonw.exe")
        check(_pip_command("tifffile") == '"' + os.path.join("C:", "Program Files", "Py", "python.exe")
              + '" -m pip install tifffile', "pythonw -> python, quoted when the path has a space")
        sys.executable = ""
        check(_pip_command("numpy") == "python3 -m pip install numpy", "no interpreter path known")
    finally:
        sys.executable = old


def test_collapse_edge_runs():
    w, h = 100.0, 80.0
    inside = [(10, 10), (60, 10), (60, 50), (10, 50)]
    check(collapse_edge_runs(inside, w, h) == inside, "outline inside the image is untouched")
    box = [(0, 0), (w, 0), (w, h), (0, h)]
    check(collapse_edge_runs(box, w, h) == box, "whole-image outline keeps its corners")

    # A trace that left the image at (0, 30), overshot to (0, 70) and came
    # back in at (0, 50): the doubled-back run on the edge is removed.
    spike = [(20, 20), (0, 30), (0, 55), (0, 70), (0, 60), (0, 50), (20, 60)]
    check(not polygon_is_simple(spike), "the clamped trace touches itself")
    out = collapse_edge_runs(spike, w, h)
    check(out == [(20, 20), (0, 30), (0, 50), (20, 60)], f"run reduced to its two ends, got {out}")
    check(polygon_is_simple(out), "no longer flagged as self-intersecting")
    check(approx(polygon_area(out), polygon_area(spike)), "area unchanged")
    check(polygon_perimeter(out) < polygon_perimeter(spike), "doubled-back length removed")
    check(collapse_edge_runs(out, w, h) == out, "second pass changes nothing")

    # A run that spans the end of the list (trace started outside the image).
    wrap = [(0, 40), (0, 20), (30, 20), (30, 60), (0, 60), (0, 70)]
    out = collapse_edge_runs(wrap, w, h)
    check(approx(polygon_area(out), polygon_area(wrap)) and polygon_is_simple(out) and len(out) == 4,
          f"run across the start of the list, got {out}")

    # Round a corner: the corner itself is kept, the overshoot beyond it is not.
    corner = [(30, 30), (0, 30), (0, 0), (40, 0), (25, 0)]
    out = collapse_edge_runs(corner, w, h)
    check(out == [(30, 30), (0, 30), (0, 0), (25, 0)], f"corner kept, got {out}")

    # Everything on one edge: nothing with an area is left.
    line = [(w, 10), (w, 30), (w, 20), (w, 60)]
    check(len(collapse_edge_runs(line, w, h)) < 3, "outline beside the image has no area")


def test_cell_measurements():
    sq = Cell(points=[(0, 0), (10, 0), (10, 10), (0, 10)])
    check(approx(sq.area_px(), 100.0) and approx(sq.perimeter_px(), 40.0), "cell area/perimeter")
    check(sq.is_simple() and sq.is_simple(), "simple cell (cached)")
    check(not Cell(points=[(0, 0), (10, 10), (10, 0), (0, 10)]).is_simple(), "bow-tie cell flagged")


def test_detect_scale_bar_synthetic():
    need("scipy")
    # White bar on a mid-gray background.
    a = np.full((300, 400), 128, dtype=np.uint8)
    a[280:290, 50:250] = 255            # 200 px wide, 10 px tall
    bar = detect_scale_bar(a)
    check(bar is not None, "should find the bar")
    x0, y0, x1, y1, length = bar
    check(abs(length - 200) <= 2, f"length ~200, got {length}")
    check(275 <= y0 <= 285, f"bar near row 280, got {y0}")

    # Dark bar on a light background.
    b = np.full((300, 400), 220, dtype=np.uint8)
    b[20:28, 100:300] = 0               # 200 px wide dark bar
    barb = detect_scale_bar(b)
    check(barb is not None and abs(barb[4] - 200) <= 2, f"dark bar length, got {barb}")

    # No bar -> None (uniform noise-free mid-gray, no saturated pixels).
    c = np.full((200, 200), 128, dtype=np.uint8)
    check(detect_scale_bar(c) is None, "no bar -> None")

    # Saturated pixels touching the ends of the bar must not lengthen it.
    d = np.full((300, 400), 128, dtype=np.uint8)
    d[280:290, 50:250] = 255
    d[283, 47:50] = 255                 # 3 stray pixels on the left end
    d[287, 250:252] = 255               # 2 on the right end
    bard = detect_scale_bar(d)
    check(bard is not None and bard[4] == 200, f"stray end pixels ignored, got {bard}")
    check(bard is not None and bard[0] == 47 and bard[2] == 252, "highlight box still covers them")

    # Gaps inside the bar (a label drawn on it), a tilt or a step must not
    # shorten it.
    g = np.full((300, 400), 128, dtype=np.uint8)
    g[270:290, 50:250] = 255            # 200 x 20 px
    g[275:285, 100:160] = 0             # a 60 x 10 px hole
    barg = detect_scale_bar(g)
    check(barg is not None and barg[4] == 200, f"bar with a hole keeps its length, got {barg}")
    t = np.full((300, 600), 128, dtype=np.uint8)
    t[280:286, 100:250] = 255           # 300 x 6 px, the right half 2 px lower
    t[282:288, 250:400] = 255
    bart = detect_scale_bar(t)
    check(bart is not None and bart[4] == 300, f"tilted bar keeps its length, got {bart}")

    # A 1-px rule or streak must not outrank a real bar ...
    e = np.full((1024, 1280), 128, dtype=np.uint8)
    e[900:908, 100:300] = 255           # the bar: 200 x 8 px
    e[960, 20:1200] = 255               # a 1-px rule, 1180 px long
    e[300, 400:460] = 255               # a 1-px saturated streak
    bare = detect_scale_bar(e)
    check(bare is not None and bare[4] == 200 and bare[1] == 900,
          f"bar preferred to thin lines, got {bare}")
    # ... but a thin bar on its own is still found.
    f = np.full((300, 400), 128, dtype=np.uint8)
    f[280:282, 50:250] = 255            # 2 px high
    barf = detect_scale_bar(f)
    check(barf is not None and barf[4] == 200, f"lone 2-px bar still found, got {barf}")


def test_load_image_modes():
    from PIL import Image
    with tempfile.TemporaryDirectory(prefix="cat_img_") as tmp:
        ramp = np.tile(np.arange(64, dtype=np.uint8) * 4, (48, 1))

        gray = os.path.join(tmp, "gray.png")
        Image.fromarray(ramp).save(gray)
        img = load_image(gray)
        check(img.mode == "L" and img.size == (64, 48),
              f"8-bit grayscale stays single-channel, got {img.mode}")
        check(np.array_equal(np.asarray(img), ramp), "grayscale pixel values unchanged")
        check(np.array_equal(np.asarray(img.convert("RGB"))[..., 1], ramp),
              "same values once converted to RGB")

        rgb = os.path.join(tmp, "rgb.png")
        Image.fromarray(np.dstack([ramp, ramp // 2, ramp // 4])).save(rgb)
        check(load_image(rgb).mode == "RGB", "colour image stays RGB")

        pal = os.path.join(tmp, "pal.png")
        Image.fromarray(ramp).convert("P").save(pal)
        check(load_image(pal).mode == "RGB", "palette image -> RGB")

        rgba = os.path.join(tmp, "rgba.png")
        Image.fromarray(np.dstack([ramp, ramp, ramp, np.full_like(ramp, 255)])).save(rgba)
        check(load_image(rgba).mode == "RGB", "RGBA image -> RGB")

        deep = os.path.join(tmp, "deep.tif")
        Image.fromarray(ramp.astype(np.uint16) * 200).save(deep)
        img16 = load_image(deep)
        check(img16.mode == "L" and int(np.asarray(img16).max()) == 255,
              f"16-bit image is contrast-stretched to 8 bit, got {img16.mode}")

        # A file that is not an image reports Pillow's error, not one about TIFF.
        junk = os.path.join(tmp, "notes.png")
        with open(junk, "wb") as fh:
            fh.write(b"this is not an image")
        try:
            load_image(junk)
        except Exception as exc:  # noqa: BLE001
            check(type(exc).__name__ == "UnidentifiedImageError",
                  f"non-image reports Pillow's error, got {type(exc).__name__}: {exc}")
        else:
            check(False, "non-image file should raise")


def test_load_image_tiff_stacks():
    need("tifffile")
    import tifffile
    with tempfile.TemporaryDirectory(prefix="cat_tif_") as tmp:
        # 64-bit float frames, which Pillow may not read itself.  Every frame
        # differs, and only the first one may be shown.
        ramp = np.tile(np.linspace(0.0, 1.0, 64), (48, 1))
        frames = np.stack([ramp * (k + 1) + k for k in range(5)])
        frames[1:, 10:20, 10:20] = 0.0
        single = os.path.join(tmp, "one.tif")
        tifffile.imwrite(single, frames[0])
        first = np.asarray(load_image(single))
        check(first.shape == (48, 64) and int(first.max()) == 255, "float TIFF is stretched to 8 bit")
        for n in (2, 3, 4, 5):
            stack = os.path.join(tmp, f"stack{n}.tif")
            tifffile.imwrite(stack, frames[:n], photometric="minisblack")
            got = np.asarray(load_image(stack))
            check(got.shape == first.shape and np.array_equal(got, first),
                  f"{n}-frame stack shows its first frame")


def test_detect_scale_bar_on_test_image():
    need("scipy")
    from PIL import Image
    import make_test_image
    with tempfile.TemporaryDirectory(prefix="cat_demo_") as tmp:
        make_test_image.generate(tmp)
        gray = np.asarray(Image.open(os.path.join(tmp, "test_tem_image.png")).convert("L"))
    bar = detect_scale_bar(gray)
    check(bar is not None, "should find the 500 nm bar in the test image")
    # make_test_image draws a 200 px (= 500 nm at 2.5 nm/px) bar.
    check(bar is not None and bar[4] == 200, f"test-image bar is 200 px, got {bar[4] if bar else None}")


def test_test_image_truth_areas():
    import make_test_image as m
    # The rectangle is exact; the circle and blob are counted from the pixels
    # drawn, which is slightly more than the ideal geometric area.
    check(m.drawn_area(m.draw_rectangle) == 300 * 180, "rectangle covers 300 x 180 px")
    circle = m.drawn_area(m.draw_circle)
    check(0 < circle / (math.pi * m.CIRCLE_R ** 2) - 1 < 0.02, f"circle a little over pi r^2, got {circle}")
    blob = m.drawn_area(m.draw_blob)
    check(0 < blob / polygon_area(m.BLOB) - 1 < 0.02, f"blob a little over its shoelace area, got {blob}")


def run_all():
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_") and callable(v)]
    failed = skipped = 0
    for t in tests:
        try:
            t()
            print(f"  PASS  {t.__name__}")
        except Skip as exc:
            skipped += 1
            print(f"  SKIP  {t.__name__} ({exc})")
        except AssertionError as exc:
            failed += 1
            print(f"  FAIL  {t.__name__}: {exc}")
        except Exception as exc:  # noqa: BLE001
            failed += 1
            print(f"  ERROR {t.__name__}: {type(exc).__name__}: {exc}")
    print(f"\n{len(tests)} tests, {_checks} checks, {failed} failed, {skipped} skipped.")
    return failed


if __name__ == "__main__":
    import sys
    sys.exit(1 if run_all() else 0)
