#!/usr/bin/env python
"""Generate a synthetic TEM-like test image with a scale bar and known shapes.

Writes test_tem_image.png and test_tem_image_truth.txt (the true areas) next to
this script, for trying out cell_area_tool without real data.

    python3 make_test_image.py        (py make_test_image.py on Windows)
"""
import os

import numpy as np
from PIL import Image, ImageDraw, ImageFont

W, H = 1200, 900
NM_PER_PX = 2.5                 # ground-truth pixel size
SCALE_BAR_NM = 500.0            # label on the bar
SCALE_BAR_PX = SCALE_BAR_NM / NM_PER_PX   # = 200 px

FONT_CANDIDATES = (
    "arial.ttf",                                        # Windows
    "/System/Library/Fonts/Supplemental/Arial.ttf",    # macOS
    "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",  # Linux
    "DejaVuSans.ttf",
)

CIRCLE_R = 120
CIRCLE_C = (350, 380)
RECT = (700, 250, 1000, 430)    # 300 x 180 px
BLOB = [(500, 600), (640, 560), (760, 640), (720, 780), (560, 800), (470, 700)]


def load_font(size):
    for name in FONT_CANDIDATES:
        try:
            return ImageFont.truetype(name, size)
        except Exception:
            continue
    try:
        return ImageFont.load_default(size)
    except Exception:   # Pillow < 10.1, or built without FreeType
        return ImageFont.load_default()


# Each shape is drawn by one function, used both for the picture and for
# counting the pixels the shape really covers.  Pillow includes the boundary
# pixels, so a "radius 120" circle is 241 px wide and covers slightly more than
# pi * 120^2; the counted area is the one a careful trace should reproduce.
def draw_circle(draw, fill, outline):
    cx, cy = CIRCLE_C
    draw.ellipse([cx - CIRCLE_R, cy - CIRCLE_R, cx + CIRCLE_R, cy + CIRCLE_R],
                 fill=fill, outline=outline, width=3)


def draw_rectangle(draw, fill, outline):
    # Pillow's rectangle includes both end pixels, hence the -1.
    draw.rectangle((RECT[0], RECT[1], RECT[2] - 1, RECT[3] - 1),
                   fill=fill, outline=outline, width=3)


def draw_blob(draw, fill, outline):
    draw.polygon(BLOB, fill=fill, outline=outline)


def drawn_area(draw_shape):
    """Number of pixels ``draw_shape`` covers."""
    mask = Image.new("L", (W, H), 0)
    draw_shape(ImageDraw.Draw(mask), 255, 255)
    return int(np.count_nonzero(np.asarray(mask)))


def generate(out_dir):
    """Write the test image and its truth file into ``out_dir``; return the text."""
    rng = np.random.default_rng(0)

    # Mottled gray background resembling EM noise.
    bg = rng.normal(150, 18, size=(H, W))
    bg = np.clip(bg, 0, 255).astype(np.uint8)
    img = Image.fromarray(bg).convert("RGB")
    draw = ImageDraw.Draw(img)

    draw_circle(draw, (70, 70, 80), (30, 30, 40))
    draw_rectangle(draw, (90, 80, 70), (30, 30, 40))
    draw_blob(draw, (60, 75, 65), (25, 30, 25))

    # Scale bar.
    bar_x0, bar_y = 80, H - 70
    draw.rectangle([bar_x0, bar_y, bar_x0 + SCALE_BAR_PX - 1, bar_y + 11], fill=(255, 255, 255))
    draw.text((bar_x0, bar_y - 30), f"{SCALE_BAR_NM:.0f} nm", fill=(255, 255, 255),
              font=load_font(24))

    image_name = "test_tem_image.png"
    img.save(os.path.join(out_dir, image_name))

    lines = []
    for label, shape in ((f"circle  r={CIRCLE_R}px", draw_circle),
                         ("rectangle 300x180", draw_rectangle),
                         ("blob (6-gon)", draw_blob)):
        area_px = drawn_area(shape)
        area_nm2 = area_px * NM_PER_PX ** 2
        lines.append(f"  {label:<20} -> {area_px} px^2   = {area_nm2:.0f} nm^2 "
                     f"= {area_nm2 / 1e6:.4f} um^2")

    info = f"""Synthetic TEM test image: {image_name}
Ground-truth pixel size: {NM_PER_PX} nm/px
Scale bar: {SCALE_BAR_NM:.0f} nm spanning {SCALE_BAR_PX:.0f} px

True areas, counted from the pixels each shape covers (trace these and compare):
""" + "\n".join(lines) + """

Note: tracing a circle/blob by hand will only approximate these values.
"""
    with open(os.path.join(out_dir, "test_tem_image_truth.txt"), "w", encoding="utf-8") as f:
        f.write(info)
    return info


def main():
    print(generate(os.path.dirname(os.path.abspath(__file__))))


if __name__ == "__main__":
    main()
