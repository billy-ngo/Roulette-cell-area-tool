# Differentiation project data

Areas and perimeters of feeder and seeker cells, outlined in TEM images with
the cell area tool in this repository.

| Cell type | Image set | Images | Cells |
| --- | --- | ---: | ---: |
| Feeder | `zap` | 94 | 259 |
| Feeder | `grid3` | 6 | 14 |
| Seeker | `grid7` | 1 | 225 |

In total: 273 feeders and 225 seekers.

## Files

- `all_cells.csv`: all 498 cells in one table. The first two columns give the
  cell type and the image set; the rest are the columns of the tool's export
  (see *Output* in the main README).
- `feeders/<image set>/` and `seekers/<image set>/`: the tool's export for each
  image set, i.e. `cell_measurements.csv` and one `<image>_outlines.jpg` per
  image showing the numbered outlines.

Areas are in nm² and perimeters in nm. The pixel size is 17.8571 nm/px for
the `zap` set and 13.4228 nm/px for `grid3` and `grid7`.

## Notes

- Image `zap042` is excluded because its scale bar was annotated incorrectly.
  Its three cells were numbers 107–109 of the `zap` set, so those numbers are
  absent from that table.
- Cell numbers restart in each image set and match the labels in the outline
  images.
- For the `self_intersecting` column and the accuracy of freehand outlines,
  see *Measurement notes* in the main README.
