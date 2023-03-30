# Spectroscopic Binary Star System
Small script to find radial velocities of spectroscopic binaries fro spectra and then solves various orbital parameters with the library BinaryStarSolver and plot it with matplotlib. 

## Requirements 
```bash
pip install -r requirements.txt # mac / linux
py -m pip install -r requirements.txt # windows
```

## Quickstart

Customize your configuration file (bss.config.yaml).

```bash
---
line_width: 0.6
font_size: 9
title_font_size: 9
font_family: monospace
fig_size_x: 9
fig_size_y: 5
object_name: HD123299
period: 51.41891
period_guess: 
title: "α Dra - Phased radial-velocities - observations collected from April 2022 to Mars 2023"
subtitle: "StarEx HR - G. Bertrand & A. Leduc"
binary_star_type: SB1 # SB1 or SB2
spec_file_regex: '_(.+)_(\d+)_(\d+)(.*).fits' # Fits file pattern
line_color: black # Define the color of the fitted velocity curve
points_color: red,black # Define the color of the dots ex 'red' or a color cylcle for each observer ex with 3 observers 'red,black,yellow'
radial_velocity_correction: barycentric # The kind of velocity correction. Must be ‘barycentric’ or ‘heliocentric’.
flux_threshold: -0.5 # The threshold a pixel must be above to be considered part of a line. The threshold is positive for emission lines and negative for absorption lines.
dpi: 150
```
