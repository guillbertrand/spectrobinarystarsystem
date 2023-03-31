# Spectroscopic Binary Star System
Small script to find radial velocities of spectroscopic binaries spectra and then solves various orbital parameters with the library BinaryStarSolver and plot it with matplotlib. 

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

And run

```bash
python bss.py sample/alphadra/ # mac / linux
py .\bss.py .\sample\alphadra\ # windows

# output
🚀 BinaryStarSystem 0.2 - Start 🚀
✨ Load configuration file 🔧  bss.config.yaml
📁 13 spectra files found !
📈 Process spectrum _alpdra_20230305_847.fits
      - Observation date : 2023-03-05T20:20:03.1803571 - 2460009.3473
      - Radial velocity correction (barycentric) : -1.1909445364983446 km / s 
      - Centroid : 6562.084591464368 Angstrom ± 0.0 Angstrom
      - Radial velocity : -33.59378019376873 km / s ± 0.0 km / s
📈 Process spectrum _hd123299_20220513_979.fits
      - Observation date : 2022-05-13T23:30:26.8146777 - 2459713.4795
      - Radial velocity correction (barycentric) : -11.41529213109784 km / s
      - Centroid : 6562.994614212781 Angstrom ± 0.0 Angstrom
      - Radial velocity : 7.976452813187491 km / s ± 0.0 km / s
📈 Process spectrum _hd123299_20220527_899.fits
      - Observation date : 2022-05-27T21:33:55.0977123 - 2459727.3986
      - Radial velocity correction (barycentric) : -11.770828725714242 km / s 
      - Centroid : 6562.225274820505 Angstrom ± 0.0 Angstrom
      - Radial velocity : -27.16730359743505 km / s ± 0.0 km / s
📈 Process spectrum _hd123299_20220602_912.fits
      - Observation date : 2022-06-02T21:52:40.6113366 - 2459733.4116
      - Radial velocity correction (barycentric) : -11.759725353579157 km / s
      - Centroid : 6561.9931332292845 Angstrom ± 0.0 Angstrom
      - Radial velocity : -37.77163195566003 km / s ± 0.0 km / s
      ...
      ...
      ...
Finding initial guesses...
Minimizing...
[γ, K, ω, e, T0, P, a, f(M)]
[-16.185, 46.319, 21.3106, 0.4521, 60028.9, 51.41891, 29.212, 0.37572]
[0.31154, 0.42081, 1.68182, 0.0066431, 0.14047, 0, 0.28739, 0.011089]
```


![results](https://github.com/guillbertrand/spectrobinarystarsystem/raw/master/sample/alphadra/bss_phased_result.png)



![debug results](https://github.com/guillbertrand/spectrobinarystarsystem/raw/master/sample/alphadra/bss_debug_result.png)

