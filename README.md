# Spectroscopic Binary Star System

!!! WORK IN PROGRESS !!!

Small script to find radial velocities of spectroscopic binaries spectra and then solves various orbital parameters with the library BinaryStarSolver and plot it with matplotlib. 

BinaryStarSolver : https://github.com/NickMilsonPhysics/BinaryStarSolver

## Requirements 
```bash
pip install -r requirements.txt # mac / linux
py -m pip install -r requirements.txt # windows
```

## Quickstart

Customize your configuration file (bss.config.yaml).

```bash
---
debug_mode: 1                                   # plot all spectra
object_name: HD123299                           # Target name (useful for Simbad query)
period: 51.41891                                # If the period of the orbit is already known use this param (period in days). 
period_guess:                                   # If the period is uncertain use this param (period in days).
title: "α Dra - Phased radial-velocities - observations collected from April 2022 to Mars 2023"
subtitle: "StarEx HR - G. Bertrand & A. Leduc"
binary_star_type: SB1                           # SB1 or SB2
spec_file_regex: '_(.+)_(\d+)_(\d+)(.*).fits'   # Fits file pattern
line_color: black                               # Define the color of the fitted velocity curve
points_color: red,black                         # Define the color of the dots ex 'red' or a color
                                                # cycle for each observer ex with 3 observers 'red,black,yellow'

radial_velocity_correction: barycentric         # The kind of velocity correction. Must be ‘barycentric’ or ‘heliocentric’.  

lambda_ref: 6562.82                             
model: voigt                                    # gaussian, voigt or lorentz
spectral_region: 10
fwhm: .5
window_width: 1.5

# result options
font_size: 9
title_font_size: 9
font_family: monospace
fig_size_x: 9
fig_size_y: 5
dpi: 150

```

And run

```bash
python bss.py sample/alphadra/ # mac / linux
py .\bss.py .\sample\alphadra\ # windows

# output
🚀 BinaryStarSystem 0.3 - Start 🚀
✨ Load configuration file 🔧  bss.config.yaml
📁 15 spectra files found !
📈 Process spectrum _alpdra_20230305_847.fits
      - Observation date : 2023-03-05T20:20:03.1803571 - 2460009.3473
      - Radial velocity correction (barycentric) : -1.1909445364983446 km / s 
      - Center of line : 6562.132045227511 Angstrom ± 0
      - Radial velocity : -32.61701629180039 km / s ± 0.0 km / s
📈 Process spectrum _hd123299_20220513_979.fits
      - Observation date : 2022-05-13T23:30:26.8146777 - 2459713.4795
      - Radial velocity correction (barycentric) : -11.41529213109784 km / s 
      - Center of line : 6563.429836848203 Angstrom ± 0
      - Radial velocity : 16.442319033294275 km / s ± 0.0 km / s
📈 Process spectrum _hd123299_20220611_884.fits
      - Observation date : 2022-06-11T21:12:37.8633742 - 2459742.3838
      - Radial velocity correction (barycentric) : -11.471399031920294 km / s
      - Center of line : 6562.151080152883 Angstrom ± 0
      - Radial velocity : -42.02794715789537 km / s ± 0.0 km / s
...
...
...
📈 Process spectrum _hd123299_20230330_806.fits
      - Observation date : 2023-03-30T19:21:12.2940714 - 2460034.3064
      - Radial velocity correction (barycentric) : -6.031714175635342 km / s 
      - Center of line : 6562.767835181292 Angstrom ± 0
      - Radial velocity : -8.414625671230416 km / s ± 0.0 km / s
📈 Process spectrum _hd123299_20230402_822.fits
      - Observation date : 2023-04-02T19:43:37.0589709 - 2460037.322
      - Radial velocity correction (barycentric) : -6.548743404662491 km / s
      - Center of line : 6562.376810950983 Angstrom ± 0
      - Radial velocity : -26.793810976748894 km / s ± 0.0 km / s
Finding initial guesses...
Minimizing...
_alpdra_20230305_847.fits phase : 0.6225235035130674
_hd123299_20220513_979.fits phase : 0.8684575382850279
_hd123299_20220527_899.fits phase : 0.13915755895947202
_hd123299_20220602_912.fits phase : 0.25609897215294897
_hd123299_20220608_885.fits phase : 0.37227413027429757
_hd123299_20220611_884.fits phase : 0.4305911968953118
_hd123299_20220613_888.fits phase : 0.4695632404468806
_hd123299_20220616_872.fits phase : 0.5276080336967461
_hd123299_20221021_784.fits phase : 0.9957992108347902
_hd123299_20230322_818.fits phase : 0.9525771744271951
_hd123299_20230326_864.fits phase : 0.03126028925719604
_hd123299_20230327_858.fits phase : 0.05058780904101878
_hd123299_20230328_814.fits phase : 0.06918408033301393
_hd123299_20230330_806.fits phase : 0.1079305259508172
_hd123299_20230402_822.fits phase : 0.1665782102374438
[γ, K, ω, e, T0, P, a, f(M)]
[-14.728, 47.275, 21.3121, 0.42084, 59874.5, 51.41891, 30.322, 0.4202]
[0.30664, 0.40488, 1.39587, 0.0069263, 0.125772, 0, 0.28103, 0.011683]
```


![results](https://github.com/guillbertrand/spectrobinarystarsystem/raw/master/sample/alphadra/bss_phased_result.png)



![debug results](https://github.com/guillbertrand/spectrobinarystarsystem/raw/master/sample/alphadra/bss_debug_result.png)

