
# Spectroscopic Binary System

[![PyPI version](https://badge.fury.io/py/spectroscopicbinarysystem.svg?1.1.42)](https://badge.fury.io/py/spectroscopicbinarysystem)

**Spectroscopic Binary System** is a package intended to contain functionality and some common tools needed for performing astrophysics on spectroscopic binary stars with Python. It allows, among other things, to automatically measure the radial velocity of SB1 type systems and to find their orbital solution with **BinaryStarSolver** (https://github.com/NickMilsonPhysics/BinaryStarSolver)


## Installation

To install **spectroscopicbinarysystem** with pip, run:
```bash
  # mac / unix
  pip install spectroscopicbinarysystem

  # windows
  py -m pip install spectroscopicbinarysystem
```

If you want to make sure none of your existing dependencies get upgraded, you can also do:
```bash
  # mac / unix
  pip install spectroscopicbinarysystem --no-deps 

  # windows
  py -m pip install spectroscopicbinarysystem --no-deps
```

## Prerequisites

Your spectra must be in fit(s) format with (at minimum) the following fields in the header.
The geographical coordinates of the observer will allow to automatically correct the heliocentric/baricentric velocity.
```
01 | SIMPLE = T                                   / File does conform to FITS standard
02 | BITPIX = -32                                 / Number of bits per data pixel
03 | NAXIS = 1                                    / Number of data axes
04 | NAXIS1 = 1212                                / Length of data axis 1
05 | CRVAL1 = 3780.17883300781                    / Coordinate at reference pixel
06 | CDELT1 = 3.267211914                         / Coordinate increment
09 | DATE-OBS= '2021-03-24T19:45:00'              / Date of observation start
12 | BSS_INST= 'SW72ED + StarEx 2400 + ASI183MM'  / Instrument
16 | OBSERVER= 'gbertrand'                        / Observer name or alias
17 | CUNIT1 = 'Angstrom'                          / Wavelength unit
18 | CTYPE1 = 'Wavelength'                        / Axis type
20 | CRPIX1 = 1                                   / Reference pixel
21 | BSS_VHEL= 0                                  / [km/s] Heliocentric speed
26 | JD-OBS = 0                                   / JD start observation
29 | GEO_LONG= 0                                  / Obs. geographic longitude
30 | GEO_LAT = 0                                  / Obs. geographic latitude
31 | GEO_ELEV= 0                                  / Obs. geographic elevation
```


## Usage/Examples

Download sample data (see **/examples/alphadra** directory) or from the **STAROS** database (https://alphadra.staros-projects.org/)

And run this code :

```python
from spectroscopicbinarysystem import SpectroscopicBinarySystem

# create SpectroscopicBinarySystem object
sbs = SpectroscopicBinarySystem(
    object_name='hd123299',
    spectra_path='./examples/alphadra/',
    t0=2451441.804,
    period_guess=51,
    conf={
        "LAMBDA_REF": 6562.82,
        "LINE_FIT_MODEL": "voigt",
        "LINE_FIT_WINDOW_WIDTH": 10,
        "LINE_FIT_CONT_NORM_EXCLUDE_WIDTH": 1.5,
        "LINE_FIT_FWHM": .5,
        "RV_CORR_TYPE": "barycentric",
        "SB_TYPE": 1
    },
    debug=True)

# plot result with matplotlib and save the results
sbs.plotRadialVelocityCurve(title="α Dra - HD123299 - Phased radial velocities", 
                            subtitle=f"{sbs.getObservationCount()} observations collected from april 2022 to april 2023",
                            savefig=True)
 
# display result with plotly
sbs.plotlyRadialVelocityCurve(
    title="α Dra - HD123299 - Phased radial velocities")

# plot 2d dynamic spectrum
sbs.plotSpec2DFlux(
    title="α Dra - HD123299 - 2d dynamic spectrum",
    subtitle=f"{sbs.getObservationCount()} observations collected from april 2022 to april 2023",
    savefig=False
)
```

![results](https://github.com/guillbertrand/spectrobinarystarsystem/blob/master/examples/alphadra/hd123299_phased_result.png)

![results](https://github.com/guillbertrand/spectrobinarystarsystem/blob/master/examples/alphadra/hd123299_2d_spectrum_result.png)

![results](https://github.com/guillbertrand/spectrobinarystarsystem/blob/master/examples/alphadra/sbs_debug_result.png)

