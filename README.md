
# Spectroscopic Binary System

[![PyPI version](https://badge.fury.io/py/spectroscopicbinarysystem.svg)](https://badge.fury.io/py/spectroscopicbinarysystem)

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

## Usage/Examples

Download sample data (see **/examples/alphadra** directory)

And run this code :

```python
from spectroscopicbinarysystem import SpectroscopicBinarySystem

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
```

![results](https://github.com/guillbertrand/spectrobinarystarsystem/blob/master/examples/alphadra/sbs_phased_result.png)

![results](https://github.com/guillbertrand/spectrobinarystarsystem/blob/master/examples/alphadra/sbs_debug_result.png)

