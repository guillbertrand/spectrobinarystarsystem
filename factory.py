import os

import warnings
import pickle

from numba import jit

# numpy
import numpy as np

# astropy
from astropy.io import fits
import astropy.units as u
from astropy import constants as const
import astropy.wcs as fitswcs
from astropy.modeling import models, fitting

# astroquery
from astroquery.simbad import Simbad

# specutils
from specutils import Spectrum1D, SpectralRegion
from specutils.manipulation import extract_region, LinearInterpolatedResampler
from specutils.fitting import fit_lines, fit_generic_continuum, fit_continuum

import matplotlib.pyplot as plt


def syntheticBinaryStarSpectra(output_dir="synthetic_sb1", contrast_A=.68, contrast_B=.13, fwhm_A=50, fwhm_B=55, k_A=60, k_B=70, snr=300, center=0.0, region=(-200, 200), n=0.1, n_samples=150):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    np.random.seed(0)

    x = np.arange(region[0], region[1], n)
    vr_A = np.linspace(center - k_A/2,  center + k_A/2, n_samples)
    vr_A = np.concatenate((vr_A, vr_A[::-1]))
    vr_B = np.linspace(center + k_B/2,  center - k_B/2, n_samples)
    vr_B = np.concatenate((vr_B, vr_A[::-1]))
    my_wcs = fitswcs.WCS(
        header={'CDELT1': n, 'CRVAL1': x[0], 'CUNIT1': 'km/s', 'CRPIX1': 1.})

    with open(os.path.join(output_dir, 'synth_data.txt'), 'w') as f:
        for i, vr in enumerate(vr_A):
            line_A = models.Gaussian1D(-contrast_A, vr_A[i],
                                       fwhm_A / 2 * np.sqrt(2 * np.log(2)))
            y_A = line_A(x) + np.random.normal(1., 1/snr, x.shape)
            A = Spectrum1D(flux=y_A*u.Jy, wcs=my_wcs)

            line_B = models.Gaussian1D(-contrast_B, vr_B[i],
                                       fwhm_B / 2 * np.sqrt(2 * np.log(2)))
            y_B = line_B(x) + np.random.normal(1., 1/snr, x.shape)
            B = Spectrum1D(flux=y_B*u.Jy, wcs=my_wcs)

            AB = Spectrum1D(flux=A.flux + B.flux - 1 * u.Jy, wcs=my_wcs)
            AB.write(os.path.join(
                output_dir, f'synth_{i}.fits'), overwrite=True, format="wcs1d-fits")

            output = f"{-contrast_A} {-contrast_B} {fwhm_A} {fwhm_B} {vr_A[i]} {vr_B[i]}"
            f.write(output + '\n')
            # plt.plot(x, A.flux)
            # plt.plot(x, B.flux)
            # plt.plot(x, A.flux + B.flux)
            # plt.grid(True)
            # plt.show()


syntheticBinaryStarSpectra()
