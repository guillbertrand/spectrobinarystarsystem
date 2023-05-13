# matplotlib
from matplotlib import pyplot as plt

# numpy
import numpy as np

# astropy
from astropy.modeling import models
import astropy.units as u
from astropy import constants as const

# specutils
from specutils import Spectrum1D
from specutils.fitting import fit_lines

#


def computeConstA_Method1(rv, fwhm, snr, contrast, n):
    """
    return a constant for a given radial velocity error in km.s-1 (J.W. BRAULT)
    """
    return rv / (fwhm / (contrast * snr * np.sqrt(n)))


def computeConstB_Method2(rv, fwhm, p, snr, contrast, center):
    """
    return a constant for a given radial velocity error in km.s-1 (BOUCHY)
    """
    return rv / ((np.sqrt(fwhm) * np.sqrt(p)) / (contrast * snr))


def findCenterOfLine(spectrum1d, fwhm, model_type='gaussian'):
    """
    METHOD EXCTRACTED FROM spectroscopicbinarysystem MODULE
    Find the center of the line in a spectrum using a fit (models available : Gaussian1D, Lorentz1D, Voigt1D)
    """
    ipeak = spectrum1d.flux.argmin()
    xpeak = spectrum1d.spectral_axis[ipeak].to(u.AA)

    s = spectrum1d - 1

    match model_type:
        case 'gaussian':
            g_init = models.Gaussian1D(
                mean=xpeak, amplitude=s.flux.argmin())
        case 'voigt':
            g_init = models.Voigt1D(
                x_0=xpeak, amplitude_L=2 / (np.pi * fwhm))
        case 'lorentz':
            g_init = models.Lorentz1D(
                amplitude=s.flux.argmin(), x_0=xpeak, fwhm=fwhm)
    #

    g_fit = fit_lines(s, g_init, window=[
                      spectrum1d.spectral_axis.argmin(), spectrum1d.spectral_axis.argmax()])
    y_fit = g_fit(s.spectral_axis)

    match model_type:
        case 'gaussian':
            center = g_fit.mean
        case _:
            center = g_fit.x_0

    return (center.value, y_fit + 1*u.Jy)

#


def runMonteCarlo(n_samples, contrast, snr, center, fwhm, model, ax, first):
    """
    Generate n_samples synthetic spectra with random noise and return standard deviation
    """

    wavelength_window = (6540, 6580)
    scale = contrast / snr

    match model:
        case 'gaussian':
            line = models.Gaussian1D(-contrast, center,
                                     fwhm / 2 * np.sqrt(2 * np.log(2)))
        case 'lorentz':
            line = models.Lorentz1D(-contrast, x_0=center, fwhm=fwhm)
    #

    x = np.arange(wavelength_window[0], wavelength_window[1], .05)

    errors = []

    for i in range(n_samples):
        y = line(x) + np.random.normal(1., scale, x.shape)
        spectrum = Spectrum1D(flux=y*u.Jy, spectral_axis=x*u.AA)
        r = findCenterOfLine(spectrum, fwhm, model)
        errors.append(r[0]-center)
        if i == 0 and first:
            ax.plot(spectrum.spectral_axis, spectrum.flux, 'k-')
            ax.plot(spectrum.spectral_axis[100:-100],
                    r[1][100:-100], 'r-')
            ax.set_xlabel('Spectral Axis ({})'.format(
                spectrum.spectral_axis.unit))
            ax.set_ylabel('Flux Axis({})'.format(spectrum.flux.unit))

        std = np.std(errors)
        print(f'{i+1}/{n_samples} std={std}')
    return std


def getRV(std, lambda_ref=6562.82):
    c = const.c.to('km/s')
    return (c * std / lambda_ref)


def run(model, n_samples, snr):
    contrast = .68
    center = 6560.123
    fwhm = np.arange(2, 8, .5)
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(12, 7))
    aa, bb = [], []
    for i in range(len(fwhm)):
        # run monteCarlo
        std = runMonteCarlo(
            n_samples, contrast, snr, center, fwhm[i], model=model, ax=ax4, first=i == 1)

        # compute A constant with first method
        a = computeConstA_Method1(
            getRV(std), getRV(fwhm[i]), snr, contrast, fwhm[i]/.05)

        # compure B constant with second method
        b = computeConstB_Method2(
            getRV(std), getRV(fwhm[i]), snr, getRV(.05), contrast, center)

        print(f'rms : {std}  a : {a}  b : {b}')
        ax1.plot(a, fwhm[i], "ko")
        ax2.plot(b, fwhm[i], "ko")
        ax3.plot(std, fwhm[i], "ko")
        aa.append(a.value)
        bb.append(b.value)

    ax1.set_title(f'Method 1 : Constant A', fontname='monospace', size=8)
    ax1.set_xlabel('Constant A', fontname='monospace', size=8)
    ax1.set_ylabel('FWHM', fontname='monospace', size=8)
    ax2.set_title(f'Method 2 : Constant B', fontname='monospace', size=8)
    ax2.set_xlabel('Constant B', fontname='monospace', size=8)
    ax2.set_ylabel('FWHM', fontname='monospace', size=8)
    ax3.set_title(f'Standard deviation in km/s', fontname='monospace', size=8)
    ax3.set_xlabel('RV error (km/s)', fontname='monospace', size=8)
    ax3.set_ylabel('FWHM', fontname='monospace', size=8)
    ax1.grid(True)
    ax2.grid(True)
    ax3.grid(True)
    ax4.grid(True)
    fig.suptitle(
        f'Profil: {model}, SNR: {snr}, Contrast: {contrast}, n_samples: {n_samples}', fontsize=10)
    plt.tight_layout(pad=2, w_pad=2, h_pad=2)
    plt.savefig(
        f'rv_error_const_{model}_snr_{snr}_sample_{n_samples}.png', dpi=150)
    plt.show()

# -------------------------------------------------


plt.rcParams['font.size'] = '8'
plt.rcParams['font.family'] = 'monospace'

# 1 run with gaussian model
run('gaussian', 1000, snr=50)
run('gaussian', 1000, snr=200)
run('gaussian', 1000, snr=400)

# 2 run with lorentz model
run('lorentz', 1000, snr=50)
run('lorentz', 1000, snr=200)
run('lorentz', 1000, snr=400)
