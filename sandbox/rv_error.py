# matplotlib
from matplotlib import pyplot as plt

# numpy
import numpy as np

# astropy
from astropy.modeling import models
import astropy.units as u
from astropy import constants as const

# specutils
from specutils import Spectrum1D, SpectralRegion
from specutils.manipulation import extract_region, LinearInterpolatedResampler
from specutils.fitting import fit_lines, fit_generic_continuum

#


def computeConstA_Method1(rv, fwhm, snr, contrast, n):
    return rv / (fwhm / (contrast * snr * np.sqrt(n)))


def computeConstB_Method2(rv, fwhm, snr, contrast, center):
    return rv / ((np.sqrt(fwhm) * np.sqrt((fwhm/center)*const.c)) / (contrast * snr))


def findCenterOfLine(spectrum1d, window_width, norm_exclude_width, fwhm, model_type='gaussian'):
    """
    >> METHOD EXCTRACTED FROM spectroscopicbinarysystem MODULE
    Find the center of the line in a spectrum using a fit (models available : Gaussian1D, Lorentz1D, Voigt1D)
    :return: None
    """
    w1 = float(window_width) / 2
    w2 = float(norm_exclude_width) / 2
    fwhm = float(fwhm)

    ipeak = spectrum1d.flux.argmin()
    xpeak = spectrum1d.spectral_axis[ipeak].to(u.AA)

    invert_s = extract_region(Spectrum1D(flux=spectrum1d.flux * -1, wcs=spectrum1d.wcs),
                              SpectralRegion(xpeak - (w1 * u.AA), xpeak + (w1 * u.AA)))

    s_fit = fit_generic_continuum(invert_s,
                                  exclude_regions=[SpectralRegion(xpeak - (w2 * u.AA), xpeak + (w2 * u.AA))])
    invert_s = invert_s / s_fit(invert_s.spectral_axis)
    invert_s = invert_s - 1

    match model_type:
        case 'gaussian':
            g_init = models.Gaussian1D(
                mean=xpeak, amplitude=invert_s.flux.argmax())
        case 'voigt':
            g_init = models.Voigt1D(
                x_0=xpeak, amplitude_L=2 / (np.pi * fwhm))
        case 'lorentz':
            g_init = models.Lorentz1D(x_0=xpeak, fwhm=fwhm)
        case _:
            g_init = models.Gaussian1D(
                mean=xpeak, amplitude=invert_s.flux.argmax())
    #

    g_fit = fit_lines(invert_s, g_init, window=SpectralRegion(
        xpeak - (w2 * u.AA), xpeak + (w2 * u.AA)))
    y_fit = g_fit(invert_s.spectral_axis)

    match model_type:
        case 'gaussian':
            center = g_fit.mean
        case _:
            center = g_fit.x_0

    return (center.value, invert_s, y_fit)

#


def runMonteCarlo(n_samples, contrast, snr, center, fwhm, model, ax, first):
    """
    Generate n_samples synthetic spectra with random noise and return standard deviation
    """

    wavelength_window = (6520, 6600)
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
        r = findCenterOfLine(spectrum, fwhm*8, fwhm*4, fwhm, model)
        errors.append(r[0]-center)
        # plot the first synthetic spectrum for debug purpose
        if i == 1 and first:
            ax.plot(spectrum.spectral_axis, spectrum.flux)
            ax.plot(r[1].spectral_axis, r[2]+1)
            ax.set_xlabel('Spectral Axis ({})'.format(
                spectrum.spectral_axis.unit))
            ax.set_ylabel('Flux Axis({})'.format(spectrum.flux.unit))

        std = np.std(errors)
        print(f'{i}/{n_samples} std={std}')
    return std


def getRV(std, lambda_ref=6562.82):
    c = const.c.to('km/s')
    return (c * std / lambda_ref)


def showPlot(fig, snr, contrast, n_samples, model):

    ax1.set_title(f'Method 1')
    ax1.set_xlabel('Constant A')
    ax1.set_ylabel('FWHM')
    ax2.set_title(f'Method 2')
    ax2.set_xlabel('Constant B')
    ax2.set_ylabel('FWHM')
    ax3.set_title(f'Standard deviation in km/s')
    ax3.set_xlabel('RV incertitude (km/s)')
    ax3.set_ylabel('FWHM')
    ax1.grid(True)
    ax2.grid(True)
    ax3.grid(True)
    ax4.grid(True)
    fig.suptitle(
        f'Profil: {model}, SNR: {snr}, Contrast: {contrast}, n_samples: {n_samples}', fontsize=10)
    plt.tight_layout(pad=2, w_pad=1, h_pad=1)
    plt.savefig(f'rv_incertitude_const_{model}_snr_{snr}.png', dpi=150)
    plt.show()


#
n_samples = 2
contrast = .68
snr = 200
center = 6560.123
fwhm = np.arange(4, 8, .5)

# first run with gaussian model
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(8, 7))
aa, bb = [], []
for i in range(len(fwhm)):
    std = runMonteCarlo(
        n_samples, contrast, snr, center, fwhm[i], model='gaussian', ax=ax4, first=i == 1)
    a = computeConstA_Method1(
        getRV(std), fwhm[i], snr, contrast, .05*fwhm[i])
    b = computeConstB_Method2(
        getRV(std), fwhm[i], snr, contrast, center)
    print(f'rms : {std}  a : {a}  b : {b}')
    ax1.plot(a, fwhm[i], "ko")
    ax2.plot(b, fwhm[i], "ko")
    ax3.plot(std, fwhm[i], "ko")
    aa.append(a.value)
    bb.append(b.value)

showPlot(fig, snr, contrast, n_samples, 'gaussian')

# # second run with lorentz model
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(8, 7))
aa, bb = [], []
for i in range(len(fwhm)):
    std = runMonteCarlo(
        n_samples, contrast, snr, center, fwhm[i], model='lorentz', ax=ax4, first=i == 1)
    a = computeConstA_Method1(
        getRV(std), fwhm[i], snr, contrast, .05*fwhm[i])
    b = computeConstB_Method2(
        getRV(std), fwhm[i], snr, contrast, center)
    print(f'rms : {std}  a : {a}  b : {b}')
    ax1.plot(a, fwhm[i], "ko")
    ax2.plot(b, fwhm[i], "ko")
    ax3.plot(std, fwhm[i], "ko")
    aa.append(a.value)
    bb.append(b.value)

showPlot(fig, snr, contrast, n_samples, 'lorentz')
