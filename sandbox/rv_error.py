# matplotlib
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

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


def computeConstB_Method2(rv, fwhm, p, snr, contrast):
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

    g_init = models.Gaussian1D(mean=xpeak, amplitude=s.flux.argmin()) if model_type == 'gaussian' else models.Lorentz1D(amplitude=s.flux.argmin(), x_0=xpeak, fwhm=fwhm)
    g_fit = fit_lines(s, g_init, window=[spectrum1d.spectral_axis.argmin(), spectrum1d.spectral_axis.argmax()])
    y_fit = g_fit(s.spectral_axis)

    center = g_fit.mean if model_type == 'gaussian' else g_fit.x_0

    return (center.value, y_fit + 1*u.Jy)

#


def runMonteCarlo(n_samples, spectral_sample_value, contrast, snr, center, fwhm, model):
    """
    Generate n_samples synthetic spectra with random noise and return standard deviation
    """
    wavelength_window = (6530, 6590)
    scale = 1 / snr

    match model:
        case 'gaussian':
            line = models.Gaussian1D(-contrast, center,
                                     fwhm / 2 * np.sqrt(2 * np.log(2)))
        case 'lorentz':
            line = models.Lorentz1D(-contrast, x_0=center, fwhm=fwhm)
    #

    x = np.arange(wavelength_window[0],
                  wavelength_window[1], spectral_sample_value)

    errors = []
    synthetic_spectrum, y_fit = None, None

    for i in range(n_samples):
        y = line(x) + np.random.normal(1., scale, x.shape)
        spectrum = Spectrum1D(flux=y*u.Jy, spectral_axis=x*u.AA)
        r = findCenterOfLine(spectrum, fwhm, model)
        errors.append(r[0]-center)
        if i == 0:
            synthetic_spectrum = spectrum
            y_fit = r[1]

        std = np.std(errors)
        #print(f'{i+1}/{n_samples} std={std}')
    return std, synthetic_spectrum, y_fit


def getRV(std, lambda_ref=6562.82):
    c = const.c.to('km/s')
    return (c * std / lambda_ref)


def run(model, n_samples, snr=150, contrast=.68, spectral_sample_value=.1, fwhm=1):
    center = 6560.123

    # run monteCarlo
    std, synthetic_spectrum, y_fit = runMonteCarlo(
                                        n_samples=n_samples, 
                                        spectral_sample_value=spectral_sample_value, 
                                        contrast=contrast, 
                                        snr=snr, 
                                        center=center, 
                                        fwhm=fwhm, 
                                        model=model)

    # compute A constant with first method
    a = computeConstA_Method1(
        rv=getRV(std), 
        fwhm=getRV(fwhm, center), 
        snr=snr, 
        contrast=contrast, 
        n=fwhm/spectral_sample_value)

    # compure B constant with second method
    b = computeConstB_Method2(
        rv=getRV(std), 
        fwhm=getRV(fwhm, center), 
        p=getRV(spectral_sample_value,center), 
        snr=snr,
        contrast=contrast)

    return (a, b, std, synthetic_spectrum, y_fit)


def testWithMultipleFWHMandSNR(model='gaussian', snr=150, n_samples=100):
    plt.rcParams['font.size'] = '8'
    plt.rcParams['font.family'] = 'monospace'
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(12, 7))
    #
    a_arr, b_arr, std_arr = [],[],[]
    fwhm_arr = np.arange(1,9,1)
    snr_arr = np.linspace(50,400, 8)
    c = .68
    spectral_sample_value = .1

    cmap = plt.get_cmap('tab20')
    colors = cmap((np.arange(8)).astype(int), alpha=1)

    for ii, f in enumerate(fwhm_arr):
        for i, s in enumerate(snr_arr):
            print(f'fwhm:{f} snr:{s}')
            a, b, std, spectrum, y_fit = run(model, n_samples=n_samples, contrast=c, spectral_sample_value=spectral_sample_value, snr=s, fwhm=f)
            a_arr.append(a)
            b_arr.append(b)
            std_arr.append(std)
            ax1.plot(a, f, color=colors[i], marker="o")
            ax2.plot(b, f, color=colors[i], marker="o")
            ax3.plot(std, f, color=colors[i], marker="o") 

            if ii==5:
                ax4.plot(spectrum.spectral_axis , spectrum.flux+ (i * .1)*u.Jy, color=colors[i], label=f'SNR {s}')
                ax4.plot(spectrum.spectral_axis[100:-100] , y_fit[100:-100] + (i * .1)*u.Jy, 'r-')
                ax4.set_xlabel('Spectral Axis ({})'.format(
                    spectrum.spectral_axis.unit))
                ax4.set_ylabel('Flux Axis({})'.format(spectrum.flux.unit))

    ax4.legend()
    ax1.set_title(
        f'Method 1 : Constant A, mean={np.mean(a_arr)}', fontname='monospace', size=8)
    ax1.set_xlabel('Constant A', fontname='monospace', size=8)
    ax1.set_ylabel('FWHM', fontname='monospace', size=8)
    ax2.set_title(
        f'Method 2 : Constant B, mean={np.mean(b_arr)}', fontname='monospace', size=8)
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
        f'Profil: {model}, Contrast: {c}, n_samples: {n_samples}', fontsize=10)
    plt.tight_layout(pad=2, w_pad=2, h_pad=2)
    plt.savefig(
        f'rv_error_estimator_{model}_sample_{n_samples}.png', dpi=150)
    plt.show()

 

if __name__ == '__main__':
    testWithMultipleFWHMandSNR('gaussian', n_samples=2000)
    testWithMultipleFWHMandSNR('lorentz', n_samples=2000)

