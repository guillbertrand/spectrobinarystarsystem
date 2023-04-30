import re
import os
import copy
import math
import warnings

# numpy
import numpy as np

# astropy
from astropy.io import fits
from astropy.time import Time
import astropy.units as u
import astropy.wcs as fitswcs
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, SpectralCoord
from astropy.modeling import models
from astropy import constants as const

# astroquery
from astroquery.simbad import Simbad

# specutils
from specutils import Spectrum1D, SpectralRegion
from specutils.manipulation import extract_region, LinearInterpolatedResampler
from specutils.fitting import fit_lines, fit_generic_continuum

# matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator)

# plotly
import plotly.graph_objects as go

# orbitalpy
from orbital import utilities

# binarystarsolve
from binarystarsolve.binarystarsolve import StarSolve

#


class SBSpectrum1D(Spectrum1D):
    """
    Class which extends the Spectrum1D class to automatically calculate the radial velocity of an absorption line
    """

    def __init__(self, filename, skycoord, conf):
        self._filename = filename
        self._basename = os.path.basename(filename)
        self._skycoord = skycoord
        self._conf = conf

        spectrum_file = fits.open(filename)

        self._header = spectrum_file[0].header
        self._jd = self._header['JD-OBS']
        self._observer = self._header['OBSERVER']
        with warnings.catch_warnings():  # Ignore warnings
            warnings.simplefilter('ignore')

            wcs_data = fitswcs.WCS(header={'CDELT1': self._header['CDELT1'], 'CRVAL1': self._header['CRVAL1'],
                                           'CUNIT1': 'Angstroms', 'CTYPE1': 'WAVE',
                                           'CRPIX1': self._header['CRPIX1']})
            flux = spectrum_file[0].data * u.Jy

            # init Spectrum1D
            super().__init__(flux=flux, wcs=wcs_data)

            # analyse spectrum
            self.findCenterOfLine()
            self.findRVCorrection()
            self.findRV()

    def getObserver(self):
        """
        Return 'OBSERVER' header field
        :return: string corresponding to the observer
        """
        return self._observer.lower()

    def getInstrument(self):
        """
        Return 'BSS_INST' header field
        :return: string corresponding to the instrument
        :rtype: string
        """
        return self._header['BSS_INST']

    def getRV(self):
        """
        Return the final computed radial velocity of the line (already corrected of heliocentric/barycentric velocity)
        :return: radial velocity
        :rtype: float
        """
        return self._rv.value

    def getJD(self):
        """
        Return 'JD-OBS' header field
        :return: float corresponding to the julian date of the observation
        :rtype: float
        """
        return float(self._header['JD-OBS'])

    def getDate(self):
        """
        Return 'DATE-OBS' header field
        :return: string corresponding to the date of the observation
        :rtype: string
        """
        return self._header['DATE-OBS']

    def setPhase(self, phase):
        """
        Set phase of the observation
        :param phase: float corresponding to the phase of the system
        """
        self._phase = phase

    def getPhase(self):
        """
        Return the phase of the observation
        :return: float corresponding to the phase of the system
        :rtype: float
        """
        return self._phase

    def getBaseName(self):
        """
        Return base filename of the spectra
        :return: string corresponding to the name of the spectra
        :rtype: string
        """
        return self._basename

    def findRVCorrection(self):
        """
        Compute radial velocity correction in function of the target and the location of the observer
        :return: None
        """
        t = Time(self._header['JD-OBS'], format='jd', scale='utc')
        loc = EarthLocation(
            self._header['GEO_LONG'], self._header['GEO_LAT'], self._header['GEO_ELEV'] * u.m)
        vcorr = self._skycoord.radial_velocity_correction(
            kind=self._conf["RV_CORR_TYPE"], obstime=t, location=loc)
        self._rv_corr = vcorr.to(u.km / u.s)

    def getRVCorrection(self):
        return self._rv_corr

    def findRV(self):
        """
        Compute the radial velocity from a line position and a radial velocity correction
        Reference line rest position in Angstroms can be customize in conf["LAMBDA_REF"]
        :return: None
        """
        c = const.c.to('km/s')
        self._rv = (c * ((self._center_of_line -
                    self._conf["LAMBDA_REF"]) / self._conf["LAMBDA_REF"])) + self._rv_corr

    def getCenterOfLine(self):
        """
        Return the center of the line computed for the spectrum using a fit (non corrected from heliocentric/barycentric velocity)
        :return: Float
        """
        return self._center_of_line

    def getDebugLineFitting(self):
        return self._debug_line_fitting

    def findCenterOfLine(self):
        """
        Find the center of the line in a spectrum using a fit (models available : Gaussian1D, Lorentz1D, Voigt1D)
        :return: None
        """
        w1 = float(self._conf['LINE_FIT_WINDOW_WIDTH']) / 2
        w2 = float(self._conf['LINE_FIT_CONT_NORM_EXCLUDE_WIDTH']) / 2
        fwhm = float(self._conf['LINE_FIT_FWHM'])

        # > Hack to replace zero values of the spectra.
        # Sometimes first value for the intensity is equal to 0.00 and prevents finding the true minima.
        # Ex : sample/alphadra/_alphadra_20230406_856.fits
        self.flux[self.flux == 0] = 1 * u.Jy
        ipeak = self.flux.argmin()
        xpeak = self.spectral_axis[ipeak].to(u.AA)

        invert_s = extract_region(Spectrum1D(flux=self.flux * -1, wcs=self.wcs),
                                  SpectralRegion(xpeak - (w1 * u.AA), xpeak + (w1 * u.AA)))

        s_fit = fit_generic_continuum(invert_s,
                                      exclude_regions=[SpectralRegion(xpeak - (w2 * u.AA), xpeak + (w2 * u.AA))])
        invert_s = invert_s / s_fit(invert_s.spectral_axis)
        invert_s = invert_s - 1

        match self._conf['LINE_FIT_MODEL']:
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

        match self._conf['LINE_FIT_MODEL']:
            case 'gaussian':
                center = g_fit.mean
            case _:
                center = g_fit.x_0

        self._debug_line_fitting = (invert_s, y_fit)
        self._center_of_line = center.value

    def __str__(self):
        return f"Spectrum : {self._basename}\n- obs: {self._observer}\n- jd: {self._jd}\n- center: {self._center_of_line} A\n- {self._conf['RV_CORR_TYPE']}: {self._rv_corr}\n- rv: {self._rv}\n "

#


class SpectroscopicBinarySystem:
    """
    Class allowing to dynamically load spectra in fit(s) format and to
    obtain the orbital solution of a spectroscopic binary system.

    Also allows to:
    - Plot the solution and the measured radial velocity points.
    - Plot a dynamic 2D spectrum

    :param object_name: Common name of the target for the Simbad query (ex : alpha dra, hd123299)
    :param spectra_path: Path of all fit(s) spectra
    :param t0: Allow to fix Periastron epoch T0 if known (julian date)
    :param period: If the period of the orbit is already known use this param (period in days).
    :param perdiod_guess: If the period is uncertain use this param (period in days).
    :param conf: Allow to customize additionnal parameters (see self._conf default value below)
    :param debug: Allow to activate log and some additional plots for debug purpose
    """

    def __init__(self, object_name, spectra_path, t0=None, period=None, period_guess=None, conf=None, debug=False):

        self._conf = {"LAMBDA_REF": 6562.82,
                      "LINE_FIT_MODEL": "voigt",
                      "LINE_FIT_WINDOW_WIDTH": 10,
                      "LINE_FIT_CONT_NORM_EXCLUDE_WIDTH": 1.5,
                      "LINE_FIT_FWHM": .5,
                      "RV_CORR_TYPE": "barycentric",
                      "SB_TYPE": 1}

        self._sb_spectra = []
        self._spectra_filename = []
        self._orbital_solution = None
        self._spectra_path = spectra_path
        self._object_name = object_name
        self._type = type
        self._t0 = t0
        self._period = period
        self._period_guess = period_guess
        self._debug = debug

        # load user configuration or defaults
        if conf:
            self._conf.update(conf)

        print('** SpectroscopicBinarySystem **')
        self.__findObjectCoordinate()
        print('Load spectra...')
        self.__loadSpectra()

    def __findObjectCoordinate(self):
        # try to find the coordinates with the common name of the object
        self._skycoord = None
        if result_table := Simbad.query_object(self._object_name):
            ra = result_table[0]['RA']
            dec = result_table[0]['DEC']
            self._skycoord = SkyCoord(f'{ra} {dec}', unit=(u.hourangle, u.deg))

    def __loadSpectra(self):
        for root, dirs, files in os.walk(self._spectra_path):
            for file in files:
                regex = re.compile('(.*).fit')
                if (re.match(regex, file)):
                    spectrum_filename = os.path.join(self._spectra_path, file)
                    sbSpec1D = SBSpectrum1D(
                        spectrum_filename, self._skycoord, self._conf)
                    self._sb_spectra.append(sbSpec1D)
                    if self._debug:
                        print(sbSpec1D)
        if self._debug:
            print(f'{len(self._sb_spectra)} processed spectra')
            plt.rcParams['font.size'] = '6'
            plt.rcParams['font.family'] = 'monospace'
            grid_size = math.ceil(len(self._sb_spectra)/6)
            fig, axs = plt.subplots(6, grid_size, figsize=(
                13, 7), sharex=True, sharey=True)
            for i, s in enumerate(self._sb_spectra):
                ax = axs.flat[i]
                extracted_profil, line_fitting = s.getDebugLineFitting()
                ax.set_title(
                    f'{s.getBaseName()}\n{s.getObserver()} JD={s.getJD()}', fontsize="6")
                ax.grid(True)
                ax.tick_params(axis='both', which='major', labelsize=6)
                ax.tick_params(axis='both', which='minor', labelsize=6)
                ax.plot(extracted_profil.spectral_axis.to(
                    u.AA), extracted_profil.flux, color="k")
                ax.plot(extracted_profil.spectral_axis.to(
                    u.AA), line_fitting, color="r")
                ax.axvline(x=s.getCenterOfLine(),
                           color='r', linestyle='-', lw=0.7)
            plt.tight_layout(pad=0.8, w_pad=2, h_pad=1)
            plt.savefig(f'{self._spectra_path}/sbs_debug_result.png', dpi=150)
            plt.show()

    def getObservationCount(self):
        """
        Return the count of processed spectra
        :return: count
        :rtype: int
        """
        return len(self._sb_spectra)

    def __getPhase(self, jd0, period, jd):
        """
        Compute the phase of a given JD
        :param jd0: JD of the first observation
        :param period: period of the orbit
        :param jd: JD to compute the phase
        :return: phase
        """
        return (jd - jd0) / period % 1

    def __computeRadialVelocityCurve(self, t, t0, K, e, w, v0):
        """
        Compute the radial velocity curve
        :param t: time
        :param t0: time of the first observation
        :param K: semi-amplitude
        :param e: eccentricity
        :param w: argument of periastron
        :param v0: systemic velocity
        :return: radial velocity
        """
        w = math.radians(w)
        # Mean anomaly
        M = 2 * np.pi * ((t - t0) % 1)
        # Eccentric anomaly
        E = utilities.eccentric_anomaly_from_mean(e, M, tolerance=0.00001)
        # True anomaly
        f = utilities.true_anomaly_from_eccentric(e, E)
        return (K * (e * np.cos(w) + np.cos(w + f)) + v0)

    def __solveSystem(self):
        """
        Compute the orbital solution with BinaryStarSolver
        """
        # write result file for BinaryStarSolver
        with open(f'{self._spectra_path}/sbs_results.txt', 'w') as f:
            for s in self._sb_spectra:
                output = f"{float(s.getJD()) - 2400000.0} {round(s.getRV(), 3)}"
                f.write(output + '\n')

        # [γ, K, ω, e, T0, P, a, f(M)]
        params, err, cov = StarSolve(
            data_file=f"{self._spectra_path}/sbs_results.txt",
            star="primary",
            Period=self._period,
            Pguess=self._period_guess,
            covariance=True,
            graphs=False,
        )

        self._orbital_solution = (params, err, cov)

        # If self._t0 compute phase delta between T0 of the model and the fixed value
        if self._t0:
            t0 = self._t0
            phase1 = self.__getPhase(
                t0, self._orbital_solution[0][5], self._orbital_solution[0][4]+2400000)
            phase2 = 1.0
            self._v0 = phase1 - phase2
        # Else use t0 compute by the model
        else:
            self._t0 = self._orbital_solution[0][4] + 2400000
            self._v0 = 0

        period = self._orbital_solution[0][5]
        for s in self._sb_spectra:
            # compute phase of the sytem
            jd = s.getJD()
            phase = self.__getPhase(float(t0), period, jd)
            s.setPhase(phase)
            if self._debug:
                print(f"{s.getBaseName()} phase : {phase}")

        print(
            f'{self._object_name} orbital solution with {len(self._sb_spectra)} spectra',
            f'- γ = {params[0]} ± {err[0]}',
            f'- K = {params[1]} ± {err[1]}',
            f'- ω = {params[2]} ± {err[2]}',
            f'- e = {params[3]} ± {err[3]}',
            f'- T0 = {params[4]} ± {err[4]}',
            f'- P = {params[5]} ± {err[5]}',
            f'- a = {params[6]} ± {err[6]}',
            f'- f(M) = {params[7]} ± {err[7]}',
            sep='\n')

    def getOrbitalSolution(self):
        self.__solveSystem()
        return self._orbital_solution

    def __findNearest(self, array, value):
        """
        Find the nearest value in an array
        :param array: array for search
        :param value: value to find
        :return: index of the nearest value
        """
        array = np.asarray(array)
        return (np.abs(array - value)).argmin()

    def __plotRadialVelocityDots(self, axs, t0):
        """
        Plot the radial velocity dots from the data
        :param axs: axes to plot
        :model t0: Periastron epoch T0 (julian date)
        :return: None
        """
        observers = {}
        instruments = {}
        marker_index = {}
        color_number = 0

        cmap = plt.get_cmap('tab20')
        markers_style = ["o", "v", "^", "s", "D", "P", "X"]

        # define colors (max 60 distinct observers)
        colors = cmap((np.arange(20)).astype(int), alpha=1)
        + cmap((np.arange(20)).astype(int), alpha=.75)
        + cmap((np.arange(20)).astype(int), alpha=.5)

        # sort sb spectra by observer name
        self._sb_spectra.sort(key=lambda x: x.getObserver())

        for s in self._sb_spectra:

            # get the observer
            obs = s.getObserver()
            if (obs not in observers.keys()):
                observers[obs] = colors[color_number]
                color_number += 1
            color = observers[obs]

            # get the instrument

            label = f"{obs} - {s.getInstrument()[:30]}…"
            if label not in instruments.keys():
                if obs not in marker_index:
                    marker_index[obs] = 0
                instruments[label] = markers_style[marker_index[obs]]
                marker_index[obs] += 1
                axs[0].errorbar(s.getPhase(), s.getRV(
                ), yerr=0, label=label, ecolor='k', capsize=0, fmt=instruments[label], color=color, lw=0.7, markersize=5)
            else:
                axs[0].errorbar(s.getPhase(), s.getRV(), yerr=0,
                                fmt=instruments[label], ecolor='k', capsize=0, color=color, lw=.7, markersize=5)

            xindex = self.__findNearest(self._model_x, s.getPhase())
            axs[1].errorbar(s.getPhase(), s.getRV() - self._model_y[xindex],
                            yerr=0, fmt=instruments[label], ecolor='k', capsize=0, color=color, lw=.7, markersize=5)

    def plotRadialVelocityCurve(self, title="", subtitle="", rv_y_multiple=10, residual_y_multiple=0.5, savefig=False, dpi=150, font_family='monospace', font_size=9):
        if not self._orbital_solution:
            self.__solveSystem()

        plt.rcParams['font.size'] = font_size
        plt.rcParams['font.family'] = font_family
        fig, axs = plt.subplots(2, 1, figsize=(11, 7), gridspec_kw={
                                'height_ratios': [4, 1]}, sharex=True)
        axs[1].set_xlabel('Phase', fontdict=None,
                          labelpad=None, fontname='monospace', size=8)
        axs[0].set_ylabel(
            'Radial velocity [km $s^{-1}$]', fontdict=None, labelpad=None, fontname='monospace', size=8)
        axs[1].set_ylabel('RV residual', fontdict=None,
                          labelpad=None, fontname='monospace', size=8)
        axs[0].grid(color='grey', alpha=0.2, linestyle='-',
                    linewidth=0.5, axis='both', which='both')
        axs[1].grid(color='grey', alpha=0.2, linestyle='-',
                    linewidth=0.5, axis='both', which='both')

        # plot orbital solution
        self._model_x = np.arange(0, 1.011, 0.001)
        self._model_y = list(map(lambda x: self.__computeRadialVelocityCurve(
            x, self._v0, self._orbital_solution[0][1], self._orbital_solution[0][3], self._orbital_solution[0][2], self._orbital_solution[0][0]), self._model_x))
        axs[0].plot(self._model_x, self._model_y, 'k',
                    alpha=0.7, lw=0.7, label='Orbital solution')

        # plot dots
        self.__plotRadialVelocityDots(axs, self._t0)

        split_oname = title.split(' ')
        t = ''.join(r"$\bf{%s}$ " % (w) for w in split_oname)
        p = f'{self._orbital_solution[0][5]} ± {round(self._orbital_solution[1][5],4)} days'
        subtitle = f'{subtitle}\nT0={self._t0} P={p}' if subtitle else f'T0={self._t0} P={p}'
        axs[0].set_title("%s\n%s" % (t, subtitle), fontsize=9,
                         fontweight="0", color='black')

        axs[0].yaxis.set_major_locator(MultipleLocator(rv_y_multiple))
        axs[0].axhline(0, color='black', linewidth=0.7, linestyle="--")

        axs[1].yaxis.set_major_locator(MultipleLocator(residual_y_multiple))
        axs[1].axhline(0, color='black', linewidth=0.7, linestyle="--")

        axs[0].legend(bbox_to_anchor=(1, 1), loc="upper left",
                      frameon=False, prop={'size': 8})
        plt.tight_layout(pad=1, w_pad=0, h_pad=1)
        plt.xticks(np.arange(0, 1.01, 0.1))
        if savefig:
            plt.savefig(
                f'{self._spectra_path}/{self._object_name}_phased_result.png', dpi=dpi)
        plt.show()

    def plotlyRadialVelocityCurve(self, title="", font_family='monospace', font_size=9, show=True, group_by_instrument=True):
        """
        Plot the radial velocity curve using plotly
        # Todo : update parameters and link to yaml config file
        :param title:
        :param font_family:
        :param font_size:
        :param show:
        :return: fig
        """

        if not self._orbital_solution:
            self.__solveSystem()

        fig = go.Figure()

        # plot orbital solution
        self._model_x = np.arange(0, 1.011, 0.001)

        self._model_y = list(map(lambda x: self.__computeRadialVelocityCurve(x, self._v0, self._orbital_solution[0][1],
                                                                             self._orbital_solution[0][3],
                                                                             self._orbital_solution[0][2],
                                                                             self._orbital_solution[0][0]),
                                 self._model_x))
        fig.add_trace(go.Scatter(x=self._model_x, y=self._model_y,
                      mode='lines', name='Orbital solution', line=dict(color='black', width=1)))

        observers = {}
        marker_index = {}
        instruments = {}
        color_number = 0
        period = self._orbital_solution[0][5]

        cmap = plt.get_cmap('tab20')

        # define colors (max 60 distinct observers)
        colors = cmap((np.arange(20)).astype(int), alpha=1)
        + cmap((np.arange(20)).astype(int), alpha=.75)
        + cmap((np.arange(20)).astype(int), alpha=.5)

        markers_style = ['circle', 'square',
                         'diamond', 'triangle-up', 'triangle-down']

        for s in self._sb_spectra:
            # compute phase of the sytem
            jd = s.getJD()
            phase = self.__getPhase(float(self._t0), period, jd)
            s.setPhase(phase)

            # get the observer
            obs = s.getObserver()
            if (obs not in observers.keys()):
                rgb = colors[color_number][:3] * 255
                str_rgb = ",".join([str(rgb[0]), str(rgb[1]), str(rgb[2])])
                observers[obs] = f'rgba({str_rgb}, {colors[color_number][3]})'
                color_number += 1
            color = observers[obs]

            if group_by_instrument:
                # get the instrument
                label = f"{obs} - {s.getInstrument()[:30]}…"
                if label in instruments:
                    fig.add_trace(
                        go.Scatter(x=[phase],
                                   y=[s.getRV()],
                                   mode='markers',
                                   marker_symbol=instruments[label],
                                   marker=dict(color=color,
                                               size=8),
                                   showlegend=False))
                else:
                    if obs not in marker_index:
                        marker_index[obs] = 0
                    instruments[label] = markers_style[marker_index[obs]]
                    marker_index[obs] += 1
                    fig.add_trace(
                        go.Scatter(x=[phase],
                                   y=[s.getRV()],
                                   mode='markers',
                                   name=label,
                                   marker_symbol=instruments[label],
                                   marker=dict(color=color,
                                               size=8),
                                   showlegend=True))
            else:  # no grouping
                # set label to date
                label = f"{obs} - {s.getDate()}"
                fig.add_trace(
                    go.Scatter(x=[phase],
                               y=[s.getRV()],
                               mode='markers',
                               name=label,
                               marker_symbol='circle',
                               marker=dict(color=color,
                                           size=8),
                               showlegend=False))

            # set hover text size
            fig.update_traces(hovertemplate=None,
                              hoverlabel=dict(font_size=16))

            # set hover text config
            fig.update_layout(hovermode="x unified",
                              hoverlabel=dict(bgcolor="white",
                                              font_size=10))

        p = f'{self._orbital_solution[0][5]} ± {round(self._orbital_solution[1][5],4)} days'
        title += f' T0={self._t0} P={p}'

        fig.update_layout(
            title=title,
            plot_bgcolor='white',
            xaxis=dict(
                ticks='outside',
                showline=True,
                tickmode='array',
                tickvals=np.arange(0, 1.01, 0.1),
                gridcolor='lightgrey',
                linecolor='black',
                ticktext=[f'{round(i, 2)}' for i in np.arange(0, 1.01, 0.1)],
                mirror=True
            ),
            xaxis_title="Phase",
            yaxis_title="Radial Velocity (km/s)",
            yaxis=dict(
                ticks='outside',
                showline=True,
                linecolor='black',
                gridcolor='lightgrey',
                mirror=True
            ),
            font=dict(
                family=font_family,
                size=int(font_size)+2,
                color="black",
            )
        )

        # plot
        if show:
            fig.show()

        return fig

    def addFootNote(self, ax, footnote):
        # Add a footnote below and to the right side of the chart
        ax.annotate(footnote,
                    xy=(1.0, -0.2),
                    xycoords='axes fraction',
                    ha='right',
                    va="center",
                    fontsize=8)

    def plotSpec2DFlux(self, title="", subtitle="", savefig=False, dpi=150, font_family='monospace', font_size=9):
        """
        Plot the 2d dynamic spectra
        :param title:
        :param font_family:
        :param font_size:
        :param show:
        :return: fig
        """
        plt.rcParams['font.size'] = font_size
        plt.rcParams['font.family'] = font_family

        if not self._orbital_solution:
            self.__solveSystem()

        # sort sb spectra by phase
        self._sb_spectra.sort(key=lambda x: x.getPhase())

        # create y axis range (phase)
        y_phase = np.arange(0, 1.01, 0.015)
        # create x axis range (wavelength) and convert to km/s
        resample_grid = np.arange(6550, 6576, 0.01)
        sc = SpectralCoord(resample_grid, unit='AA')
        wv_to_kms = sc.to(u.km / u.s, doppler_convention='optical',
                          doppler_rest=6562.82 * u.AA)
        spec2d = np.zeros((len(wv_to_kms), len(y_phase)))

        # resample each spectrum
        for s in self._sb_spectra:
            ss = copy.copy(s)
            # apply heliocentric/barycentric correction
            ss.shift_spectrum_to(
                radial_velocity=s.getRVCorrection())
            fluxc_resample = LinearInterpolatedResampler()
            output_spectrum1D = fluxc_resample(ss, sc)
            phase = s.getPhase()
            indice = int(round(phase, 2) * len(y_phase))
            spec2d[:, indice] = output_spectrum1D.flux

        # prepare imshow
        plt.rcParams["figure.figsize"] = (8, 7)
        fig, ax = plt.subplots()
        ax.imshow(np.rot90(spec2d), extent=[wv_to_kms.min().value,
                                            wv_to_kms.max().value, 0, 1], aspect='auto')

        ax.set_ylabel('Phase', fontdict=None,
                      labelpad=None, fontname='monospace', size=9)
        ax.set_xlabel(
            'Radial velocity [km $s^{-1}$]', fontdict=None, labelpad=None, fontname='monospace', size=9)

        split_oname = title.split(' ')
        t = ''.join(r"$\bf{%s}$ " % (w) for w in split_oname)
        p = f'{self._orbital_solution[0][5]} ± {round(self._orbital_solution[1][5],4)} days'
        subtitle = f'{subtitle}\nT0={self._t0} P={p}' if subtitle else f'T0={self._t0} P={p}'
        ax.set_title("%s\n%s" % (t, subtitle), fontsize=9,
                     fontweight="0", color='black')

        plt.tight_layout(pad=3, w_pad=0, h_pad=1)
        plt.yticks(np.arange(0, 1.01, 0.1))
        if savefig:
            plt.savefig(
                f'{self._spectra_path}/{self._object_name}_2d_spectrum_result.png', dpi=dpi)
        plt.show()

        return spec2d
