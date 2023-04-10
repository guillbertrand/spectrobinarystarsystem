import logging
import re
import os
import math
import argparse
import yaml
import warnings

# numpy
import numpy as np

# astropy
from astropy.io import fits
from astropy.time import Time
import astropy.units as u
import astropy.wcs as fitswcs 
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation
from astropy.modeling import models

# astroquery
from astroquery.simbad import Simbad

# specutils
from specutils import Spectrum1D,SpectralRegion
from specutils.manipulation import noise_region_uncertainty, extract_region
from specutils.analysis import centroid
from specutils.fitting import fit_lines, fit_generic_continuum

# matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)


# orbitalpy
from orbital import utilities

# binarystarsolve
from binarystarsolve.binarystarsolve import StarSolve

#

__version__ = 0.5
__debug_mode__ = 0

# 

def getPhase(jd0, period, jd):
    return (jd-jd0) / period % 1 

#

def extractObservations(specs, period = None):
    sc = None
    result_table = Simbad.query_object(conf['object_name'])
    if result_table:
            ra = result_table[0]['RA']
            dec = result_table[0]['DEC']
            sc = SkyCoord(ra+' '+dec, unit=(u.hourangle, u.deg))

    if __debug_mode__:
        plt.rcParams['font.size'] = '6'
        plt.rcParams['font.family'] = 'monospace'
        fig, axs = plt.subplots(5,5, figsize=(11,7), sharex=True, sharey=True)
    
    obs = {}
    i = 0
    for s in specs:
        logging.info('\U0001F4C8 Process spectrum %s'% os.path.basename(s))
        f = fits.open(s)
        header = f[0].header
        logging.info('      - Observation date : %s - %s'% (header['DATE-OBS'], header['JD-OBS']))
        rv = None
        if (sc):
            ax = None
            if __debug_mode__:
                ax = axs.flat[i]
            if('BSS_VHEL' in header and header['BSS_VHEL']):
                rvc = None
            else:
                rvc = radialVelocityCorrection(sc, header['JD-OBS'], header['GEO_LONG'], header['GEO_LAT'], header['GEO_ELEV'])
            logging.info('      - Radial velocity correction (%s) : %s '% (conf['radial_velocity_correction'],rvc))
            center = findCenterOfLine(f, ax, header['CDELT1'])
            rv = getRadialVelocity(center, rvc, __lambda_ref__) 
            logging.info('      - Center of line : %s ± %s'% (center[0], center[1]))
            logging.info('      - Radial velocity : %s ± %s'% rv)

            if __debug_mode__:
                ax.set_title('%s\nJD=%s  \n%s' % (os.path.basename(s), header['JD-OBS'],header['OBSERVER']), fontsize="7")
                ax.grid(True)
                ax.tick_params(axis='both', which='major', labelsize=6)
                ax.tick_params(axis='both', which='minor', labelsize=6)
        

        obs[float(header['JD-OBS'])] = {'fits':s, 'radial_velocity_corr':rvc, 'centroid': centroid, 'radial_velocity':rv, 'header':header }
        i+=1
    
    if __debug_mode__:
        plt.tight_layout(pad=0.8, w_pad=1.5, h_pad=1)
        plt.savefig(wdir+'/bss_debug_result.png', dpi=conf['dpi'])
        plt.show()
    return obs

#

def getRadialVelocity(position, radial_velocity_correction = 0, lambda0 = 6562.82 * u.AA):
    c = 299792.458 * u.km / u.s
    return  ((c  * ((position[0]-lambda0)/lambda0)) + radial_velocity_correction , (c * (position[1])/lambda0.value))

#

def computeRadialVelocityCurve(t, t0, K, e, w, v0):
    w = math.radians(w)
    # Mean anomaly
    M = 2 * np.pi *  ((t - t0) %1)    
    # Eccentric anomaly
    E = utilities.eccentric_anomaly_from_mean(e,M, tolerance=0.00001)   
    # True anomaly
    f = utilities.true_anomaly_from_eccentric(e, E) 
    return (K * (e * np.cos(w) + np.cos(w + f)) + v0 ) 

#

def radialVelocityCorrection(skycoord, jd, longitude, latitude, elevation):
    t = Time(jd, format='jd', scale='utc')
    loc = EarthLocation(longitude,latitude,elevation*u.m)
    vcorr = skycoord.radial_velocity_correction(kind=conf['radial_velocity_correction'], obstime=t, location=loc)  
    return vcorr.to(u.km / u.s)

def findNearest(array, value):
    array = np.asarray(array)
    return (np.abs(array - value)).argmin()

def findCenterOfLine(spectrum,ax,dispersion):
    specdata = spectrum[0].data
    header = spectrum[0].header
    with warnings.catch_warnings():  # Ignore warnings
        warnings.simplefilter('ignore')

        wcs_data = fitswcs.WCS(header={'CDELT1': header['CDELT1'], 'CRVAL1': header['CRVAL1'],
                                    'CUNIT1': 'Angstroms', 'CTYPE1': 'WAVE',
                                    'CRPIX1': header['CRPIX1']})
        flux= specdata * u.Jy
        s = Spectrum1D(flux=flux, wcs=wcs_data)
        
 
        # [2:-2] > Hack to excude first and last value of the spectra.
        # Sometimes fist value for the intensity is equal to 0.00 and prevents finding the true minima.
        # Ex : sample/alphadra/_alphadra_20230406_856.fits
        ipeak = s.flux[2:-2].argmin() 
        xpeak = s.spectral_axis[ipeak+2].to(u.AA)

        w1 = float(conf['spectral_region']) / 2
        w2 = float(conf['window_width']) /2
        fwhm = float(conf['fwhm']) 

        invert_s = extract_region(Spectrum1D(flux=s.flux*-1, wcs=wcs_data),SpectralRegion(xpeak-(w1*u.AA), xpeak+(w1*u.AA)))
        s_fit = fit_generic_continuum(invert_s, exclude_regions=[SpectralRegion(xpeak-(w2*u.AA), xpeak+(w2*u.AA))])
        invert_s   = invert_s / s_fit(invert_s.spectral_axis)
        invert_s = invert_s - 1

        match conf['model']:
            case 'gaussian':
                g_init = models.Gaussian1D(mean=xpeak, amplitude=invert_s.flux.argmax())
            case 'voigt':
                g_init = models.Voigt1D(x_0=xpeak, amplitude_L=2/(np.pi*fwhm))
            case 'lorentz':
                g_init = models.Lorentz1D(x_0=xpeak, fwhm=fwhm)
            case _:
                g_init = models.Gaussian1D(mean=xpeak, amplitude=invert_s.flux.argmax())
        #
        
        g_fit = fit_lines(invert_s, g_init,window=SpectralRegion(xpeak-(w2*u.AA), xpeak+(w2*u.AA)))
        y_fit = g_fit(invert_s.spectral_axis) 

        match conf['model']:
            case 'gaussian':
                center = g_fit.mean
            case _:
                center = g_fit.x_0
        #
        if __debug_mode__:
            ax.plot(invert_s.spectral_axis.to(u.AA), invert_s.flux, color="k")
            ax.plot(invert_s.spectral_axis.to(u.AA), y_fit, color="r")
            ax.axvline(x=center.value, color='r', linestyle='-',lw=0.7)
     
        region = SpectralRegion(6555*u.AA, 6570*u.AA)
 
        #error = (299792.458*conf['precision']*dispersion/__lambda_ref__.value) 
        return (center.value*u.AA, 0, extract_region(s,region))
    
def initPlot():
    plt.rcParams['font.size'] = conf['font_size']
    plt.rcParams['font.family'] = conf['font_family']
    fig, axs =  plt.subplots(2,1,figsize=(conf['fig_size_x'],conf['fig_size_y']),gridspec_kw={'height_ratios': [4, 1]}, sharex=True)
    axs[1].set_xlabel('Phase', fontdict=None, labelpad=None, fontname = 'monospace',size=8)
    axs[0].set_ylabel('Radial velocity [km $s^{-1}$]', fontdict=None, labelpad=None, fontname = 'monospace',size=8)
    axs[1].set_ylabel('RV residual', fontdict=None, labelpad=None, fontname = 'monospace',size=8)
    axs[0].grid(color='grey', alpha=0.2, linestyle='-', linewidth=0.5, axis='both', which='both')
    axs[1].grid(color='grey', alpha=0.2, linestyle='-', linewidth=0.5, axis='both', which='both')
    return (fig, axs)

def plotRadialVelocityCurve(ax, v0, K, e, w, jd0,color="red", lw=0.5, alpha=1, label=""):
    model_x = np.arange(0,1.011, 0.001)
    model_y = list(map(lambda x: computeRadialVelocityCurve(x,jd0,K,e,w,v0), model_x))
    ax.plot(model_x, model_y, color, alpha=alpha, lw=lw, label=label)
    return (model_x, model_y)

def plotRadialVelocityDotsFromData(specs, period, jd0, error, axs, model):  
    if(period):
        for jd in specs:
            phase = getPhase(float(jd0), period, float(jd))
            logging.info('%s phase : %s' % (os.path.basename(specs[jd]['fits']), phase))
            specs[jd]['phase'] = phase
    observers = {}
    i = 0
    for jd, s in specs.items():
        obs = s['header']['OBSERVER'].lower()   
        if(obs not in observers.keys()):
            c = [x.strip() for x in conf["markers_colors"].split(',')][i]
            observers[obs] = {'color': c, 'instruments':{}}
            i+=1
    
    for jd, s in specs.items():
        obs = s['header']['OBSERVER'].lower()
        label = "%s - %s…" % (obs, s['header']['BSS_INST'][0:30])
        if(label not in observers[obs]['instruments'].keys()):
            observers[obs]['instruments'][label] =  [x.strip() for x in conf["markers_styles"].split(',')][len(observers[obs]['instruments'])]
            axs[0].errorbar(s['phase'], s['radial_velocity'][0].value,yerr = 0, label= label, ecolor='k', capsize=0,fmt =observers[obs]['instruments'][label], color=observers[obs]['color'], lw=0.7)
        else:
            axs[0].errorbar(s['phase'], s['radial_velocity'][0].value,yerr = 0, fmt =observers[obs]['instruments'][label], ecolor='k', capsize=0,color=observers[obs]['color'], lw=.7)    
        print(s['header']['DATE-OBS'], observers[obs]['instruments'][label])
        xindex = findNearest(model[0], s['phase'])
        axs[1].errorbar(s['phase'], s['radial_velocity'][0].value- model[1][xindex],yerr = 0, fmt =observers[obs]['instruments'][label], ecolor='k', capsize=0,color=observers[obs]['color'], lw=.7)            
    

def saveAndShowPlot(ax, t0, p):
    t = ''
    split_oname = conf['title'].split(' ')
    for w in split_oname:
        t += r"$\bf{%s}$ " % (w)
    
    ax[0].set_title("%s\n%s\nT0=%s P=%s" % (t,conf['subtitle'],t0,p),fontsize=conf['title_font_size'],fontweight="0", color='black' )

    ax[0].yaxis.set_major_locator(MultipleLocator(conf['fig_rv_y_multiple'])) 
    ax[0].axhline(0, color='black', linewidth=0.7, linestyle="--")

    ax[1].yaxis.set_major_locator(MultipleLocator(conf['fig_residual_y_multiple']))
    ax[1].axhline(0, color='black', linewidth=0.7, linestyle="--")
    
    ax[0].legend(bbox_to_anchor=(1, 1), loc="upper left", frameon=False,prop={'size': 8})
    plt.tight_layout(pad=1, w_pad=0, h_pad=1)
    plt.xticks(np.arange(0, 1.01, 0.1))
    plt.savefig(wdir+'/bss_phased_result.png', dpi=conf['dpi'])
    plt.show()  

if __name__ == '__main__':
    #
    default_conf_filename = 'bss.config.yaml'

    FORMAT = '%(message)s'
    logging.basicConfig(level=logging.INFO, format=FORMAT)

    logging.info('\U0001F680 BinaryStarSystem %s - Start \U0001F680' % __version__)

    #
    parser = argparse.ArgumentParser()
    parser.add_argument("path", type=str, help="Path to your fits directory")
    parser.add_argument("-c", "--config", type=str, default=default_conf_filename, help="Custom config file name")
    
    args = parser.parse_args()
    path = args.path
    conf_filename = args.config

    wdir = path

    FORMAT = '- %(message)s'
    logging.basicConfig(level=logging.INFO, format=FORMAT)

    # load yaml configuration file
    confpath = os.path.join(wdir,conf_filename)
    cust_confpath = os.path.join(wdir, conf_filename)
    if not (os.path.exists(confpath)):
        confpath = conf_filename
        if not (os.path.exists(confpath)):
            confpath = None

    if confpath:
        logging.info('\U00002728 Load configuration file \U0001F527  %s' % (confpath))
        with open(confpath, 'r', encoding='utf8') as f:
            conf = yaml.load(f,Loader=yaml.FullLoader)
    else :
        logging.info('\U0001F4C1 Error : %s not found !' % default_conf_filename)

    __debug_mode__ = conf['debug_mode']
    __lambda_ref__ = float(conf['lambda_ref']) * u.AA

    # find spec files 
    specs = []
    for root, dirs, files in os.walk(wdir):
        for file in files:
            regex = re.compile(conf['spec_file_regex'])
            if(re.match(regex, file)):
                specs.append(os.path.join(wdir, file))
    if not len(specs):
        logging.info('\U0001F4C1 Error : 0 spectrum file found !')
    else:
        logging.info('\U0001F4C1 %d spectra files found !' % (len(specs)))

        p = conf['period_guess']
        if(conf['period']):
            p = conf['period']
        data = extractObservations(specs, p)
        
        # write bss result file for BinaryStarSolver
        with open(wdir+'/bss_results.txt', 'w') as f: 
            for key, value in data.items():
                output = '%s %s' % (float(key)-2400000., round(value['radial_velocity'][0].value,3))
                f.write(output+'\n')

        
        #[γ, K, ω, e, T0, P, a, f(M)]
        params, err, cov = StarSolve(data_file = wdir+"/bss_results.txt", star = "primary", Period= conf['period'], Pguess=conf['period_guess'], covariance = True, graphs=False)
        
        # If conf['T0'] Compute phase delta between T0 of the model and the fixed value 
        if(conf['T0']):
            t0 = conf['T0'] 
            phase1 =  getPhase(t0,params[5],params[4]+2400000)
            phase2 =  1.0
            v0 = phase1 - phase2 
        # Else use t0 compute by the model
        else:
            t0 = params[4] + 2400000
            v0 = 0

        # plot results
        (fig, axs) = initPlot()
        #plotRadialVelocityCurve(axs[0], params[0], 48.512, 0.4229, 21.28,1,'k--', 0.8, 0.8, 'K.Pavlovski, et al. 2021')
        model = plotRadialVelocityCurve(axs[0], params[0], params[1], params[3], params[2], v0, conf['line_color'], 0.8, 0.8, 'Orbital solution')
        plotRadialVelocityDotsFromData(data, params[5], t0, err, axs, model)
        saveAndShowPlot(axs,t0,'%s±%s days' % (params[5],round(err[5],4)))
        print('[γ, K, ω, e, T0, P, a, f(M)]')
        print(params)
        print(err)


    
    