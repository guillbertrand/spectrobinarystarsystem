
from PIL import Image
from spectroscopicbinarysystem import SpectroscopicBinarySystem
from astropy.coordinates import SkyCoord, EarthLocation, SpectralCoord
from skimage.transform import radon, iradon
import matplotlib.pyplot as plt
from matplotlib import image
# numpy
import numpy as np
import astropy.units as u


def tomography(img):
    resample_grid = np.arange(6550, 6576, 0.01)
    sc = SpectralCoord(resample_grid, unit='AA')
    wv_to_kms = sc.to(u.km / u.s, doppler_convention='optical',
                      doppler_rest=6562.82 * u.AA)

    def crop_center(img, cropx, cropy):
        y, x = img.shape
        startx = x//2 - cropx//2
        starty = y//2 - cropy//2
        return img[starty:starty+cropy, startx:startx+cropx]
    rr = iradon(img)
    plt.rcParams["figure.figsize"] = (8, 7)
    fig, ax = plt.subplots()
    ax.imshow(crop_center(rr, 600, 600), norm="linear", aspect='auto')
    plt.show()


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
    debug=False)

# sbs.plotRadialVelocityCurve(
#     title="α Dra - HD123299 - Phased radial velocities",
#     subtitle=f"{sbs.getObservationCount()} observations collected from april 2022 to april 2023\nhttps://alphadra.staros-projects.org/\n",
#     savefig=True)

# spec2d = sbs.plotSpec2DFlux(
#     title="α Dra - HD123299 - Hα line 2d dynamic spectrum",
#     subtitle=f"{sbs.getObservationCount()} observations collected from april 2022 to april 2023\nhttps://alphadra.staros-projects.org/\n",
#     savefig=True
# )

# tomography(spec2d)

srcImage = Image.open("examples/mizar.jpg")
grayImage = srcImage.convert('L')
array = np.array(grayImage)


def crop_center(img, cropx, cropy):
    y, x = img.shape
    startx = x//2 - cropx//2
    starty = y//2 - cropy//2
    return img[starty:starty+cropy, startx:startx+cropx]


rr = iradon(array)
plt.rcParams["figure.figsize"] = (8, 7)
fig, ax = plt.subplots()
ax.imshow(crop_center(rr, 300, 300), norm="linear", aspect='auto')
plt.show()
