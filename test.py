from spectroscopicbinarysystem import SpectroscopicBinarySystem

sbs = SpectroscopicBinarySystem(
    object_name='hd123299',
    spectra_path='./examples/alphadra/',
    t0=2451441.804,
    period_guess=51,
    conf={
        "LAMBDA_REF": 6562.82,
        "LINE_FIT_MODEL": "gaussian",
        "LINE_FIT_FWHM": 1.0,
        "RV_CORR_TYPE": "barycentric",
        "SB_TYPE": 1
    },
    verbose=True,
    debug=True)

# # plot result with matplotlib and save the results
sbs.plotRadialVelocityCurve(
    title="α Dra - HD123299 - Phased radial velocities",
    subtitle=f"{sbs.getObservationCount()} observations collected from april 2022 to may 2023\nhttps://alphadra.staros-projects.org/\n",
    savefig=True,
    residual_y_multiple=2)


# # plot 2d dynamic spectrum
sbs.plotSpec2DFlux(
    title="α Dra - HD123299 - Hα line 2d dynamic spectrum",
    subtitle=f"{sbs.getObservationCount()} observations collected from april 2022 to may 2023\nhttps://alphadra.staros-projects.org/\n",
    savefig=True
)

# # display result with plotly
# sbs.plotlyRadialVelocityCurve(
#     title="α Dra - HD123299 - Phased radial velocities",
#     group_by_instrument=False)
