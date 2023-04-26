from spectroscopicbinarysystem import SpectroscopicBinarySystem

sbs = SpectroscopicBinarySystem(
    object_name='hd123299',
    spectra_path='./examples/alphadra/',
    t0=2451441.804,
    period_guess=51,
    debug=True)

sbs.plotRadialVelocityCurve(
    title="α Dra - HD123299 - Phased radial velocities", savefig=True)

# sbs.plotlyRadialVelocityCurve(
#     title="α Dra - HD123299 - Phased radial velocities",
#     group_by_instrument=False)
