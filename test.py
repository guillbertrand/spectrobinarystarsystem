from spectroscopicbinarysystem.spectroscopicbinarysystem import SpectroscopicBinarySystem

sbs = SpectroscopicBinarySystem('hd123299','./sample/alphadra/', t0=2451441.804, period_guess=51, debug=True)
#sbs.plotRadialVelocityCurve(title="α Dra - HD123299 - Phased radial velocities",savefig=True)
sbs.plotlyRadialVelocityCurve(title="α Dra - HD123299 - Phased radial velocities")
