from spectroscopicbinarysystem.spectroscopicbinarysystem import SpectroscopicBinarySystem

sbs = SpectroscopicBinarySystem('hd123299','./sample/alphadra/', t0=2451441.804, period_guess=51, debug=True)
sbs.getOrbitalSolution()

title = "α Dra - HD123299 - Phased radial velocities"
subtitle = "April 2022 to April 2023 - Star'Ex HR (2400 l/mm)"

sbs.plotRadialVelocityCurve(title,subtitle,savefig=True)
