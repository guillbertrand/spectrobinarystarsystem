from spectroscopicbinarysystem import SpectroscopicBinarySystem, printPhaseEphem

sbs = SpectroscopicBinarySystem(
    object_name='Sigma Gem',
    spectra_path='C:\\Users\\g-ber\\Documents\\ASTRO\\starex\\studies\\siggem',
    t0=2450388.36853,
    period_guess=19.60447, #Duemmler et al. (1997)
    conf={
        "LAMBDA_REF": 6562.82,
        "LINE_FIT_MODEL": "gaussian",
        "LINE_FIT_FWHM": 3.0,
        "RV_CORR_TYPE": "barycentric",
        "SB_TYPE": 1
    },
    verbose=True,
    debug=True)

# print basic phase ephemeris for the next 20 days
# printPhaseEphem(jd0=2451441.804, period=51.4719,
#                 start_date='2023-06-24T23:00:00')

# # plot result with matplotlib and save the results
sbs.plotRadialVelocityCurve(
    title="Sigma Gem - Phased radial velocities",
    subtitle=f"{sbs.getObservationCount()} spectra\n",
    savefig=True,
    residual_y_multiple=2
)

# sbs.plotSTAROSRadialVelocityCurve(
#     savefig=True,
#     color="k",
#     linecolor="k",
#     residual_y_multiple=2
# )

# obscode = { 'a.maetz':'Alain Maetz',
#            'amaetz':'Alain Maetz',
#            'e.bryssinck':'Erik Bryssinck',
#            'f.weil':'Franck Weil ',
#            'jj broussat':'Jean-Jacques Broussat',
#            'jlecomte':'Julien Lecomte',
#            'leadbeater':'Robin Leadbeater ',
#            'lmn':'Magnus Larsson',
#            'sdv':'Simon De Visscher',
#            'simon de visscher':'Simon De Visscher',
#            'vdesnoux':'Valerie Desnoux',
#            'a. garrigos':'Antonio Garrigos',
#            'albert stiewing':'Albert Stiewing',
#            'cbuil':'Christian Buil',
#            'ebe':'Etienne Bertrand',
#            'mlelain':'Matthieu Le Lain',
#            'olivier garde':'Olivier Garde',
#            'xdu':'Xavier Dupont',
#            'aleduc':'Arthur Leduc',
#            'gbertrand':'Guillaume Bertrand',
#            'rdiz':'Rick Diz',
#            'v. lecocq':'Vincent Lecocq',
#            'll. rib?':'Luis Ribe de Pont'
#            }

# for o in sbs.getObservers():
#     sbs.plotRadialVelocityCurve(
#         title=f"α Dra - HD123299 - {obscode[o]}",
#         subtitle="",
#         savefig=True,
#         observer=o,
#         figsize=(6, 5))

# plot 2d dynamic spectrum
# sbs.plotSpec2DFlux(
#     title="α Dra - HD123299 - Hα line 2d dynamic spectrum",
#     subtitle=f"{sbs.getObservationCount()} observations collected from april 2022 to july 2023\nhttps://alphadra.staros-projects.org/\n",
#     savefig=True
# )
sbs.plotSTAROSSpec2DFlux(
    title="Sig Gem - Hα line 2d dynamic spectrum",
    subtitle=f"{sbs.getObservationCount()} observations collected from april 2022 to july 2023\nhttps://alphadra.staros-projects.org/\n",
    savefig=True
)

# # display result with plotly
# sbs.plotlyRadialVelocityCurve(
#      title="α Dra - HD123299 - Phased radial velocities",
#      group_by_instrument=False)
