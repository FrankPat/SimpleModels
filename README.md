# SimpleModels
Simple climate and cryosphere models in MatLab

ClimateFeedbackModel.m

Simple radiation balance model of incoming and outgoing radiation. Controls are the fraction of the solar constant and a control on the emissivity of the system (crude representation of the greenhouse effect). When albedo is made a function of temperature (ice), hysteresis occurs (ice free and frozen Earth). The model displace the energy balance as a function of temperature, Solar constant versus temperature and emissivity versus temperature.

MarineIceSheetModel.m:

Analytical steady state solution of a 1d marine ice sheet following Schoof (2007), extended with an analytical profile of a steady-state ice shelf. Each set of parameters results in a steady-state ice sheet geometry. The model can be forced by either changes in sea level (it lowers or rises the bedrock) or mean surface mass balance. The model comes with two bed profiles, a linearly downward sloping bed and a overdeepened bed. The latter exhibits hysteresis. Besides the ice sheet geometry, other panels show ice flux and the phase space between the forcing variable and grounding line position.

MinimalGlacierModel.m:

Minimal glacier model (Oerlemans) on a linearly-sloping bedrock. Forcing is done by changing the equilibrium line altitude (ELA), here defined as Er being the difference in elevation between the highest point of the bedrock profile and the ELA. The slope of the bed can also be controlled. The linear model considers ice thickness derived from the yield stress (perfect plasticity). The nonlinear model makes ice thickness also a function of the length of the glacier. The latter exhibits hysteresis. Besides the glacier profile, plots of length versus slope and length versus Er are provided.

WeertmanIceSheetModel.m:

Plastic ice sheet model based on Weertman (1974). The model is controoled by the size of the accumulation area, the slope of the ELA and the latitudinal position of the ELA. Bedrock adjustment is incorporated. The model exhibits hysteresis as a function of ELA. Both the ice sheet profile and the equilibria can be displayed.


Simple Climate Model (SCM) and Simple Climate and Sea Level Model (SCM-SLRM)

Climate model emulator based on a basic energy balance model and extended with different feedback parameters.
Simple radiative model based that predicts global mean temperature as a function of CO2 equivalent atmosphere concentrations. The latter is the basic input for the model. CO2-eq. is related to radiative forcing as:

  dQ = 5.35 ln([CO2]/[CO2_o])

where [CO2_o] is the preindustrial level (280 ppmv). The black body temperature is then:

  dTs0 = -dQ/L0, where L0 = -4 \sigma Te^3

where \sigma is the Stefan-Boltzmann constant and Te the temperature at the top of the troposphere (255 K). Atmospheric temperartures are then related to the black body temperature via a feedback factor Ff:

  dTs = Ff . dTs0

where Ff = La0/Laf and Laf is the sum of all feedbacks. Ff ~ 1.74. Ff can also be determined analytically by correlating dTs0 directly to observed temperatures (OptimizeFeedback = 1). Future scenarios rely on the relationship between carbon emissions and CO2-eq. concentrations in the atmosphere, taken from the observed changes in the last 50 years. The relation is automatically calculated.

Datasets used are:
(1) observed CO2-eq. atmosphere concentration, taken into account CO2, CH4, N20 emissions as well as negative forcing effects from aerosols.
(2) Historical CO2 emissions, including land-use changes (Gt CO2)
(3) Historical temperature changes (for optimization of Ff)
(4) Future carbon emissions scenarios (SSPs)

Mean global surface temperatures (yearly). This file contains a brief summary of the changes in Earth's global average surface temperature estimated by combining the Berkeley Earth land-surface temperature field with a reinterpolated version of the HadSST4 ocean temperature field. The current citation for this dataset is: Rohde, R. A. and Hausfather, Z.: The Berkeley Earth Land/Ocean Temperature Record, Earth Syst. Sci. Data, 12, 3469-3479, https://doi.org/10.5194/essd-12-3469-2020, 2020. The dataset differs slightly from the dataset as described in the citation as HadSST3 has been replaced with the newer HadSST4, and associated interpolation parameters have been refit accordingly.  No other changes in methods were needed when moving to the new version of HadSST. http://berkeleyearth.lbl.gov/auto/Global/Land_and_Ocean_summary.txt

Observed trends in total greenhouse gas concentration levels between 1860 and 2019, considering all greenhouse gases and other forcing agents (including aerosols) https://www.eea.europa.eu/ims/atmospheric-greenhouse-gas-concentrations

Global CO2 emissions from fossil fuels https://ourworldindata.org/co2-emissions

SSP emission scenarios: https://www.ipcc.ch/data/ and https://tntcat.iiasa.ac.at/SspDb/dsd?Action=htmlpage&page=10

To convert from gigatonnes carbon to gigatonnes of carbon dioxide, you simply multiply 44 over 12. In other words, 1 gigatonne of carbon equals 3.67 gigatonnes of carbon dioxide.

Sea level rise

This code produces results of Vermeer and Rahmstorf.  "Global Sea Level Linked to Global Temperature", PNAS 2009 For a description of what it does, see that paper. Code written by Stefan Rahmstorf and Martin Vermeer. Please note that this is a scientific code to be used by scientists who know what they are doing; it is not "fool-proof" for the general user in the sense that all combinations of parameter choices, options etc. have been thoroughly tested. Please report any problems or errors to the authors. The code uses the SSAtrend of Moore et al. 2005 for smoothing the data, which can be downloaded at http://www.glaciology.net/software/ssatrend-m (As a simple matlab alternative, fitting a polynomial is available as an option below. This is not recommended; it gives only a poor fit, 
although it has little effect on the future projections.) Reference: J. C. Moore, A. Grinsted, S. Jevrejeva, Eos 86, 226 (2005).

Data: GMSL dataset at CSIRO: https://research.csiro.au/slrwavescoast/sea-level/measurements-and-data/sea-level-data/

lambda: coefficient for linear combination of T and dT/dt (MV)
dH/dt = a T + b dT/dt = a (T + lambda dT/dt), so lambda = b/a

