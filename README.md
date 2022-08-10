# SimpleModels
Simple climate and cryosphere models in MatLab

MarineIceSheetModel.m:

Analytical steady state solution of a 1d marine ice sheet following Schoof (2007), extended with an analytical profile of a steady-state ice shelf. Each set of parameters results in a steady-state ice sheet geometry. The model can be forced by either changes in sea level (it lowers or rises the bedrock) or mean surface mass balance. The model comes with two bed profiles, a linearly downward sloping bed and a overdeepened bed. The latter exhibits hysteresis. Besides the ice sheet geometry, other panels show ice flux and the phase space between the forcing variable and grounding line position.

MinimalGlacierModel.m:

Minimal glacier model (Oerlemans) on a linearly-sloping bedrock. Forcing is done by changing the equilibrium line altitude (ELA), here defined as Er being the difference in elevation between the highest point of the bedrock profile and the ELA. The slope of the bed can also be controlled. The linear model considers ice thickness derived from the yield stress (perfect plasticity). The nonlinear model makes ice thickness also a function of the length of the glacier. The latter exhibits hysteresis. Besides the glacier profile, plots of length versus slope and length versus Er are provided.

WeertmanIceSheetModel.m:

Plastic ice sheet model based on Weertman (1974). The model is controoled by the size of the accumulation area, the slope of the ELA and the latitudinal position of the ELA. Bedrock adjustment is incorporated. The model exhibits hysteresis as a function of ELA. Both the ice sheet profile and the equilibria can be displayed.

