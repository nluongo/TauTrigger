# TauTrigger

TauEtStudy.cpp
Read in signal and background data from output_Z80.root and output_MB80.root and plot average energy per cell as well as reconstructed energy vs. true energy. Also create a text file containing information on different cell subset candidates for reconstructed energy and event pre-processing.

TauEtFuncs.cpp
Holds functions called by all other scripts working with tau data

TauAllRecos.cpp
Loop over all possible reconstructed energy definitions to determine which best separates signal and background

TauFCore.cpp
Loops over all possible core and isolation energy definitions and applies the FCore algorithm to determine which best separates signal and background

TauFCoreEtHisto.cpp
Uses given FCore definition to calculate and apply a 95% FCore cut on signal and background and graph efficiencies resulting. Also plots single and di-tau rates with and without the FCore cut.
