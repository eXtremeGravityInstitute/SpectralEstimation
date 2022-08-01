# SpectralEstimation
Low latency spectral estimation for LIGO/Virgo/Kagra

The package includes codes for reading LIGO/Virgo data and performing low latency spectral estimation using Bayesian inference. There are two ways to access the data. For members of the LVK collaboration with access top the data grid, there is a code called "getspec.c" that finds and reads data frames, then launches the SpecFit code that performs the spectral estimation. For those outside the LVK collaboration there is the code called "getspecgwosc.c" that processes the publically available data from the GWOSC website https://www.gw-openscience.org.

The core algorithm, SpecFit.c, is based on a variant of the approach described in the paper https://journals.aps.org/prd/abstract/10.1103/PhysRevD.103.104057 (e-Print: 2101.01188).

On the LVK data grid the codes can be compile using the Makefile. The full analysis can then be run via a command line call. For example, to get an 8 second PSD for the Hanford detector at the time of GW150914, with the segment extending from -6 to +2 seconds from the trigger time, execute the command

./getspec 8 1126259462 2 2048 H C02 0

Here the arguments are observation length (s), GPS trigger time, time from trigger to end of segment (s), Nyquist frequency (Hz), detector label, frame type, cleaned or not cleaned. 

