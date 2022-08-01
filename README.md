# SpectralEstimation
Low latency spectral estimation for LIGO/Virgo/Kagra

The package includes codes for reading LIGO/Virgo data and performing low latency spectral estimation using Bayesian spectral estimation (with a simple MCMC). There are two ways to import data. For members of the LVK collaboration with access to the LVK data grid, there is a code called "getspec.c" that finds and reads data frames, then launches the SpecFit code that performs the spectral estimation. For those outside the LVK collaboration there is the code called "getspec_gwosc.c" that processes the publically available data from the GWOSC website https://www.gw-openscience.org.

The core algorithm, SpecFit.c, is based on a variant of the approach described in the paper https://journals.aps.org/prd/abstract/10.1103/PhysRevD.103.104057 (e-Print: 2101.01188).

On the LVK data grid the codes can be compile using Makefile_Grid. The full analysis can then be run via a command line call. For example, to get an 8 second PSD for the Hanford detector at the time of GW150914, with the segment extending from -6 to +2 seconds from the trigger time, execute the command

./getspec 8 1126259462 2 2048 H C02 0

Here the arguments are observation length (s), GPS trigger time, time from trigger to end of segment (s), Nyquist frequency (Hz), detector label, frame type, cleaned or not cleaned. 

To perforn the same analysis using the GWOSC data, compile the codes using either Makefile_GWOSC_OSXor Makefile_GWOSC_linux, then download the relevant data file:

curl --output H-H1_GWOSC_16KHZ_R1-1126259447-32.txt.gz https://www.gw-openscience.org/eventapi/html/GWTC-1-confident/GW150914/v3/H-H1_GWOSC_16KHZ_R1-1126259447-32.txt.gz

Unpack:

gunzip H-H1_GWOSC_16KHZ_R1-1126259447-32.txt.gz

Then run the analysis:

./getspec_gwosc H-H1_GWOSC_16KHZ_R1-1126259447-32.txt 8 1126259462 2 2048 H

Both versions of the code produce a point estimate of the PSD (the mean PSD after burn-in) and a copy of the cleaned and whitend data. Note that the cleaning will remove any glitches *and* any GW signals. The cleaned/whitend data is used for checking the Gaussianity of the residuals. This gets done automatically with a call to the code Gauss_Test_Freq.c. The default setting is to check the Gaussianity using 8 Hz wide windows. A gnuscript is then called to make a plot of the resdiuals, along with the PSD and the results of the Anderson-Darling Gaussianity test. Points above the red line have failed the test at 95% confidence. For the examples shown above, we expect 12 or 13 regions to fail the test even if the data is perfectly Gaussian. 







