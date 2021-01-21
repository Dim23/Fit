# Quick information with step-by-step instruction.

To start the fitting process, run `GammaFit.C` :

    root -b GammaFit.C(fileadres, current_mult, outadres, minNch)

Where the arguments are:

`fileadres` - input root file.

`current_mult` - TH1D Histogram with multiplisity.

`outadres` - output file from fitting process.

`minNch` - lower value of the fitting area ( `25` by default).

### OUTPUT

Resulting file contains TCanvas showing fit results and data-to-fit ratio - `data fit result` .

The resulting fit function of the multiplicity distribution - `fit_func` .

TGraphErrors of impact parametr as a function of centrality - `fit_B_Mean` .
