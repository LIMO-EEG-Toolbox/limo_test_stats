# PCP simulations

This subfolder contains the scripts for generating the simulated data and analyzing results related to the PCP algorithm.

## Dependencies

Simulations are based on the signal+noise ERP generator of Yeung et al. (2018). https://data.mrc.ox.ac.uk/data-set/simulated-eeg-data-generator.
Results reuse stats functions from https://github.com/CPernet/Robust_Statistical_Toolbox.

## Simulation_PCP

[Simulation_PCP.m](https://github.com/LIMO-EEG-Toolbox/limo_test_stats/blob/master/PCP_simulations/Simulation_PCP.m) is the main scripts that creates the simulated data, classify results and compute distances (Pearson and Kolmogorov-Svirnov). It creates the various data calling the [generate_SNtrials.m](https://github.com/LIMO-EEG-Toolbox/limo_test_stats/blob/master/PCP_simulations/generate_SNtrials.m) function, which wraps around the ERP generator from Yeung et al. (2018). 

[Results.m](https://github.com/LIMO-EEG-Toolbox/limo_test_stats/blob/master/PCP_simulations/Results_PCP.m) performs the analysis. This reproduces what is in the article, loading the csv files.

## Simulation_limits

By default, a PCA is peformed on well conditioned data, i.e. with more trials and time or freqency frames. With real life data, this is not always possible and this scripts thus tests the limits of what can be done varying the number of trials (1500 to 126) and sampling rate (1000Hz, 500Hz, 250Hz) and using white noise outliers. This allowed to compare the same data at different sampling rates with the same trials/frames ratio.

The summary of all the results is on the wiki [here](https://github.com/LIMO-EEG-Toolbox/limo_test_stats/wiki/Defining-outlier-trials-with-PCP)
