# Serial integration of sensory evidence for perceptual decisions and oculomotor responses

This repo contains the code and the data. To run the analyses, first unzip `./data/raw_data.zip`, then run the scripts `build_dataset.R`, `prepare_data.R`, and `prepare_data_latency_split.R` (for the analysis split by saccadic latency). Be prepared that this is quite slow to run (the models are fit in Stan using MCMC sampling). The R-functions folder contain R code to perform 1D cluster test, and the stan code of the Bayesian reverse-correlation analyses with random-walk prior is in the stan-code folder. Finally, to visualize the results use the function that begins with `plot_`. 

The analyses and plotting require several R packages to run: `rstan`, `rethinking`, `ggplot2` and `ggpubr`.

For any question, please contact: m.lisi [at] essex.ac.uk
