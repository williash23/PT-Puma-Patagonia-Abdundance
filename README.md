# PT-Puma-Patagonia-Abundance
All scripts and data required to run two comparative methods for estimating abundance of puma in Torres Del Paine NP in Chile. 

Method 1: Space-to-Event model (STE) from Moeller et al. 2018 fit via bayesian implementation that varies sample window length and generates multiple data sets by shifting overall study start time by 1 second. Only needed script is "run_scenarios_STE.R"

Method 2: Partial genotype partial Identity Model (gPSIM) from Augustine et al. 2020, fit via 'nimble'. Scripts need to be run in order: "1_prep_data_gSPIM.Rdata", "2_prep_nimble_gSPIM.Rdata", "3_run_gSPIM.R", "4_plot_outputs_gSPIM.R".
