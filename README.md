# 2019_eLife_Desender
Repo for "Confidence predicts speed-accuracy tradeoff for subsequent decisions", authored by Kobe Desender, Annika Boldt, Tom Verguts, & Tobias H donner. Published in eLife.

The raw data from this paper are available here https://osf.io/83x7c/

#BEHAVIOR
1. "preprocess_and_behavior.R" loads and preprocesses raw data in R 3.4.2, saves the cleaned data for hddm fitting, and creates the figures and analysis for the behavior.

2. "hddm_fit.py" loads the preprocessed data, and fits this using the HDDM (version 0.6.0) in python 2.7.

3. "plot_ddm_output.R" loads the hddm output and makes beautiful plots from these.

#EEG

4. "preprocessing"

5. "hddm_eegRegression_fit.py" loads preproceseed EEG data, and fits a model estimating regression coefficients relating previous trial Pe/ERN to bound/drdift.

5. "hddm_eegBinned_fit.py" loads preproceseed EEG data, and fits a model estimating the influence of previous trial binned amplitude Pe on bound/drift.

6. "Plotting":


