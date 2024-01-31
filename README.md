Author: Grace Yu, Ting Li, Josh Speagle, Gustavo Medina (U of Toronto)

<ting.li@astro.utoronto.ca>

Code developed for Yu, Li, Speagle et al. 2024.

Notebook1: Given the DES dataset
1)Process the data (by making simple cuts)
2)develop a mixture model (hierarchical bayesian model) using the photometries  (Dynesty Sampling takes around 5 hours using 100 CPUs)
3)output probabilities of BHB, BS, and outlier for each object.

INPUT: "DES_DR2_hot_stars.fits" (A datafile contains dereddened photometry from DES DR2 with blue star sample. This file is not available on github as it is too big for GitHub. Please contact Ting for the file if you need, but you could also download from DES or NOIRLab Datalab.)

OUTPUT: "processed_DES.fits" (from step 2), "data_with_probability.fits" (from step 3), "mixture_model_coefficients.npy" (from step 3)


Notebook2: Given the probability of BHB, BS, outlier for each star and probability threshold
1) select BHB candidates that meet the threshold
2) Make region cut to remove overdensities
3) Derive a density profile by inhomogeneous Poisson point process (Dynesty Sampling takes around 1 minute)
4) Compare with other literature values

INPUT:"data_with_probability.fits"

OUTPUT: "BHB_catalog.fits" (from step 1), "entire_catalog.csv" (which contains all sources with their BHB probability and distances, but with fewer columns and the numbers are rounded, same as Table 3 in the paper.)

NOTE: 
1) transform_coordinates_to_sgr.py is useful to obtain the Sgr coordinate, required in Notebook 2 if we want to remove overdensities in Sgr Stream.
2) In Notebook 2, some steps in deriving the density profile require the computation of an analytically-intractable integral. Hence I precomputed the values at certain points, which are stored in precomputed_integral_values.txt. I used these points to make interpolation to deal with that integral in Notebook 2 when deriving the density profile. precompute_integral.ipynb shows how I computed it. You don't have to use it if your integral is friendly to deal with.
