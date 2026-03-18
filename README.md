# Effects of tree diversity and mycorrhizal type on the spatio-temporal variability of leaf litter production and quality in temperate forests

**Authors:** Elisabeth Bönisch, Benjamin M. Delory, Šárka Angst, Gerrit Angst, Abderrahim Diane, Peter Dietrich, Olga Ferlian, Andreas Fichtner, Rahel Gormanns, Claudia Guimarães-Steinicke, Stephan Hättenschwiler, Julius Quosh, Tama Ray, Anika Walter, Nico Eisenhauer*, Rémy Beugnon*

**Journal:** Journal of Ecology

**Article DOI:** TBA

**Data DOIs:**  

*Leaf litterfall biomass & wood biomass increment*  Bönisch, E., Beugnon, R., Gormanns, R., & Eisenhauer, N. (2026). Monthly species-level leaf litterfall biomass in the MyDiv tree diversity experiment (2023–2024) [Dataset]. MyDiv database. https://doi.org/10.25829/SAHC-C606

*Leaf litterfall C, N, P* Bönisch, E., Beugnon, R., Walter, A., Diane, A., & Eisenhauer, N. (2026). Monthly species-level leaf litter C, N, and P concentrations in the MyDiv tree diversity experiment (2023–2024) [Dataset]. MyDiv database. https://doi.org/10.25829/PH0D-S656

*Precipitation data:* data from DWD (Deutscher Wetterdienst) for Bad Lauchstädt (Wetterstation 02878; 2018-2024, monthly data) https://opendata.dwd.de/climate_environment/CDC/observations_germany/climate/

This folder contains all scripts used to build the models and figures of our paper.

**Analysis scripts:**

*Main results*

Fig_2_1_Main-monthly-cumulative-effects-wood-biomass_inventory_data_kg_ha.R - Testing the Hypothesis I (species-rich tree communities would exhibit increased annual productivity) & Hypothesis II (highest wood biomass increment and leaf litterfall in AM tree communities, intermediate productivity in tree communities of both mycorrhizal types and lowest productivity in pure EM tree communities; regardless of tree species richness) and plotting Fig. 2.

Fig_2_2_TUKEYTEST_Main-monthly-cumulative-effects-wood-biomass_inventory_data.R - Tukey Test between mycorrhiza treatments (AM vs. EM, AM vs. AM+EM, EM vs. AM+EM).

Fig_3_1_Temporal-spatial-variability-annual and monthly-and-number of months of litterfall.R - Testing the Hypothesis III (extended and more variable litterfall over time in species-rich communities containing both mycorrhizal types) & Hypothesis IV (combining higher tree species richness with trees associated with both mycorrhizal types would enhance spatial variation in litter production ) and plotting Fig. 3.

Fig_3_2_TUKEYTEST_spatio-temporal-variability.R - Tukey Test between mycorrhiza treatments (AM vs. EM, AM vs. AM+EM, EM vs. AM+EM).

Fig_3_3_TUKEYTEST_number-of-months-of-litterfall.R - Tukey Test between mycorrhiza treatments (AM vs. EM, AM vs. AM+EM, EM vs. AM+EM).

Fig_4_1_CNP-annual-fluxes-concentration-ratios.R - Testing the Hypothesis V (leaf litterfall quality to increase with the abundance of AM trees and in species rich communities) and plotting Fig. 4.

Fig_4_2_TUKEYTEST-CNP-annual-fluxes-concentration-ratios.R - Tukey Test between mycorrhiza treatments (AM vs. EM, AM vs. AM+EM, EM vs. AM+EM).

Fig_5_1_CNP-monthly-fluxes-concentrations-ratios-cumulative-time-series-plots.R - Plotting Fig. 5.


*Methods*

Increment-calculation.R - Calculation of the wood biomass increment of MyDiv trees from 2022 to 2022 (Table S9, S10, S11).

Main-effects-cross-plot-litterfall.R - Analysis and plotting of cross-plot leaf litterfall (leaf litter from adjacent plots) (Table S6; Fig. S2, S3)

Precipitation_2018to2024.R - Plotting of precipitation patterns at the MyDiv site across years (highlights the drought period) (Fig. S18).

Propagation error from the NIRS-inferred nutrient concentrations.R - Sensitivity analyses (Fig. S13 and S14) with the estimation error per sample (for C, N, and P concentrations and fluxes) and the consequences of the error on the model results (regression for each simulation of C, N, and P concentrations and fluxes).

Sample selection-optimisation-NIRS.R - Selection of leaf litter samples for the analysis with NIRS.

Tree inventory data 2023-dead-trees.R - Dead trees in MyDiv (Table S4).
