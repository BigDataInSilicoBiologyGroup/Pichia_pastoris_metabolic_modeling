# Pichia_pastoris_metabolic_modeling
Script for Matlab for uprgulation and downregulation purposes for Product of intereset increased production

Script to calculate upregulation, downregulation and deletion candidates to upregulate high value chemical bioproduction. 

Before use the script user must initialize CobraToolbox 3.0 toolbox "InitCobraToolbox", because script uses toolbox functionality

Script imput parameters :
1) Filename_model - the metabolic model file name and path (example: MODEL1612130000_v3_Glycerol_version_4_methanol_cultivation.xls)
2) Objective_calc - the high value chemical ID from metabolic model (usually an exchange reaction)
3) Biomass_name - the biomass reaction name ID from metabolic model.

Script output parameters:
1) Table_positive_contra_proportional - reactions with flux value from left to right and are inversely proportionally to forced product (heme) changing fluxes;
2) Table_negative_contra_proportional - reactions with flux value from right to left and are inversely proportionally to forced product (heme) changing fluxes;
3) Table_negative_directly_proportional - reactions with flux value from left to right and are directly proportionally to forced product (heme) changing fluxes;
4) Table_positive_directly_proportional - reactions with flux value from right to left and are directly proportionally to forced product (heme) changing fluxes.

Additionally script creates "Optimisation_result.xlsx" file with tables, which names can be found in Supplementary materials 2.

Script uses CobraToolbox 3.0 functionality thus there is need to pre-install CobraToolbox 3.0 from https://opencobra.github.io/cobratoolbox/stable/

Author : Agris Pentjuss
license : GNU public license 3.0
