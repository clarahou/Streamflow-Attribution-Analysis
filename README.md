# Streamflow-Attribution-Analysis
This repository contains MATLAB scripts for streamflow attribution analysis based on the Budyko framework and LASSO regression. The code is designed to calculate key parameters, perform multivariate linear regression, and evaluate streamflow elasticity to environmental variables. These analyses focus on understanding the impacts of precipitation, potential evapotranspiration, snowfall, vegetation cover, and precipitation intensity on streamflow across catchments.

Author: Ying Hou

Included Scripts

Budyko_11year_moving.m
Description: This script calculates the Budyko model parameter n using Choudhury’s formulation over 11-year sliding windows for each of 931 catchments.
Key Steps:
Load annual precipitation, streamflow, and evapotranspiration data (1982–2016).
Apply an 11-year sliding window to calculate the multi-year water balance.
Iteratively solve for parameter n for each catchment and sliding window.
Outputs:
A matrix (n_11year) containing parameter n for each catchment across all sliding windows.

LASSO_regression_11year_moving.m
Description: This script uses the Least Absolute Shrinkage and Selection Operator (LASSO) regression to evaluate how changes in snowfall fraction (fs), leaf area index (LAI), and precipitation intensity (Pi) contribute to changes in the Budyko parameter n.
Key Steps:
Calculate differences (Δ) for each variable across sliding windows.
Perform LASSO regression for each catchment with five-fold cross-validation to select the optimal regularization parameter (lambda).
Extract regression coefficients (a1, a2, a3, a0), coefficient of determination (R²), and statistical significance (p-value) for each catchment.
Outputs:
Regression coefficients for each variable and catchment.
R² and p-values for the regression models.

Elasticity_analysis.m
Description: This script computes the elasticity of streamflow (Q) to various environmental variables, including precipitation (P), potential evapotranspiration (E0), snowfall fraction (fs), leaf area index (LAI), and precipitation intensity (Pi).
Key Steps:
Use regression coefficients from LASSO and the Budyko model to derive sensitivities.
Convert sensitivities into elasticities for each variable.
Compute contributions of each variable to streamflow changes during the study period.
Outputs:
Sensitivities and elasticities of Q to each variable for each catchment.
Quantified contributions of each variable to streamflow changes.

How to Use

Preparation:
Ensure MATLAB is installed with the necessary toolboxes (e.g., Statistics and Machine Learning Toolbox for LASSO regression).
Place the data files (e.g., n_11year_sliding, fs_11year_sliding) in the appropriate directory.

Run the Scripts:
Start with BudykoParameterCalculation.m to generate n_11year.
Use LASSORegression.m to evaluate the relationships between n and environmental variables.
Finally, run ElasticityCalculation.m to calculate the elasticities and contributions.

Data Requirements:
The scripts assume the input data is formatted as matrices with rows representing catchments and columns representing time.

Key Applications
Catchment Hydrology: Understand how environmental factors drive hydrological partitioning in cold regions.
Sensitivity Analysis: Quantify the relative importance of precipitation intensity, snowfall, vegetation cover, and climate variables on streamflow.
Impact Assessment: Evaluate the impacts of climate change on water resources in snow-dominated catchments.

Citation:
Hou, Y. (2025). Streamflow-Attribution-Analysis [Dataset]. Github. https://github.com/clarahou/Streamflow-Attribution-Analysis

References:
The authors will add appropriate references later.
