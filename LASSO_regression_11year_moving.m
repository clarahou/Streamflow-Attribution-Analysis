%*************************************************************************%
% Description: This code performs LASSO multiple linear regression on the
% changes of parameter n using three explanatory variables: snow fraction (fs),
% leaf area index (LAI), and precipitation intensity (Pi). The regression
% equation is Δn = a1 * Δfs + a2 * ΔLAI + a3 * ΔPi + a0
% Author: Ying Hou
% Date: 24th Jan, 2025
%*************************************************************************%

%% Load data
% Load the data for n, fs, and Pi from 931 catchments, with 25 time steps for each
% matrix represents the 11-year sliding average of each parameter for the period 1982-2016
load('n_11year_sliding', 'n_11year');  % Load n (931 x 25 matrix)
load('fs_11year_sliding', 'fs_11year'); % Load snow fraction (fs) (931 x 25 matrix)
load('LAI_11year_sliding', 'LAI_11year'); % Load leaf area index (LAI) (931 x 25 matrix)
load('Pi_11year_sliding', 'Pi_11year'); % Load precipitation intensity (Pi) (931 x 25 matrix)

Catchment_Number = 931; % Number of catchments

%% Calculate Δn, Δfs, ΔLAI, ΔPi series
% Calculate the differences for each catchment from the 25 time steps
% Δn = n(t+1) - n(t), Δfs = fs(t+1) - fs(t),  ΔLAI = LAI(t+1) - LAI(t), ΔPi = Pi(t+1) - Pi(t)
% These will result in 24 time steps for each catchment
Delta_n = diff(n_11year, 1, 2);  % 931 x 24 matrix
Delta_fs = diff(fs_11year, 1, 2);  % 931 x 24 matrix
Delta_LAI = diff(LAI_11year, 1, 2);  % 931 x 24 matrix
Delta_Pi = diff(Pi_11year, 1, 2);  % 931 x 24 matrix

%% Initialize the matrices to store regression results
a1_all = zeros(Catchment_Number , 1);  % Regression coefficients for Δfs
a2_all = zeros(Catchment_Number , 1);  % Regression coefficients for ΔLAI
a3_all = zeros(Catchment_Number , 1);  % Regression coefficients for ΔPi
a0_all = zeros(Catchment_Number , 1);  % Constant term (a0)
R2_all = zeros(Catchment_Number , 1);  % R^2 for each catchment
pvalue_all = zeros(Catchment_Number , 1);  % p-value for the regression for each catchment

%% Perform LASSO regression
for i = 1:Catchment_Number   % Loop over all catchments
    % Prepare the data for the regression
    X = [Delta_fs(i, :)', Delta_LAI(i, :)', Delta_Pi(i, :)'];  % Independent variables
    Y = Delta_n(i, :)';  % Dependent variable

    % Perform LASSO regression with 5-fold cross-validation to find the optimal lambda
    [B, FitInfo] = lasso(X, Y, 'CV', 5);  % Perform LASSO with 5 cross-validation

    % Choose the optimal lambda based on the smallest MSE
    optimal_lambda_index = FitInfo.IndexMinMSE;  % Index of the best lambda
    optimal_B = B(:, optimal_lambda_index);  % Coefficients for the optimal lambda
    optimal_intercept = FitInfo.Intercept(optimal_lambda_index);  % Optimal intercept
    
    % Store the regression coefficients and statistics
    a1_all(i) = B(1);   % Coefficient for Δfs
    a2_all(i) = B(2);   % Coefficient for ΔLAI
    a3_all(i) = B(3);   % Coefficient for ΔPi
    a0_all(i) = optimal_intercept;  % Constant term (a0)

    % Calculate R² (coefficient of determination)
    Y_pred = X * optimal_B + optimal_intercept;  % Predicted Y
    R2_all(i) = 1 - sum((Y - Y_pred).^2) / sum((Y - mean(Y)).^2);  % R² calculation

    % Perform standard linear regression to compute p-value
    mdl = fitlm(X, Y);  % Fit a linear model
    pvalue_all(i) = mdl.Coefficients.pValue(1);  % Extract p-value for the intercept
end