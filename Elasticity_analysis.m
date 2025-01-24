%*************************************************************************%
% Description: This script calculates the sensitivity and elasticity of 
% streamflow (Q) to precipitation (P), potential evapotranspiration (E0), 
% snowfall fraction (fs), leaf area index (LAI), and precipitation intensity (Pi). 
% The calculations are based on Budyko model and regression coefficients 
% obtained through previous LASSO regression. 
% Author: Ying Hou
% Date: 24th Jan, 2025
%*************************************************************************%

%% Load Data
% Regression coefficients from LASSO regression
load('a1_all','a1_all');  % Regression coefficients for ¦¤fs
load('a2_all','a2_all'); % Regression coefficients for ¦¤LAI
load('a3_all','a3_all');  % Regression coefficients for ¦¤Pi
load('a0_all','a0_all');  % Constant term (a0)
% Catchment average values for variables
load('P_aveg', 'P_aveg');  % Average precipitation
load('E0_aveg', 'E0_aveg');  % Average potential evapotranspiration
load('Qobs_aveg', 'Qobs_aveg');  % Average streamflow
load('n_aveg', 'n_aveg');  % Average Budyko parameter n
load('fs_aveg', 'fs_aveg');  % Average snowfall fraction
load('LAI_aveg', 'LAI_aveg');  % Average leaf area index
load('Pi_aveg', 'Pi_aveg');  % Average precipitation intensity

%% Loop through all catchments
for s =1:size(a1,1)
    % Extract catchment-specific average values and regression coefficients
    Q = Qobs_aveg(s, 1);  % Average streamflow
    P = P_aveg(s, 1);  % Average precipitation
    E0 = E0_aveg(s, 1);  % Average potential evapotranspiration
    n = n_aveg(s, 1);  % Average Budyko parameter
    a = a1_all(s, 1);  % Coefficient for ¦¤fs
    b = a2_all(s, 1);  % Coefficient for ¦¤LAI
    c = a3_all(s, 1);  % Coefficient for ¦¤Pi
    e = a0_all(s, 1);  % Constant term
    
    % Sensitivity Calculations (based on Budyko model and LASSO regression)
    sensitivity_var(s,1)= (E0*P*P^(n - 1))/(E0^n + P^n)^(1/n + 1) - E0/(E0^n + P^n)^(1/n) + 1; % Q to P
    sensitivity_var(s,2)= (E0*E0^(n - 1)*P)/(E0^n + P^n)^(1/n + 1) - P/(E0^n + P^n)^(1/n); % Q to E0
    sensitivity_var(s,3)= (-E0*P*(log(E0^n + P^n)/(n^2*(E0^n + P^n)^(1/n)) - (E0^n*log(E0) + P^n*log(P))/(n*(E0^n + P^n)^(1/n + 1)))).*b; % Q to fs
    sensitivity_var(s,4)= (-E0*P*(log(E0^n + P^n)/(n^2*(E0^n + P^n)^(1/n)) - (E0^n*log(E0) + P^n*log(P))/(n*(E0^n + P^n)^(1/n + 1)))).*a; % Q to LAI
    sensitivity_var(s,5)= (-E0*P*(log(E0^n + P^n)/(n^2*(E0^n + P^n)^(1/n)) - (E0^n*log(E0) + P^n*log(P))/(n*(E0^n + P^n)^(1/n + 1)))).*c; % Q to Pi
    deltaQ_deltan(s,1) = (-E0*P*(log(E0^n + P^n)/(n^2*(E0^n + P^n)^(1/n)) - (E0^n*log(E0) + P^n*log(P))/(n*(E0^n + P^n)^(1/n + 1))));% Q to n
    
    % Elasticity Calculations
    elastivity_var(s,1) = P_aveg(s,1).*sensitivity_var(s,1)./Qobs_aveg(s,1);% Q to P
    elastivity_var(s,2) = E0_aveg(s,1).*sensitivity_var(s,2)./Qobs_aveg(s,1);% Q to E0
    elastivity_var(s,3) = fs_aveg(s,1).*sensitivity_var(s,3)./Qobs_aveg(s,1);% Q to fs
    elastivity_var(s,4) = LAI_aveg(s,1).*sensitivity_var(s,4)./Qobs_aveg(s,1);% Q to LAI
    elastivity_var(s,5) = Pi_aveg(s,1).*sensitivity_var(s,5)./Qobs_aveg(s,1);% Q to Pi
end
