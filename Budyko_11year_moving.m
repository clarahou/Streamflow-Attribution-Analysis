%*************************************************************************%
% Description: This is the code for calculating parameter n  Choudhury¡¯s
% formulation of the Budyko model based on 11-year sliding window for each
% catchment in the analysis
% Author: Ying Hou
% Date: 24th Jan, 2025
%*************************************************************************%

%% Load data
% Load yearly precipitation series data for all catchments from 1982 to 2016
load ('Pre_yearly_1982_2016','Pre_yearly');  
% Load yearly streamflow series data for all catchments from 1982 to 2016
load ('Q_yearly_1982_2016','Q_yearly');      
% Load yearly evapotranspiration series data for all catchments from 1982 to 2016
load ('E0_yearly_1982_2016','E0_yearly');    

Catchment_Number = 931; % Number of catchments

%% Parameter n calculation for each sliding window
for moving = 1:25 % There are 25 sliding windows from 1982 to 2016
    % Select the relevant data for the sliding window (11 years)
    Q_ob_all = Q_yearly(:, moving:moving+10);   % Streamflow for the 11-year window
    P_all = Pre_yearly(:, moving:moving+10);    % Precipitation for the 11-year window
    E0_all = E0_yearly(:, moving:moving+10);    % Evapotranspiration for the 11-year window
    
    % Initialize parameter n and Budyko model variables
    n_initial_3 = ones(Catchment_Number, 1) * 0.05;
    n_initial_4 = n_initial_3 + 0.01;
    E_Budyko_Initial_3 = ones(Catchment_Number, 1) * -9999;
    E_Budyko_Initial_4 = ones(Catchment_Number, 1) * -9999;
    
    % Calculate the modeled evapotranspiration (E) using water balance
    for i = 1:Catchment_Number
        E_model_all(i, 1) = nanmean(P_all(i, :)) - nanmean(Q_ob_all(i, :));  % Calculate E based on water balance
    end

    % Calculate the initial evapotranspiration (E) using Budyko's model
    for i = 1:Catchment_Number
        P_all(i) = nanmean(P_all(i, :)); 
        E0_all(i) = nanmean(E0_all(i, :));
        
        % Calculate the initial E for parameter n = 0.05 and n = 0.06 (n_initial_4)
        E_Budyko_Initial_3(i) = P_all(i) * E0_all(i) / (P_all(i)^n_initial_3(i) + E0_all(i)^n_initial_3(i))^(1/n_initial_3(i));
        E_Budyko_Initial_4(i) = P_all(i) * E0_all(i) / (P_all(i)^n_initial_4(i) + E0_all(i)^n_initial_4(i))^(1/n_initial_4(i));
    end

    % Calculate the differences between the model results and water balance results
    Diff_3 = abs(E_Budyko_Initial_3 - E_model_all);
    Diff_4 = abs(E_Budyko_Initial_4 - E_model_all);

    % Iterative approach for refining the value of n
    for i = 1:Catchment_Number
        while Diff_4(i) < Diff_3(i)
            n_initial_3(i) = n_initial_4(i);
            n_initial_4(i) = n_initial_4(i) + 0.01;  % Increment n by 0.01 for the next iteration
            
            P_all(i) = nanmean(P_all(i, :)); 
            E0_all(i) = nanmean(E0_all(i, :));
            
            % Recalculate E based on Budyko's model with updated n
            E_Budyko_Initial_3(i) = P_all(i) * E0_all(i) / (P_all(i)^n_initial_3(i) + E0_all(i)^n_initial_3(i))^(1/n_initial_3(i));
            E_Budyko_Initial_4(i) = P_all(i) * E0_all(i) / (P_all(i)^n_initial_4(i) + E0_all(i)^n_initial_4(i))^(1/n_initial_4(i));
            
            % Update the differences
            Diff_3(i) = abs(E_model_all(i, 1) - E_Budyko_Initial_3(i));
            Diff_4(i) = abs(E_model_all(i, 1) - E_Budyko_Initial_4(i));
        end
        % Store the value of n for each catchment and window
        n_11year(i, moving) = n_initial_3(i);
    end
end