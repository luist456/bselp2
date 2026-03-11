%% Problem set 4: More nonlinearities and panel LP
clear all
close all
clc

addpath("mat_files");
addpath("routines");
%% Question 1: What other sources of nonlinearity might be relevant in this framework?

% Load necessary data
load("shock_data_linear.mat"); % Contains shock series
load("Z.mat"); % Contains data for all variables
load("LHSlabels.mat"); % Contains the list of dependent variables


shock_sq = [shock_data_linear, shock_data_linear.^2]; % Build a two columns 
% array containing the shocks: linear shock in the first column, squared
% shock in the second column

% Estimate LP with squared shocks

% --- Common inputs ---
p_y = 12; % Maximum number of lags for depedent variable
p_x = 12; % Maximum number of lags for other controls
H =  25; % Maximum number of horizons
hStart = 0;
seType = "ols";
lpMode = "cum";
c = 0; % Constant already added by make_regressions_func


for v = 1:length(LHSlabels) % Loop over all variables

    Y_LABEL = LHSlabels(v);   
    Yname   = Y_LABEL{:};     
    disp(['Making LP for: ', Yname]) % Pick the name of the current variable

    [YY, X] = make_regressors_func(Z,Y_LABEL,p_y,p_x); % Extract regressor matrix

    res = lp_ols(YY, X, shock_sq, H, hStart, lpMode, seType, c);

    LPcoeffsAsy_LHS.("beta").(Yname) = res.beta; % Store shocks coefficients

    LPcoeffsAsy_LHS.("se").(Yname)  = res.se; % Store corresponding standard errors
end

disp('--------- Done estimating asymmetric LPs ---------')


%% Question 3: State-dependence

load("unc_shock.mat"); % Load uncertainty shock 

ebp = Z.EBP; % Save 'state' variable

% --- First stage regression ---
X = [ones(length(ebp),1) lagmatrix( [Z.EBP Z.URATE Z.Dl_EMP Z.Dl_IP Z.TBILL1Y], 1:1)]; %
% Set of controls

T = [ebp X unc_KS_tot];
good = all(~isnan(T), 2);
EBP_reg = ebp(good);
X_reg = X(good,:);
unc_KS_tot_reg = unc_KS_tot(good); % Get rid of all rows with NaN values

[bhat, ~, ~, ~, ~] = regress(EBP_reg, [unc_KS_tot_reg X_reg]); % Run regresssion

EBPhat = unc_KS_tot_reg * bhat(1); % Get fitted value, using only uncertainty
% shock

ebp_avg = mean(EBPhat); % Average of fitted values

interaction = shock_data_linear(good) .* (EBPhat - ebp_avg); % Interaction term for LP

perc = 16:6:84;                          
p_vals = prctile(EBPhat(:), perc);  % Create a grid of percentiles of EBPhat
% to use as conditioning state

results = cell(length(p_vals),1); % To store the results

% General inputs for LP 
p_y = 2;
p_x = 2;
hStart = 0;
c = 0;
H = 25;
lpMode = "cum";
S = [shock_data_linear(good) interaction]; % Linear shock & interaction

[YY, X] = make_regressors_func(Z,"EMP",p_y,p_x); % Extract regressor matrix

for i = 1:length(p_vals)
    
    p = p_vals(i);
    R = [1, p - ebp_avg]; % Condition on a particular deviation from the mean
    
    results{i} = lp_ols_interaction(YY(good), X(good,:), S, H, hStart, "cum", c, R);
    
end



% --- Plot ---

figure;
hold on;

baseColor = [0.8 0 0];

% Horizon vector (adjust if needed)
irf0 = results{1}.("irf");
irf0 = irf0(:);   % ensure column
horizons = 0:length(irf0)-1;

% Define line width range
minLW = 0.8;
maxLW = 3.5;

for i = 1:length(results)
    
    irf = results{i}.("irf");
    irf = irf(:);
    
    % Line width increases with percentile
    lw = minLW + (i-1)/(length(results)-1) * (maxLW - minLW);
    
    plot(horizons, irf, ...
         'Color', baseColor, ...
         'LineWidth', lw);
end

% ---- Horizontal zero line ----
yline(0, 'k', 'LineWidth', 1.5);

xlabel('Horizon');
ylabel('IRF');
title("Employment");
grid on;

hold off;

