%% Question 2.1: Replication Using Standard Local Projections
close all
clear all
clc

%%
% Load necessary data
load("shock_data_linear.mat"); % Contains shock series
load("Z.mat"); % Contains data for all variables
load("LHSlabels.mat"); % Contains the list of dependent variables


% Build a Tx2 array with in the first column the linear shock and in the 
% second column only the positive values
shock_pos = shock_data_linear;
shock_pos(shock_data_linear < 0) = 0; % First column (postive values)

shock_array = [shock_data_linear, shock_pos]; % Final array


% Estimate LP --- Employment ---
p_y = 12; % Maximum number of lags for depedent variable
p_x = 12; % Maximum number of lags for other controls
H =  25; % Maximum number of horizons
hStart = 0;
lpMode = "cum";
c = 0; % Constant already added by make_regressors_func
 
[YY, X] = make_regressors_func(Z,"EMP",p_y,p_x); % Extract regressor matrix
R_pos = [1,1]; % Linear combination to have positive shock
R_neg = [1,0]; % Linear combination to have negative shock
res_pos = lp_ols_interaction(YY, X, shock_array, H, hStart, lpMode, c, R_pos);
res_neg = lp_ols_interaction(YY, X, shock_array, H, hStart, lpMode, c, R_neg);

w = 0.05;

beta_pos = res_pos.irf .* w;
se_pos = res_pos.se_irf .* w;
h = res_pos.h;

beta_neg = -res_neg.irf .* w; % Flipped sign for comparison
se_neg = res_neg.se_irf .* w;

z = 1;
% 68% confidence bands
lo_pos = beta_pos - z*se_pos;
hi_pos = beta_pos + z*se_pos;

lo_neg = beta_neg - z*se_neg;
hi_neg = beta_neg + z*se_neg;


figure('Color','w'); hold on; box on;

c_neg = [0 0.447 0.741];      % blue  -> negative shock
c_pos = [0.850 0.325 0.098];  % red   -> positive shock

% CI areas 
fill([h; flipud(h)], [lo_pos; flipud(hi_pos)], c_pos, ...
    'FaceAlpha',0.2,'EdgeColor','none','HandleVisibility','off');
fill([h; flipud(h)], [lo_neg; flipud(hi_neg)], c_neg, ...
    'FaceAlpha',0.2,'EdgeColor','none','HandleVisibility','off');

% IRFs 
p1 = plot(h, beta_pos, 'LineWidth',2, 'Color',c_pos);
p2 = plot(h, beta_neg, 'LineWidth',2, 'Color',c_neg);

yline(0,'k','LineWidth',1,'HandleVisibility','off');

title("EMP",'Interpreter','none');
xlabel('Horizon');
ylabel('Response');

legend([p1 p2], {'IRF (+)','IRF (-)'}, 'Location','best');

grid on;


%% Question 2.2: Replication Using Standard Local Projections


% Estimate LP --- Employment ---
p_y = 12; % Maximum number of lags for depedent variable
p_x = 12; % Maxiimum number of lags for other controls
H =  25; % Maximum number of horizons
hStart = 0;
lpMode = "cum";
c = 0; % Constant already added by make_regressors_func
 
[YY, X] = make_regressors_func(Z,"EMP",p_y,p_x); % Extract regressor matrix
R_diff = [0,1]; % Linear combination for difference
res_diff = lp_ols_interaction(YY, X, shock_array, H, hStart, "cum", c, R_diff);


beta_diff = res_diff.irf .* w;
se_diff = res_diff.se_irf .* w;
h = res_diff.h;

z = 1;
% 68% confidence bands
lo_diff = beta_diff - z*se_diff;
hi_diff = beta_diff + z*se_diff;


figure('Color','w'); hold on; box on;

c_diff = [0.000 0.500 0.000];  % dark green

% CI dashed borders only
plot(h, lo_diff, '--', 'Color', c_diff, 'LineWidth', 1.2, ...
    'HandleVisibility','off');
hold on
plot(h, hi_diff, '--', 'Color', c_diff, 'LineWidth', 1.2, ...
    'HandleVisibility','off');

% IRFs 
p1 = plot(h, beta_diff, 'LineWidth',2, 'Color',c_diff);


yline(0,'k','LineWidth',1,'HandleVisibility','off');

title("EMP",'Interpreter','none');
xlabel('Horizon');
ylabel('Response');

legend(p1, {'Difference'}, 'Location','best');

grid on;
