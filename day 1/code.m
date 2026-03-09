%% Problem set 1: Monetary Policy Shocks and Transmission Mechanisms
clear all
close all
clc


%% Question 1: Data Visualization and Preliminary Analysis

load("Data_RR_MP.mat"); % Load data

% Extract variables
dates = DatesUSED;                     % Quarterly dates (1969Q1–2007Q3)
ff_rate = FedFundrate;                 % Federal Funds Rate
rr_shock = MPproxyUSED;                % Romer–Romer Monetary Policy Shock
RealRate = dataYyIRFsUSED(:,1);        % Real Interest Rate
gdp = dataYyIRFsUSED(:,2);             % Real Gross Domestic Product
TotExp = dataYyIRFsUSED(:,3);          % Real Total Expenditure
DurExp = dataYyIRFsUSED(:,4);          % Durable Expenditure
NonDurExp = dataYyIRFsUSED(:,5);       % Non-Durable Expenditure


% Create figure
figure('Color', 'w', 'Position', [100, 100, 900, 400]);

% Plot both series on the same axis
plot(dates, ff_rate, 'LineWidth', 1.6, 'Color', [0 0.45 0.74]); hold on;
plot(dates, rr_shock, 'LineWidth', 1.4, 'LineStyle', '--', 'Color', [0.85 0.33 0.10]);

% Formatting
title('Federal Funds Rate and Romer–Romer Monetary Policy Shock');
xlabel('Time');
ylabel('Value');
grid on;
legend({'Federal Funds Rate', 'Romer–Romer MP Shock'}, 'Location', 'best');
set(gca, 'FontSize', 12);


%% Same figure but with Federal Funds Rate in first differences

% Compute first differences of the Federal Funds Rate
d_ff_rate = diff(ff_rate);

% Align dates and shock (drop first observation)
dates_diff = dates(2:end);
rr_shock_diff = rr_shock(2:end);

% Create figure
figure('Color', 'w', 'Position', [100, 100, 900, 400]);

% Plot both series on the same axis
plot(dates_diff, d_ff_rate, 'LineWidth', 1.6, 'Color', [0 0.45 0.74]); hold on;
plot(dates_diff, rr_shock_diff, 'LineWidth', 1.4, 'LineStyle', '--', 'Color', [0.85 0.33 0.10]);

% Formatting
title('First Difference of Federal Funds Rate and Romer–Romer MP Shock');
xlabel('Time');
ylabel('Quarterly Change');
grid on;
legend({'Δ Federal Funds Rate', 'Romer–Romer MP Shock'}, 'Location', 'best');
set(gca, 'FontSize', 12);

%% Question 2: A toy example

H = 10; % Number of horizons

betas_cum_gdp = nan(H+1,1); % To store shock coefficients
ses   = nan(H+1,1); % To store corresponding standard deviations

T = length(gdp); % Number of observations

for h = 0:H % Loop over horizons
 
    t0 = 2;
    t1 = T - h; 

    Y = gdp((t0+h):(t1+h)) - gdp((t0-1):(t1-1)); % Long differences of 
    % dependent variable. y_{t+h} - y_{t-1}

    X = [ones(length(Y),1), rr_shock(t0:t1)]; % Matrix of regressors: constant + 
    % shock at time t (one period-ahead of y_{t-1})

    % OLS
    b = (X' * X) \ (X' * Y); % Coefficients (X'X)^{-1}X'Y
    res = Y - X*b; % Residuals

    n = size(X,1);
    k = size(X,2);
    s2 = (res' * res) / (n - k);
    V  = s2 * inv(X' * X);

    betas_cum_gdp(h+1) = b(2); % Pick second coefficient (first is the constant)
    ses(h+1)   = sqrt(V(2,2)); % Same here
end

hgrid = (0:H)';

% 68% bands: +/- 1 * se
upper68_cum_gdp = betas_cum_gdp + ses;
lower68_cum_gdp = betas_cum_gdp - ses;

% Final plot
figure;
plot(hgrid, betas_cum_gdp, '-o'); hold on; grid on; % IRF
ci_color = [0.5 0.5 0.5];   % grey (RGB)
plot(hgrid, upper68_cum_gdp, '--', 'Color', ci_color, 'LineWidth', 1.2); % Upper bound CB
plot(hgrid, lower68_cum_gdp, '--', 'Color', ci_color, 'LineWidth', 1.2); % Lower bound CB
yline(0, '-');
xlabel('Horizon h');
ylabel('\beta_h');
title('Local Projection: GDP_{t+h} - GDP_{t-1} on shock_t (const + shock)');


%% Question 4: Estimation and Impulse Response Functions (IRFs)

%% --- Federal Fund Rate ---
y = ff_rate;
d_ff_rate = [NaN(1,1); diff(ff_rate)];
X = lagmatrix([rr_shock d_ff_rate], 1:4); % Controls
S = rr_shock; % Shock
H = 20; % Maximum horizon of IRF
hStart = 0; % Start from 0
lpMode = "cum"; % Long diff specification
seType = "hac"; % Newey West standard error
c = 1; % Add constant

res = lp_ols(y, X, S, H, hStart, lpMode, seType, c); % Run the function

w = (1 / res.beta(1)) / 2; % Impact 50bp normalization

h     = res.h;
beta  = res.beta .* w; % Normalized response
se    = res.se .* w; % Normalized standard errors

upper = beta + se; % 68% confidence bands
lower = beta - se;

figure
hold on

% Shaded confidence band
fill([h; flipud(h)], ...
     [upper; flipud(lower)], ...
     [0.8 0.8 1], ...
     'EdgeColor','none', ...
     'FaceAlpha',0.4);

% Coefficient line
plot(h, beta, 'b-', 'LineWidth', 2);

% Zero line
yline(0, 'k--', 'LineWidth', 1);

xlabel('Horizon', 'FontSize', 12)
ylabel('Coefficient', 'FontSize', 12)
title('FedFundRate -- Newey West', 'FontSize', 14)

legend('±1 SE Band', 'Estimate', 'Location', 'Best')
grid on
box on

hold off


%% --- GDP --- 

y = gdp; % Dependent variable
d_gdp = [NaN(1,1); diff(gdp)];
X = lagmatrix([rr_shock d_gdp], 1:4); % Controls
S = rr_shock; % Shock
H = 20; % Maximum horizon of IRF
hStart = 0; % Start from 0
lpMode = "cum"; % Long diff specification
seType = "hac"; % Newey West standard error
c = 1; % Add constant

res = lp_ols(y, X, S, H, hStart, lpMode, seType, c); % Run the function

h     = res.h;
beta  = res.beta .* w; % Normalized response
se    = res.se .* w; % Normalized standard errors

upper = beta + se; % 68% confidence bands
lower = beta - se;

figure
hold on

% Shaded confidence band
fill([h; flipud(h)], ...
     [upper; flipud(lower)], ...
     [0.8 0.8 1], ...
     'EdgeColor','none', ...
     'FaceAlpha',0.4);

% Coefficient line
plot(h, beta, 'b-', 'LineWidth', 2);

% Zero line
yline(0, 'k--', 'LineWidth', 1);

xlabel('Horizon', 'FontSize', 12)
ylabel('Coefficient', 'FontSize', 12)
title('GDP -- Newey West', 'FontSize', 14)

legend('±1 SE Band', 'Estimate', 'Location', 'Best')
grid on
box on

hold off


%% Question 5: Comparison
load("sw_shock.mat"); % Load Smets and Wouters shock

y = ff_rate;
X_rr = lagmatrix([rr_shock d_ff_rate], 1:4); % Control matrix for RR
X_sw = lagmatrix([sw_shock d_ff_rate], 1:4); % Control matrix for SW
hStart = 0; % Start from 0
lpMode = "cum"; % Long diff specification
seType = "hac"; % Newey West standard error
c = 1; % Add constant

res_rr = lp_ols(y, X_rr, rr_shock, H, hStart, lpMode, seType, c); % Run for RR
res_sw = lp_ols(y, X_sw, sw_shock, H, hStart, lpMode, seType, c); % Run for SW

w_rr = 1 ./ sum(res_rr.beta(1:6)); % Normalization for RR
w_sw = 1 ./ sum(res_sw.beta(1:6)); % Normalization for SW


% --- Plot ---

% IRF 1 (rr)
beta_rr  = res_rr.beta .* w_rr;
se_rr    = res_rr.se   .* w_rr;
upper_rr = beta_rr + se_rr;
lower_rr = beta_rr - se_rr;

% IRF 2 (sw)
beta_sw  = res_sw.beta .* w_sw;
se_sw    = res_sw.se   .* w_sw;
upper_sw = beta_sw + se_sw;
lower_sw = beta_sw - se_sw;

figure
hold on

% --- Shaded 68% CI: rr ---
fill([h; flipud(h)], ...
     [upper_rr; flipud(lower_rr)], ...
     [0.8 0.8 1], ...
     'EdgeColor','none', ...
     'FaceAlpha',0.35);

% --- Shaded 68% CI: sw ---
fill([h; flipud(h)], ...
     [upper_sw; flipud(lower_sw)], ...
     [1 0.85 0.85], ...
     'EdgeColor','none', ...
     'FaceAlpha',0.35);

% --- IRF lines (store handles!) ---
p1 = plot(h, beta_rr, 'b-', 'LineWidth', 2);
p2 = plot(h, beta_sw, 'r-', 'LineWidth', 2);

% Zero line
yline(0, 'k--', 'LineWidth', 1);

xlabel('Horizon', 'FontSize', 12)
ylabel('Coefficient', 'FontSize', 12)
title('Federal Funds Rate', 'FontSize', 14)

legend([p1 p2], {'Romer & Romer', 'Smets & Wouters'}, 'Location', 'Best')

grid on
box on
hold off



