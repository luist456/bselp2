%% Problem set 2: short-run and long-run reestrictions in LPs
clear all
close all
clc
addpath("routines");


%% Question 1: Short run restrictions

load("Data_RR_MP.mat"); % Load data

RealRate = dataYyIRFsUSED(:,1); % Real Interest Rate
gdp = dataYyIRFsUSED(:,2); % Real Gross Domestic Product 
DurExp = dataYyIRFsUSED(:,4); % Real Durable Expenditure
NonDurExp = dataYyIRFsUSED(:,5); % Real Non-Durable Expenditure
rr_shock = MPproxyUSED; % Romer-Romer shock series
ff_rate = FedFundrate; % Federal Funds Rate


%% Question 1: Cholesky comparison

y = [gdp DurExp NonDurExp ff_rate RealRate]; % Variables for VAR 
% Order is important! --> Slow - ff_rate - Fast

[T,N] = size(y);

p = 4; % Maximum number of lags in reduced form VAR estimation 
c = 1; % Add constant
H = 20; % Maximum horizon of IRFs
nboot = 1000; % Number of bootstrap repetitions 
prc = 68; % Confidence bands
cumulate = []; % Do not cumulate IRFs

[beta, residuals] = VAR(y,p,c); % Estimate reduced form VAR parameters and
% residuals

wold = woldirf(beta, c, p, H); % Estimate Wold IRFs (not structural)

sigma = (residuals' * residuals) ./ (T - 1 -p - N*p); % Compute variance 
% matrix of reduced form residuals

S = chol(sigma, "lower"); % Compute Cholesky lower triangular factor 

cholirf  = choleskyIRF(wold, S); % Compute Cholesky IRFs

[bootchol, upper, lower, boot_beta] = ...
    bootstrapChol(y, p, c, beta, residuals, nboot, H, prc, cumulate); % Compute
% bootrstrap confidence bands

w = 1 / cholirf(4,4,1); % Normalizing constant --> To read: 4th variable 
% (i.e. FedFundRate), 4th shock (i.e. MP shock), 1st horizon (i.e. impact)


% --- Plot ---
varNames = {'GDP','DurExp','NonDurExp','FedFundRate', 'RealRate'};

horizon = 0:H;

figure;
for v = 1:length(varNames)

    subplot(2,3,v)
    hold on

    irf = squeeze(cholirf(v,4,:)) * w;
    ub  = squeeze(upper(v,4,:))  * w;
    lb  = squeeze(lower(v,4,:))  * w;

    % Confidence band (shaded)
    fill([horizon fliplr(horizon)], ...
         [ub' fliplr(lb')], ...
         [0.75 0.75 0.75], ...
         'EdgeColor','none', ...
         'FaceAlpha',0.4);

    % IRF line
    plot(horizon, irf, 'k', 'LineWidth', 2)

    % zero line
    yline(0,'--k')

    title(varNames{v})
    xlabel('Horizon')
    ylabel('Response')
    grid on
    box on

    hold off

end

