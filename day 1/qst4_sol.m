%% Durable Expenditure
y = DurExp; % Dependent variable
d_DurExp = [NaN(1,1); diff(DurExp)];
X = [lagmatrix([rr_shock d_DurExp], 1:8)] ; % Controls
S =  rr_shock; % Shock
H = 20; % Maximum horizon of IRF
hStart = 0; % Start from 0
lpMode = "cum"; % Long diff specification
seType = "ols"; % Newey West standard error
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
title('Durable Expenditure -- OLS', 'FontSize', 14)

legend('±1 SE Band', 'Estimate', 'Location', 'Best')
grid on
box on

hold off


%% Non Durable Expenditure

y = NonDurExp; % Dependent variable
d_NonDurExp = [NaN(1,1); diff(NonDurExp)];
X = [lagmatrix([rr_shock d_NonDurExp], 1:8)]; % Controls
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
title('Non-Durable Expenditure -- OLS', 'FontSize', 14)

legend('±1 SE Band', 'Estimate', 'Location', 'Best')
grid on
box on

hold off