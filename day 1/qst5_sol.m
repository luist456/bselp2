%% Question 5: Comparison

% --- GDP ---
y = gdp;
X_rr = lagmatrix([rr_shock gdp], 1:4);
X_sw = lagmatrix([sw_shock gdp], 1:4);
hStart = 0; % Start from 0
lpMode = "cum"; % Long diff specification
seType = "hac"; % Newey West standard error
c = 1; % Add constant

res_rr = lp_ols(y, X_rr, rr_shock, H, hStart, lpMode, seType, c); % Run the function
res_sw = lp_ols(y, X_sw, sw_shock, H, hStart, lpMode, seType, c); % Run the function


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
title('GDP', 'FontSize', 14)

legend([p1 p2], {'Romer & Romer', 'Smets & Wouters'}, 'Location', 'Best')

grid on
box on
hold off

