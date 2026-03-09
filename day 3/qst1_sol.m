%% Question 1.1: LP-IV

%% --- GDP ---
y = gdp; % Dependent variable
X = lagmatrix([d_gdp rr_shock], 1:p); % Control matrix
S = d_ff_rate; % To be instrumented (Federal Funds Rate)
Ziv = rr_shock; % Instrument (Romer-Romer)
H = 20; % Maximum horizon of IRF
hStart = 0; % Start LP at h = 0
lpMode = "cum"; % Long-difference
c = 1; % Add constant

res = lp_iv(y, X, S, Ziv, H, hStart, lpMode, c);

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
title('GDP -- LP-IV', 'FontSize', 14)

legend('±1 SE Band', 'Estimate', 'Location', 'Best')
grid on
box on

hold off


%% --- Real Interest Rate ---

y = RealRate; % Dependent variable
X = lagmatrix([d_rr_rate rr_shock], 1:p); % Control matrix
S = d_ff_rate; % To be instrumented
Ziv = rr_shock; % Instrument
H = 20; % Maximum horizon of IRF
hStart = 0; % Start LP at h = 0
lpMode = "cum"; % Long-difference
c = 1; % Add constant

res = lp_iv(y, X, S, Ziv, H, hStart, lpMode, c);

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
title('Real Interest Rate -- LP-IV', 'FontSize', 14)

legend('±1 SE Band', 'Estimate', 'Location', 'Best')
grid on
box on

hold off


%% Question 1.2: Multiple instruments

% --- IRF of Real Interest Rate ---
y = RealRate; % Dependent variable

X_rr = lagmatrix([d_rr_rate rr_shock], 1:p); % Control matrix
X_sw = lagmatrix([d_rr_rate sw_shock], 1:p); % Control matrix
X_b = lagmatrix([d_rr_rate rr_shock sw_shock], 1:p); % Control matrix

H = 20; % Maximum horizon of IRF
hStart = 0; % Start LP at h = 0
lpMode = "cum"; % Long-difference
c = 1; % Add constant

res_rr = lp_iv(y, X_rr, S, Ziv_rr, H, hStart, lpMode, c);
res_sw = lp_iv(y, X_sw, S, Ziv_sw, H, hStart, lpMode, c);
res_b = lp_iv(y, X_b, S, Ziv_b, H, hStart, lpMode, c);


% --- Plot ---

% --- Horizon (shared) ---
h = res_rr.h(:);

% --- Extract IRFs and 68% bands (±1 SE) ---
b_rr  = w_rr * res_rr.beta(:);  se_rr = w_rr * res_rr.se(:);
lo_rr = b_rr - se_rr;    hi_rr = b_rr + se_rr;

b_sw  = w_sw * res_sw.beta(:);  se_sw = w_sw * res_sw.se(:);
lo_sw = b_sw - se_sw;    hi_sw = b_sw + se_sw;

b_b   = w_b * res_b.beta(:);   se_b  = w_b * res_b.se(:);
lo_b  = b_b - se_b;      hi_b  = b_b + se_b;

% --- Plot ---
figure; hold on; box on;

% Confidence bands (shaded)
plot(h, lo_rr, 'k--', 'LineWidth', 1, 'HandleVisibility','off');
plot(h, hi_rr, 'k--', 'LineWidth', 1, 'HandleVisibility','off');

patch([h; flipud(h)], [lo_sw; flipud(hi_sw)], 'b', ...
    'FaceAlpha', 0.12, 'EdgeColor', 'none', ...
    'HandleVisibility','off');

plot(h, lo_b, 'r--', 'LineWidth', 1, 'HandleVisibility','off');
plot(h, hi_b, 'r--', 'LineWidth', 1, 'HandleVisibility','off');

% IRF lines (on top of bands)
l1 = plot(h, b_rr, 'k-', 'LineWidth', 2);
l2 = plot(h, b_sw, 'b-', 'LineWidth', 2);
l3 = plot(h, b_b,  'r-', 'LineWidth', 2);

% Zero line
yline(0, '--', 'LineWidth', 1);

% Cosmetics
xlim([min(h) max(h)]);
xlabel('Horizon');
ylabel('Response');
title('Real Interest Rate');

legend([l1 l2 l3], {'Romer-Romer', 'Smets-Wouters', 'Both'}, 'Location', 'best');
set(gca, 'FontSize', 12);
hold off;
