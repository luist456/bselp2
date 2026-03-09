%% Question 3: A toy example in levels

%% GDP -- growth rates

d_gdp = [NaN(1,1); diff(gdp)];
H = 10; % Number of horizons

betas = nan(H+1,1); % To store shock coefficients
ses   = nan(H+1,1); % To store corresponding standard deviations

T = length(gdp); % Number of observations

for h = 0:H % Loop over horizons
 
    t0 = 2; % To avoid first NaN value given by first difference
    t1 = T - h;

    Y = d_gdp((t0+h):(t1+h));   % Dependent variable in levels y_{t+h}

    X = [ones(length(Y),1), rr_shock(t0:t1)]; % constant + shock

    % OLS
    b = (X' * X) \ (X' * Y); 
    res = Y - X*b;

    n = size(X,1);
    k = size(X,2);
    s2 = (res' * res) / (n - k);
    V  = s2 * inv(X' * X);

    betas(h+1) = b(2);
    ses(h+1)   = sqrt(V(2,2));
end

hgrid = (0:H)';

upper68 = betas + ses;
lower68 = betas - ses;

figure;
plot(hgrid, betas, '-o'); hold on; grid on;
ci_color = [0.5 0.5 0.5];
plot(hgrid, upper68, '--', 'Color', ci_color, 'LineWidth', 1.2);
plot(hgrid, lower68, '--', 'Color', ci_color, 'LineWidth', 1.2);
yline(0, '-');
xlabel('Horizon h');
ylabel('\beta_h');
title('Local Projection: GDP_{t+h} on shock_t (const + shock)');


%% Comparison


figure;
hold on;

% --- Confidence intervals for betas_cum_gdp ---
h_ci = plot(hgrid, upper68_cum_gdp, '--k', 'LineWidth', 1.5);
plot(hgrid, lower68_cum_gdp, '--k', 'LineWidth', 1.5, ...
     'HandleVisibility','off');   % avoid duplicate legend entry

% --- Long difference IRF ---
h_ld = plot(hgrid, betas_cum_gdp, '--k', 'LineWidth', 2);

% --- Cumulative IRF ---
h_cum = plot(hgrid, cumsum(betas), 'r', 'LineWidth', 2);

% Zero line
yline(0, 'k--');

% Legend (ordered cleanly)
legend([h_ld h_ci h_cum], ...
       {'Long difference', ...
        '68% CI (Long difference)', ...
        'Cumulative'}, ...
       'Location', 'Best');

xlabel('Horizon');
ylabel('Response');
title('Impulse Response Comparison');

box on;
hold off;

