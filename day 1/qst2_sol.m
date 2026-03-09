%% Question 2: A toy example

H = 10; % Number of horizons

betas = nan(H+1,1); % To store shock coefficients
ses   = nan(H+1,1); % To store corresponding standard deviations

T = length(RealRate); % Number of observations

for h = 0:H % Loop over horizons
 
    t0 = 2;
    t1 = T - h;

    Y = RealRate((t0+h):(t1+h)) - RealRate((t0-1):(t1-1)); % Long differences of 
    % dependent variable

    X = [ones(length(Y),1), rr_shock(t0:t1)]; % Matrix of regressors: constant + 
    % shock

    % OLS
    b = (X' * X) \ (X' * Y); % Coefficients 
    res = Y - X*b; % Residuals

    n = size(X,1);
    k = size(X,2);
    s2 = (res' * res) / (n - k);
    V  = s2 * inv(X' * X);

    betas(h+1) = b(2); % Pick second coefficient (first is the constant)
    ses(h+1)   = sqrt(V(2,2)); % Same here
end

hgrid = (0:H)';

% 68% bands: +/- 1 * se
upper68 = betas + ses;
lower68 = betas - ses;

% Final plot
figure;
plot(hgrid, betas, '-o'); hold on; grid on;
ci_color = [0.5 0.5 0.5];   % grey (RGB)
plot(hgrid, upper68, '--', 'Color', ci_color, 'LineWidth', 1.2);
plot(hgrid, lower68, '--', 'Color', ci_color, 'LineWidth', 1.2);
yline(0, '-');
xlabel('Horizon h');
ylabel('\beta_h');
title('Local Projection: rr_{t+h} - rr_{t-1} on shock_t (const + shock)');


