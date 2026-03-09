%% Step 1:

data = readtable("qst2_data.xlsx"); % Read data


%% Step 2:

prod = 100 .* log(data.OUTNFB ./ data.HOANBS); % (log) labor productivity

pch = 100 .* log(data.HOANBS ./ data.CNP16OV); % (log) per-capita hours


%% Step 3:

p = 4; % Maximum number of lags of controls

d_prod = [NaN(1,1); diff(prod)];
d_pch = [NaN(1,1); diff(pch)];
X = lagmatrix([d_prod d_pch], 1:p); % Matrix of controls for both steps


%% Step 4:

H = 20;
lpMode = "cum";
seType = "hac";
c = 1;
hStart = 0;
y = prod;
S = [d_prod d_pch]; % I want the coefficients of these two variables

res = lp_ols(y, X, S, H, hStart, lpMode, seType, c);

w_1 = res.beta(H+1,1);
w_2 = res.beta(H+1,2);

lin_comb = w_1 .* d_prod + w_2 .* d_pch; % Resulting linear combination that...


%% Step 5:

H = 20;
lpMode = "cum";
seType = "hac";
c = 1;
hStart = 0;
S = lin_comb; % Signal of technology shock 
LHS = [prod pch];
varNames = {'LaborProd','Hours'};


for v = 1:size(LHS,2) % Loop over variables

    res = lp_ols(LHS(:,v), X, S, H, hStart, lpMode, seType, c);

    results.beta.(varNames{v}) = res.beta;
    results.se.(varNames{v}) = res.se;

end


%% Step 6:

horizon = 0:H;

figure;
for v = 1:length(varNames)

    subplot(2,1,v)
    hold on
    
    % Point estimate
    beta = results.beta.(varNames{v});
    se   = results.se.(varNames{v});
    
    lo = beta - se;     % 68% band
    hi = beta + se;

    % --- Shaded confidence band ---
    patch([horizon fliplr(horizon)], ...
          [lo' fliplr(hi')], ...
          'k', ...
          'FaceAlpha', 0.12, ...
          'EdgeColor', 'none', ...
          'HandleVisibility','off');

    % --- Dashed band edges ---
    plot(horizon, lo, 'k--', 'LineWidth', 1, 'HandleVisibility','off');
    plot(horizon, hi, 'k--', 'LineWidth', 1, 'HandleVisibility','off');

    % --- IRF line ---
    plot(horizon, beta, 'k-', 'LineWidth', 2);

    yline(0,'--','HandleVisibility','off')

    title(varNames{v})
    xlabel('Horizon')
    ylabel('Response')
    grid on
    hold off

end


%% Step 7:

H = 20;
lpMode = "cum";
seType = "hac";
c = 1;
hStart = 0;
S = d_pch; % Orthogonal to technology shock
X = [lin_comb lagmatrix([d_prod d_pch], 1:p)]; % New controls --> added 
% linear combination
LHS = [prod pch];
varNames = {'LaborProd','Hours'};


for v = 1:size(LHS,2) % Loop over variables

    res = lp_ols(LHS(:,v), X, S, H, hStart, lpMode, seType, c);

    results.beta.(varNames{v}) = res.beta;
    results.se.(varNames{v}) = res.se;

end


% --- Plot ---

horizon = 0:H;

figure;
for v = 1:length(varNames)

    subplot(2,1,v)
    hold on
    
    % Point estimate
    beta = results.beta.(varNames{v});
    se   = results.se.(varNames{v});
    
    lo = beta - se;     % 68% band
    hi = beta + se;

    % --- Shaded confidence band ---
    patch([horizon fliplr(horizon)], ...
          [lo' fliplr(hi')], ...
          'k', ...
          'FaceAlpha', 0.12, ...
          'EdgeColor', 'none', ...
          'HandleVisibility','off');

    % --- Dashed band edges ---
    plot(horizon, lo, 'k--', 'LineWidth', 1, 'HandleVisibility','off');
    plot(horizon, hi, 'k--', 'LineWidth', 1, 'HandleVisibility','off');

    % --- IRF line ---
    plot(horizon, beta, 'k-', 'LineWidth', 2);

    yline(0,'--','HandleVisibility','off')

    title(varNames{v})
    xlabel('Horizon')
    ylabel('Response')
    grid on
    hold off

end
