%% Question 3: State-dependence 

%% Unemployment

% General inputs for LP 
p_y = 2;
p_x = 2;
hStart = 0;
c = 0;
H = 25;
lpMode = "cum";
S = [shock_data_linear(good) interaction]; % Linear shock & interaction

[YY, X] = make_regressors_func(Z,"URATE",p_y,p_x); % Extract regressor matrix

ebp_avg = mean(EBPhat); % Average of fitted values

results = cell(length(p_vals),1); % To store the results

for i = 1:length(p_vals)
    
    p = p_vals(i);
    R = [1, p - ebp_avg]; % Condition on a particular deviation from the mean
    
    results{i} = lp_ols_interaction(YY(good), X(good,:), S, H, hStart, lpMode, c, R);
    
end


% --- Plot ---

figure;
hold on;

% Choose ONE color (black here, but change if you want)
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
title("Unemployment Rate");
grid on;

hold off;


%% Industrial production 

% General inputs for LP 
p_y = 2;
p_x = 2;
hStart = 0;
c = 0;
H = 25;
S = [shock_data_linear(good) interaction]; % Linear shock & interaction

[YY, X] = make_regressors_func(Z,"IP",p_y,p_x); % Extract regressor matrix

ebp_avg = mean(EBPhat); % Average of fitted values

for i = 1:length(p_vals)
    
    p = p_vals(i);
    R = [1, p - ebp_avg]; % Condition on a particular deviation from the mean
    
    results{i} = lp_ols_interaction(YY(good), X(good,:), S, H, hStart, "cum", c, R);
    
end


% --- Plot ---

figure;
hold on;

% Choose ONE color (black here, but change if you want)
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
title("Industrial Production");
grid on;

hold off;