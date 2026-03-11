%% Question 2: And now?

%% Option 1: A table
z_10 = 1.645;   % 10% significance
z_32 = 1.0;     % 32% significance 

sig_10_counts = zeros(length(LHSlabels),1);
sig_32_counts = zeros(length(LHSlabels),1);
var_names     = strings(length(LHSlabels),1);

for v = 1:length(LHSlabels)

    Y_LABEL = LHSlabels(v);   
    Yname   = Y_LABEL{:};     
    
    beta = LPcoeffsAsy_LHS.("beta").(Yname)(:,2); % Retrieve coefficients for
    % current variable
    se = LPcoeffsAsy_LHS.("se").(Yname)(:,2);     % Retrieve standard errors for
    % current variable

    % t-statistics
    t_stat = beta ./ se;

    % Count significant horizons
    sig_10_counts(v) = sum(abs(t_stat) > z_10);
    sig_32_counts(v) = sum(abs(t_stat) > z_32);

    var_names(v) = Yname;
end

% Create table
SignificanceTable = table(var_names, sig_10_counts, sig_32_counts, ...
    'VariableNames', {'Variable', 'Significant_10pct', 'Significant_32pct'});

disp('--------- Done estimating asymmetric LPs ---------')
disp(SignificanceTable)

%% Option 2: A plot

% Determine number of horizons
first_var = LHSlabels{1};
H = size(LPcoeffsAsy_LHS.beta.(first_var)(:,2),1);

sig_count_horizon = zeros(H,1);

for v = 1:length(LHSlabels)

    Yname = LHSlabels{v};

    beta = LPcoeffsAsy_LHS.beta.(Yname)(:,2);
    se   = LPcoeffsAsy_LHS.se.(Yname)(:,2);

    t_stat = beta ./ se;

    % Indicator for significance at each horizon
    sig_indicator = abs(t_stat) > z_10;

    % Add to horizon counts
    sig_count_horizon = sig_count_horizon + sig_indicator;

end

% --- Plot ---
figure
bar(0:H-1, sig_count_horizon)
ylim([0 5])
xlabel('Horizon')
ylabel('Number of significant coefficients (10%)')
title('Significant Responses Across Variables')
grid on