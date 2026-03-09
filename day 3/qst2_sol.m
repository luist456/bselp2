%% Question 2.1: Replication Using Standard Local Projections

% --- Other variables ---

% --- Common inputs ---
p_y = 12; % Maximum number of lags for depedent variable
p_x = 12; % Maxiimum number of lags for other controls
H =  25; % Maximum number of horizons
hStart = 0;
c = 0;
R_pos = [1,1]; % Linear combination for positive shock
R_neg = [1,0]; % Linear combination for negative shock


for v = 1:length(LHSlabels) % Loop over all variables

    Y_LABEL = LHSlabels(v);   
    Yname   = Y_LABEL{:};     
    disp(['Making LP for: ', Yname]) % Pick the name of the current variable

    [YY, X] = make_regressors_func(Z,Y_LABEL,p_y,p_x); % Extract regressor matrix

    res_pos = lp_ols_interaction(YY, X, shock_array, H, hStart, "cum", c, R_pos);
    res_neg = lp_ols_interaction(YY, X, shock_array, H, hStart, "cum", c, R_neg);

    LPcoeffsAsy_LHS.("pos").("beta").(Yname) = res_pos.irf .* w; % Store 
    % shocks coefficients (+)
    LPcoeffsAsy_LHS.("pos").("se").(Yname)  = res_pos.se_irf .* w; % Store
    % corresponding standard errors (+)

    LPcoeffsAsy_LHS.("neg").("beta").(Yname) = res_neg.irf .* w; % Store 
    % shocks coefficients (-)
    LPcoeffsAsy_LHS.("neg").("se").(Yname)  = res_neg.se_irf .* w; % Store
    % corresponding standard errors (-)
    
end

disp('--------- Done estimating asymmetric LPs ---------')


% --- Plots ---

irf_posStruct = LPcoeffsAsy_LHS.("pos").("beta");
se_posStruct   = LPcoeffsAsy_LHS.("pos").("se");


irf_negStruct = LPcoeffsAsy_LHS.("neg").("beta");
se_negStruct   = LPcoeffsAsy_LHS.("neg").("se");

z = 1;   % 68% CI

for i = 1:numel(LHSlabels)
    
    v = LHSlabels{i};
    
    irf_pos = irf_posStruct.(v);   % H x 1
    se_pos  = se_posStruct.(v);     % H x 1

    irf_neg = - irf_negStruct.(v);   % H x 1 flipped sign for comparison
    se_neg  = se_negStruct.(v);     % H x 1
    
    H = size(irf_pos,1);
    h = (0:H-1)';

    % 68% confidence bands
    lo_pos = irf_pos - z*se_pos;
    hi_pos = irf_pos + z*se_pos;
    lo_neg = irf_neg - z*se_neg;
    hi_neg = irf_neg + z*se_neg;

    figure('Color','w'); hold on; box on;

   
    c_neg = [0 0.447 0.741];      % blue  -> negative shock
    c_pos = [0.850 0.325 0.098];  % red   -> positive shock

    % CI areas 
    fill([h; flipud(h)], [lo_pos; flipud(hi_pos)], c_pos, ...
        'FaceAlpha',0.2,'EdgeColor','none','HandleVisibility','off');
    fill([h; flipud(h)], [lo_neg; flipud(hi_neg)], c_neg, ...
        'FaceAlpha',0.2,'EdgeColor','none','HandleVisibility','off');

    % IRFs 
    p1 = plot(h, irf_pos, 'LineWidth',2, 'Color',c_pos);
    p2 = plot(h, irf_neg, 'LineWidth',2, 'Color',c_neg);

    yline(0,'k','LineWidth',1,'HandleVisibility','off');

    title(v,'Interpreter','none');
    xlabel('Horizon');
    ylabel('Response');

    legend([p1 p2], {'IRF (+)','IRF (-)'}, 'Location','best');

    grid on;
end


%% Question 2.2: Replication Using Standard Local Projections

% --- Common inputs ---
p_y = 12; % Maximum number of lags for depedent variable
p_x = 12; % Maxiimum number of lags for other controls
H =  25; % Maximum number of horizons
hStart = 0;
c = 0;
R_diff = [0,1]; % Linear combination for difference


for v = 1:length(LHSlabels) % Loop over all variables

    Y_LABEL = LHSlabels(v);   
    Yname   = Y_LABEL{:};     
    disp(['Making LP for: ', Yname]) % Pick the name of the current variable

    [YY, X] = make_regressors_func(Z,Y_LABEL,p_y,p_x); % Extract regressor matrix

    res_diff = lp_ols_interaction(YY, X, shock_array, H, hStart, "cum", c, R_diff);
    
    LPcoeffsAsy_LHS.("diff").("beta").(Yname) = res_diff.irf .* w; % Store shocks coefficients
    LPcoeffsAsy_LHS.("diff").("se").(Yname)  = res_diff.se_irf .* w; % Store corresponding standard errors
    
end

disp('--------- Done estimating asymmetric LPs ---------')


% --- Plots ---

irf_diffStruct = LPcoeffsAsy_LHS.("diff").("beta");
se_diffStruct   = LPcoeffsAsy_LHS.("diff").("se");


varNames = fieldnames(irf_diffStruct);
z = 1;   % 68% CI

for i = 1:numel(varNames)
    
    v = varNames{i};
    
    irf_diff = irf_diffStruct.(v);   % H x 1
    se_diff  = se_diffStruct.(v);     % H x 1

    
    H = size(irf_diff,1);
    h = (0:H-1)';

    % 68% confidence bands
    lo_diff = irf_diff - z*se_diff;
    hi_diff = irf_diff + z*se_diff;
    
    figure('Color','w'); hold on; box on;

    c_diff = [0.000 0.500 0.000];  % dark green
    
    % CI dashed borders only
    plot(h, lo_diff, '--', 'Color', c_diff, 'LineWidth', 1.2, ...
    'HandleVisibility','off');
    hold on
    plot(h, hi_diff, '--', 'Color', c_diff, 'LineWidth', 1.2, ...
    'HandleVisibility','off');

    % IRFs 
    p1 = plot(h, irf_diff, 'LineWidth',2, 'Color',c_diff);

    yline(0,'k','LineWidth',1,'HandleVisibility','off');

    title(v,'Interpreter','none');
    xlabel('Horizon');
    ylabel('Response');

    legend(p1, {'Difference'}, 'Location','best');

    grid on;
end
