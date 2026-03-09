%% Step 1:
d_gdp = [NaN(1,1); diff(gdp)];
d_DurExp = [NaN(1,1); diff(DurExp)];
d_NonDurExp = [NaN(1,1); diff(NonDurExp)];
d_RealRate = [NaN(1,1); diff(RealRate)];
d_ff_rate = [NaN(1,1); diff(ff_rate)];

X_slow =  [d_gdp d_DurExp d_NonDurExp];

X_fast = d_RealRate;


%% Step 2:

p = 4; % Maximum lag of regressors

RHS = [lagmatrix(X_slow, 0:p) lagmatrix(X_fast, 1:p) ...
    lagmatrix(d_ff_rate, 1:p)];


%% Step 3:

LHS = [ff_rate gdp DurExp NonDurExp RealRate]; % List of variables

start = [0, 1, 1, 1, 0]; % Starting horizon of LP estimation


%% Step 4: 

varNames = {'FedFundrate','GDP','DurExp','NonDurExp','RealRate'};

H = 20; % Maximum horizon in IRFs

lpMode = "cum";
seType = "hac";
c = 1;

for v = 1:size(LHS,2) % Loop over variables

    res = lp_ols(LHS(:,v), RHS, d_ff_rate, H, start(v), lpMode, seType, c);
    % Store results using dynamic field name
    results.beta.(varNames{v}) = res.beta;
    results.se.(varNames{v}) = res.se; % Store standard errors
    results.h.(varNames{v}) = res.h;

end


%% Step 5: 


figure;

for v = 1:length(varNames)

    subplot(2,3,v)
    hold on

    h    = results.h.(varNames{v});
    beta = results.beta.(varNames{v});
    se   = results.se.(varNames{v});

    upper68 = beta + se;
    lower68 = beta - se;

    % 68% confidence band (shaded)
    fill([h' fliplr(h')], ...
         [upper68' fliplr(lower68')], ...
         [0.75 0.75 0.75], ...
         'EdgeColor','none', ...
         'FaceAlpha',0.4);

    % IRF line
    plot(h, beta, 'k', 'LineWidth', 2)

    % zero line
    yline(0,'--k')

    title(varNames{v})
    xlabel('Horizon')
    ylabel('Response')
    grid on
    box on

    hold off

end