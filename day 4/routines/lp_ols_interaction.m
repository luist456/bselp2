function res = lp_ols_interaction(y, X, S, H, hStart, lpMode, c, R)
% Inputs:   y: Tx1 dependent variable
%           X: Txk matrix of regressors 
%           S: TxM matrix of shocks
%           R: 1xM combination of parameters of interest
%           hStart: integer where to start the LP (usually 0 or 1)
%           lpMode: "cum" --> cumulative changes (i.e. long differences)
%                   "lev" --> levels
%           c: integer 0 --> no constant added
%                      1 --> constant added 
%
% Outputs:  res = struct()
%           res.beta = (H+1-hStart)x1 coefficients of shocks
%           res.se = (H+1-hStart)x1 standard errors of coefficients of shocks
%           res.h = (H+1-hStart)x1 horizons of IRFs
     
    no_coeffIRF = size(R,2);

    beta = NaN(H+1-hStart, no_coeffIRF); % Initialize matrix to store coefficients
    % of shocks
    se   = NaN(H+1-hStart, no_coeffIRF); % Same for corresponding standard errors

    se_irf = NaN(H+1-hStart, 1); % Resulting standard errors for IRF

    for h = hStart:H % Loop over horizons

        switch lpMode % Build dependent variable for both specifications

            case "level" % Level specification
                y_lead = y(h+1:end);
                yh = [y_lead; NaN(h,1)];

            case "cum" % Long diff specification 
                y_lead = [y(1+h:end); NaN(h,1)];       % y_{t+h}
                y_lag  = [NaN(1,1); y(1:end-1)];       % y_{t-1}
                yh = y_lead - y_lag;

            otherwise
                error("lpMode must be 'level' or 'cum'."); % Sanity check 
        end

        Z = [yh X S]; 
        good = all(~isnan(Z), 2); % Get rid of all nan values

        yreg = yh(good); % Final dependent variable

        if c == 1 % Add constant term in this case
            Xreg = [ones(size(X(good,:),1),1) X(good,:)]; % Final matrix of 
            % controls (case constant)
        elseif c == 0 % Do not add constant term in this case
            Xreg = X(good,:); % Final matrix of controls (case no constant)
        else
            error('c must be 0 or 1') % Sanity check
        end

        Sreg = S(good, :); % Final matrix of shocks

        RHS = [Xreg Sreg]; % Right hand side of LP regression (controls + shocks)

        row = h + 1 - hStart; % Current element of the IRF

                [bhat, ~, r, ~,  ~] = regress(yreg, RHS); % Use built-in
                % MATLAB function regress

                beta_h = bhat(end-no_coeffIRF+1:end);
                beta(row, :) = beta_h'; % Store coeffcients of interest (shocks)
 
                XtXinv = inv(RHS' * RHS);
                uhat = r;  % Residuals

                % ---- HC0 (White) ----
                meat = RHS' * ( (uhat.^2) .* RHS );  
                V_HC0 = XtXinv * meat * XtXinv;

                se(row, :) = sqrt(diag(V_HC0(end-no_coeffIRF+1:end, end-no_coeffIRF+1:end))).';
                % Store standard errors of coefficients of interest
                
                se_irf(row) = sqrt(R * V_HC0(end-no_coeffIRF+1:end, end-no_coeffIRF+1:end) * R');
                % Store standard errors of IRF

           
    end

    res = struct(); % Store results
    res.beta = beta; % Coefficients of interest
    res.irf = beta * R'; % IRF
    res.se   = se; % Standard errors of coefficients of interest
    res.se_irf = se_irf; % Standard errors of IRF
    res.h    = (hStart:H)'; % Horizons of IRF

end