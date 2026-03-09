function [beta, residuals] = VAR(y,p,c)

% Function to estimate a VAR(p) with or without a constant using OLS
% Inputs:   y = T x N matrix of endog variables
%           p = VAR lag order
%           c = 1 if constant required
% Outputs:  beta = Np+1 x N matrix of estimated coefficient or Np x N if no
%           constant is included
%           residuals = T-p x N matrix of OLS residuals

[T, N] = size(y);
yfinal = y(p+1:end, :); % Because with lags we lose first p obs

if c == 1
    X = [ones(T-p,1) lagmakerMatrix(y,p)];
else 
    X = lagmakerMatrix(y,p);
end

beta = (X'*X) \ X' * yfinal;
residuals = yfinal - X*beta;

end