function [ynext] = bootstrapVAR(y, p ,c, beta, residuals)

% Function to bootstrap a VAR sunig the resampling technique

% Inputs:  y = original data
%          beta = Np + 1 x N matrix fo estimated coefficients
%          c = 1 if constant required
%          residuals = T-p x N matrix of OLS residuals from VAR(p)
%          estimation

% Outputs: ynext = T x N matrix of bootstrapped time series

[T, N] = size(y);
yinit = y(1:p, :); % Starting point

if c == 1 
    const = beta(1, :);
    pi = beta(2:end, :);
else
    const = 0;
    pi = beta;
end

ynext = zeros(T, N); % Initialize of bootsrpped dataset
ynext(1:p, :) = yinit;
yinit = reshape(flipud(yinit)', 1, []); % Initial data to feed the model

for i = 1:T-p
    ynext(p+i, :) = const + yinit * pi + residuals(randi(T-p),:); % Draw a random  
                    % integer to jointly select from residuals
    yinit  = reshape(flipud(ynext(i+1:p+i,:))', 1, []); % New values to feed into the model
end

end