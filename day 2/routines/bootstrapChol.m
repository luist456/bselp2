function [bootchol, upper, lower, boot_beta] = ...
    bootstrapChol(y, p, c, beta, residuals, nboot, horizon, prc, cumulate)

% Function to compute bootstrapped Chol IRFs using the resampling method

% Inputs:   y = T x N matrix of original data
%           p = integer VAR lag order
%           c = 1 if constant required
%           beta = Np+1 matrix of estimated coefficients
%           residuals = T-p x N matrix of OLS residuals from VAR(p)
%           estimation
%           nboot = integer number of bootsrap iterations
%           horizon = integer for the IRFs
%           prc = integer between 0 and 100 to select size of bands

% Outputs:  bootchol = N x N x horizon+1 x nboot array of bootrsapped Wold
%           IRFs
%           upper = N x N x horizon+1 array of upper percentiles
%           lower = N x N x horizon+1 array of lower percentiles
%           boot_beta = N x Np+1 x nboot array of bootstrapped VAR coeff
%           estimates
%           cumulate = indices of variables for which we want cumulated
%           IRFs

[T, N] = size(y);
bootchol = zeros(N, N, horizon+1, nboot); % Initialize boot Wald IRFs sample
boot_beta = zeros(size(beta, 1), N , nboot); % Initialize boot VAR(p) coeff sample

for b=1:nboot

    % Bootstrap a new data set
    varboot = bootstrapVAR(y, p, c, beta, residuals);

    % Compute the new VAR coefficients and save
    [betaloop, resloop] = VAR(varboot, p, c);
    boot_beta(:, :, b) = betaloop;

    % Compute the Wold IRF and save
    wold_loop = woldirf(betaloop, c, p, horizon);
    sigma = (resloop' * resloop) ./ (T - 1 -p - N*p);
    S = chol(sigma, 'lower');
    chol_loop = choleskyIRF(wold_loop, S);

    % Cumulate where necessary
    for i = 1:N
        chol_loop(cumulate, i, :) = ...
            cumsum(squeeze(chol_loop(cumulate, i, :)), 2);
    end

    bootchol(:, :, :, b) = chol_loop;

end

up = (50 + prc/2);
low = (50 - prc/2);

% Extract the desired percentiles from the bootstrap distribution of the
% Wold IRFs
upper = prctile(bootchol, up, 4);
lower = prctile(bootchol, low, 4);

end