function cholirf  = choleskyIRF(wold, S)

% Function to compute the point estmate of the IRF of a VAR identified
% using Cholesky

% Inputs:   wold = N x N x horizon+1 array of Wold IRFs
%           S    = lower triangular matrix Cholesky factor

% Outputs:  chol = N x N x horizon+1 array of Cholesky identified IRFs

[N, ~, horizon] = size(wold);
cholirf = zeros(N, N, horizon);

for h = 1:horizon
    cholirf(:,:,h) = wold(:,:,h) * S;
end

end