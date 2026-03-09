function irfwold = woldirf(beta, c, p, horizon)

% Function to compute the matrices of the Wold representation of a
% stationary VAR(P)

% Inputs:  beta = NP+1 x N matrix of estimated coefficients 
%          c = 1 if constant required
%          horizon = how many horizon of the Wold matrices (+1) will be
%          computed

% Outputs: irfwold = N x N x horizon + 1 array of the Wold coefficient
%          matrices

[BigA, N] = companionMatrix(beta, c, p); % Take the companion matrix 

irfwold = zeros(N, N, horizon+1); % Initialize matrix of results

for h = 1:horizon+1 % For every horizon
    
    temp = BigA^(h-1); % Compute the powers of the companion matrix
    irfwold(:, :, h) = temp(1:N, 1:N); % Take only firs N rows and N columns

end

end
