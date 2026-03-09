function [comp, N] =  companionMatrix(beta, c, p)

% Function to create the F matrix of the companion form associated with a
% VAR(p)

% Inputs:  beta = Np+1 x N matrix of estimated coeff or Np x N if no
%          constant is included
%          c = 1 if constant required
%          p = VAR lag order

% Outputs: comp = Np x Np matrix F of the companion form 
%          N = number of dependent variables in the VAR

if c == 1 
    N = (size(beta,1)- 1) / p;
    beta = beta(2:end,:);
else 
    N = size(beta, 1) / p;
end

comp = zeros(N*p, N*p); % Initialization companion matrix
comp(1:N, :) = beta'; % Fill first N rows with the estimated coeff
comp(N+1:end, 1:N*(p-1)) = eye(N*(p-1)); % Fill the remaing part except for the last block column with ident matr
% SEE PAG 7 for more explenation

end