function x = lagmakerMatrix(y,p)

% Function to create the antrix of regressors according to the SUR
% representation for a VAR

% Inputs:  y = T x N matrix of endog variables
%          p = VAR lag order
% Outputs: x = T-p x NP matrix of lagged dependent variables as regressors

[T, N] = size(y);

x = zeros(T-p, N*p); % Allocate output matrix
counter = 0;
for i = 1:p % For every lag
    for j = 1:N % For every variable
        counter = counter + 1;
        x(:,counter) = y(p+1-i:T-i, j); % Fill column counter accordingly
    end
end
end