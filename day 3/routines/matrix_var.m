% prepares matrices for estimating a VAR
function [y_cut,xx]=matrix_var(yy,pp,TT,nn)
% 	// INPUT : yy (T by n matrix of variables)
% 	//		   pp  # of lags (scalar)
% 	//         TT = size of the estimation window (scalar)
% 	//         nn = number of vars in the VAR (scalar)
% 	// OUTPUT: y_cut left hand side variables (T-p by n matrix)
% 	//		   xx   = [y_(t-1),y_(t-2),...,y_(t-p) 1] is a (T-p by n*pp+1) matrix

k = nn*pp+1; % number of parameters in each equation
xx  =zeros(TT-pp,nn*pp+1);

for ll=1:pp
    zz = yy(pp-ll+1:TT-ll,:);
    xx(:,(ll-1)*nn+1:ll*nn) = zz;
end
xx(:,k) = ones(TT-pp,1);
y_cut   = yy(pp+1:TT,:);