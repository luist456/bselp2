function [YY, REGR] = make_regressors_func(Z,Y_LABEL,p_y,p_x)

clear YY REGR

% -------------------------------------------------------------------------
% EA indicators: FC control added to right-hand side using spread(:,i)
if strcmp(Y_LABEL,'URATE')
    YY              = Z.URATE;
    XXX             = [Z.URATE Z.Dl_EMP Z.Dl_IP Z.Dl_CPI Z.TBILL1Y ...
                        Z.EBP];
    REGR            = build_REGR_2lags(XXX, p_y, p_x);

elseif strcmp(Y_LABEL,'IP')
    YY              = Z.IP;
    XXX             = [Z.Dl_IP Z.URATE Z.Dl_EMP Z.Dl_CPI Z.TBILL1Y ... 
                      Z.EBP];
    REGR            = build_REGR_2lags(XXX, p_y, p_x);

elseif strcmp(Y_LABEL,'EMP')
    YY              = Z.EMP;
    XXX             = [Z.Dl_EMP Z.URATE Z.Dl_IP Z.Dl_CPI Z.TBILL1Y ...
                       Z.EBP];
    REGR            = build_REGR_2lags(XXX, p_y, p_x);

elseif strcmp(Y_LABEL,'CPI')
    YY              = Z.CPI;
    XXX             = [Z.Dl_CPI Z.Dl_EMP Z.Dl_IP Z.TBILL1Y ...
                        Z.EBP];
    REGR            = build_REGR_2lags(XXX, p_y, p_x);

elseif strcmp(Y_LABEL,'EBP')
    YY              = Z.GZEBP;
    XXX             = [Z.GZEBP Z.urgap Z.Dl_EMP Z.Dl_IP Z.Dl_CPI ...
                      Z.TBILL1Y];
    REGR            = build_REGR_2lags(XXX, p_y, p_x);

elseif strcmp(Y_LABEL,'BAA')
    YY              = Z.BAA10YM;
    XXX             = [Z.BAA10YM Z.URATE Z.Dl_EMP Z.Dl_IP Z.Dl_CPI ...
                      Z.TBILL1Y];
    REGR            = build_REGR_2lags(XXX, p_y, p_x);


elseif strcmp(Y_LABEL,'FCI')
    YY              = Z.CHIFCI;
    XXX             = [Z.CHIFCI Z.URATE Z.Dl_EMP Z.Dl_IP Z.Dl_CPI ...
                      Z.TBILL1Y];
    REGR            = build_REGR_2lags(XXX, p_y, p_x);

elseif strcmp(Y_LABEL,'GZS')
    YY              = Z.GZSPR;
    XXX             = [Z.GZSPR Z.URATE Z.Dl_EMP Z.Dl_IP Z.Dl_CPI ...
                      Z.TBILL1Y];
    REGR            = build_REGR_2lags(XXX, p_y, p_x);

else
    error('WRONG LABEL HERE')
end

end % <-- end main function


% ========================================================================
% Helper: builds regressor matrix with different lag order for col 1 vs rest
% Output columns: [const, lags(col1,1..p_y), lags(col2,1..p_x), ..., lags(coln,1..p_x)]
function REGR = build_REGR_2lags(XXX, p_y, p_x)

[T,n] = size(XXX);

K = 1 + p_y + (n-1)*p_x;
REGR = NaN(T, K);

col = 1;
REGR(:,col) = 1;
col = col + 1;

% Lags for "dependent variable regressor" (first column)
for L = 1:p_y
    REGR((L+1):T, col) = XXX(1:(T-L), 1);
    col = col + 1;
end

% Lags for all other controls (columns 2..n)
for j = 2:n
    for L = 1:p_x
        REGR((L+1):T, col) = XXX(1:(T-L), j);
        col = col + 1;
    end
end

end
