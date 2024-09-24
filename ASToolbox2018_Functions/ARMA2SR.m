function [alfa,beta,residuals,var_residuals] = ARMA2SR(x,p,q)
% ARMA model estimate (based on a two-step regression)
%
% [alfa,beta,residuals,var_residuals] = ARMA2SR(x,p,q)
%  Estimates an ARMA(p,q) model 
%      x[t] = alfa(1)x[t-1] + ... + alfa(p)x[t-p] + e[t] + 
%              beta(1)e[t-1] + ...  + beta(q)e[t-q] 
%  by using a two-step regression. 
%  We first fit an AR(P) model to x  (P determined by using a BIC criterion)
%  and estimate the corresponding residuals; we then regress x on its p 
%  lagged values and q lagged residuals.
% 
% NOTE: Series is demeaned, i.e. we do not include a constant term in the
% model.
% 
%	INPUTS: 
%       x - a vector (the series to be modeled)
%       p - a non-negative integer, order of the AR process 
%       q - a non-negative integer, order of the MA process 
%	OUTPUTS:
%       alfa -  a column vector of length p with the alfa(AR parameters) estimates
%       beta -  a column vector with length q with the beta (MA parameters) estimates    
%       residuals - the vector of residuals 
%       (length is equal to:  length of x-max(p,q)-P, 
%        where P is number of lags used in first step)
%        var_residuals - variance of residuals 

%   Written by:
%
%   Luís AGUIAR-CONRARIA              Maria Joana SOARES                      
%   Dep. Economics                    Dep. Mathematics and Applications   
%   University of Minho               University of Minho
%   4710-057 Braga                    4710-057 Braga
%   PORTUGAL                          PORTUGAL
%                           
%   lfaguiar@eeg.uminho.pt            jsoares@math.uminho.pt 

% REFERENCES: 
% 
% [1]  E.J. Hannan and L. Kavalieris. "A Method for Autoregressive-Moving
%       Average Estimation" Biometrika, Vol 71, No2, Aug 1984.
% [2]  E.J. Hannan and A.J. McDougall. "Regression Procedures for ARMA
%       Estimation" Journal of the American Statistical Association, Vol
%       83, No 409, June 1988.
% [3]  E.J. Hannan and R. Rissanen,"Recursive Estimation of ARMA Order", 
%      Biometrika, 69, 81–94 (1986)  

if ( nargin < 3 )
    error('Must input series, p and q');
end
N = length(x);
maxpq = max(p,q);
pmax = ceil(log(N)^1.5); % Maximum value allowed for P
if (pmax <= maxpq || ~(isvector(x)) )
    error('length of x too small for orders (p,q) desired ');
elseif ( ~(p == round(p)) || p<0 )
    error('p must be a non-negative integer')
elseif ( ~(q == round(q)) || q<0 )
    error('q must be a non-negative integer')
end
x = x(:); % Make x a column vector
x = x-mean(x); % Subtract mean to x

if q == 0 % AR model; use function AROLS (given in the toolbox)
    beta = []; % 
    [alfa,residuals,var_residuals] =  AROLS(x,p); 
else % ARMA model
    % ---------------------------------------------------------------------
    % FIRST STEP:  Determine appropriate order of autoregression 
    %   (using  a Bayesian Information Criterion) and regress to compute the residuals 
    %   to be used in the second regression   
    %---------------------------------------------------------------------    
        BIC0 = Inf;
        dif_BIC = Inf;      
        P = p+q-1;           
        while(P < pmax-1 && dif_BIC >0)
            P = P+1;
            xtil = x(P+1:end);
            n = N-P;
            X = zeros(n,P);
            for jCol = 1:P
                X(:,jCol) = x(P+1-jCol:N-jCol);
            end
            alfa = pinv(X)*xtil;
            res = xtil-X*alfa;
            BIC = log(res'*res/n) + P*log(N)/N; % 
            dif_BIC = BIC0 - BIC;      
            BIC0 = BIC;
        end
    % ------------------------------------------------------------------ 
    %  SECOND STEP: Regress on p lagged values and q lagged residuals 
    %-------------------------------------------------------------------
        ytil = xtil(maxpq+1:end);
        nn = length(ytil);
        X = zeros(nn,p+q);
        for jCol = 1:p 
            X(:,jCol) = xtil(maxpq+1-jCol:N-P-jCol);
        end
        for jCol = 1:q
            X(:,p+jCol) = res(maxpq+1-jCol:N-P-jCol);
        end
        coefs = pinv(X)*ytil;
        alfa = coefs(1:p);
        beta = coefs(p+1:p+q); 
        residuals = ytil-X*coefs;
        var_residuals = var(residuals);
end
end %--------------- END OF FUNCTION ------------------------------
