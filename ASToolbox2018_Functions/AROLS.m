function [alfa,residuals,var_residuals] = AROLS(x,p)
%AR model estimate (based on Ordinary Least Squares)
%
% [alfa,residuals,var_residuals] = AROLS(x,p)
%  Computes an AR(p) model 
%       x[t] = alfa(1)x[t-1]+ ... + alfa(p)x[t-p] + e[t]
%  by using Ordinary Least Squares.
%  NOTE: We do not include a constant term in the regression.
%
%  INPUTS: 
%       x - a vector (the series to be modeled).
%   Optional input:
%       p - a non-negative integer, order of the AR process (Default: 1).
%	OUTPUTS:
%       beta - a column vector of length p with the parameters' estimates
%       residuals - the vector of residuals
%       var_residuals - variance of residuals
%       
 
%   Written by:
%
%   Luís AGUIAR-CONRARIA              Maria Joana SOARES                      
%   Dep. Economics                    Dep. Mathematics and Applications   
%   University of Minho               University of Minho
%   4710-057 Braga                    4710-057 Braga
%   PORTUGAL                          PORTUGAL
%                           
%   lfaguiar@eeg.uminho.pt            jsoares@math.uminho.pt 

%--------------- Test on inputs and default parameters ----------------------
 if ( nargin < 1 )
    error('Must input series');
 elseif ( nargin < 2 )
    p = 1;
 end
 N = length(x);
 if (N < p || ~(isvector(x)) )
   error('x must be vector with length greater than order of model');
 elseif ( ~( p == round(p)) || p<0 )
   error('Model order must be a non-negative integer')
 end
%-------------------   MAIN COMPUTATIONS        --------------------------
 x = x(:); % Make x a column vector
 x = x - mean(x);  % Subtract mean to x
 N = length(x);
 if p == 0 % AR(0)
    alfa = []; 
    residuals = x;
    var_residuals = var(x);
 else
    y = x(p+1:end); % AR(p), p>0
    X = zeros(N-p,p);
    for jCol = 1:p
        X(:,jCol) = x(p+1-jCol:end-jCol);
    end   
    alfa = pinv(X)*y; % Parameters' estimates     
                      % Uses matlab function pinv 
                      % to compute pseudo-inverse
    residuals = y - X*alfa; % Residuals
    var_residuals = var(residuals); % Estimate of variance of residuals
 end
end %----------------- END OF FUNCTION -----------------------------------

