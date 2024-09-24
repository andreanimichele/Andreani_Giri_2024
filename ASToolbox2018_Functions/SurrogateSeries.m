function matrix_sur = SurrogateSeries(x,n_sur,p,q)
% Computes surrogates of a series 
%   matrix_sur = SurrogateSeries(x,n_sur,p,q) 
%   Computes a matrix with n_sur surrogate series of a given series x.
% 
%   The surrogates are constructed by first fitting an ARMA(p,q) model to 
%   the series x and then constructing new samples by drawing errors 
%   from a Gaussian distribution.
%   
%   It uses the function ARMA2SR (ARMA estimate with a two-step regression) 
%   of this toolbox to fit the ARMA(p,q) model to x and the matlab function 
%   FILTER to construct the surrogates. 
%   
%   INPUTS: 
%       x - a vector (the raw series)
%   Optional inputs: 

%       n_sur - positive integer, the number of surrogates (Default: n_sur =1)
%       p,q - integers, orders of the ARMA process (Defaults: p = 0, q = 0)
%   OUTPUT:
%       matrix_sur - the matrix of size (length(x) x n_sur) whose columns 
%                    contain the surrogates of x.   
%   Copyright 2018, L. Aguiar-Conraria and M.J. Soares
%
%   Luís AGUIAR-CONRARIA              Maria Joana SOARES                      
%   Dep. Economics                    Dep. Mathematics and Applications   
%   University of Minho               University of Minho
%   4710-057 Braga                    4710-057 Braga
%   PORTUGAL                          PORTUGAL
%                           
%   lfaguiar@eeg.uminho.pt            jsoares@math.uminho.pt 

%%%%%%%%%%%%%%%%%%%% Test on inputs and default parameters %%%%%%%%%%%%%%%%

if ( nargin < 1 )
    error('Must input series');
end
if ( nargin < 2 || isempty(n_sur) ), n_sur = 1; end
if ( nargin < 3 || isempty(p) ), p = 0; end
if ( nargin < 4 || isempty(q) ), q = 0; end    
if ( ~(isvector(x)) )
    error('x must be a vector')
end

if ( ~(p == round(p)) || p<0 )
	 error('Model order p must be a non-negative integer')
end
if ( ~(q == round(q)) || q<0 )
	 error('Model order q must be a non-negative integer')
end
N = length(x);
if ( N < p || N < q )
	error('Length of series must be greater or equal to p and q')
end

%------------------------------------------------------------------------
x=x(:); % Make x a column vector

mx = mean(x); % Mean of x
sx = std(x);  % Standard deviation of x
Mx = repmat(mx,N,n_sur); % Matrix of size of surx, filled with the value mx
Sx = repmat(sx,N,n_sur); % Matrix of size of surx, filled with the value sx

[alfa,beta,~,sigmaSquared] = ARMA2SR(x,p,q); % Use of function ARMA2SR of this 
                                             % toolbox to fit an ARMA(p,q) model
                                             % to x
err = sqrt(sigmaSquared)*randn(N,n_sur);  
beta = [1 beta'];
alfa = [1 -alfa'];
matrix_sur = filter(beta,alfa,err);
% Make surrogates have the same mean and variance of original series
matrix_sur = ProcessMatrix(matrix_sur,0,1); % Normalize matrix of surrogates
matrix_sur = matrix_sur.*Sx+Mx; % Make surrogates have the same mean and  
                                    % standard deviation as x
end

