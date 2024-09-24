function p = ChiSquareCDF(x,n)
% Cumulative distribution function of the chi-square distribution 
%   p = ChiSquareCDF(x,n) returns the chi-square cumulative distribution
%   function with n degrees of freedom at the value x.

if ( nargin <2 )
    error ('Have to input two parameters');
end

if ( ~(n == round(n)) || n<0 )
   error('n must be a non-negative integer')
end

p = gammainc(x/2, n/2);
