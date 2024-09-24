function [peak_freq,energy_freq,inst_freq,K_peak,K_energy,...
          FF_peak,FF_energy,FF_inst,sigma_t,sigma_f,area_Heisen,duration]=...
          GMWMeasures(beta,gamma)
    
%Computes measures associated with a Generalized Morse Wavelet
%
%   [peak_freq,energy_freq,inst_freq,K_peak,K_energy,...
%     FF_peak,FF_energy,FF_inst,sigma_t,sigma_f,area_Heisen,duration]=...
%         GMWMeasures(beta,gamma)
%	Computes several quantities associated with a Generalized Morse Wavelet
%	with parameters beta and gamma.
%	These quantities are: 
%   . peak frequency, energy frequency and central instantaneous frequency; 
%   . normalizing constant for amplitude at peak frequncy equal to 2 and 
%     normalizing constant for unit energy;
%	. Fourier factors (to convert scales to Fourier periods) defined in 
%     terms of peak frequency, energy frequency or instantaneous frequency;
%	. standard deviation in time and standard deviation in frequency;
%	. Heisenberg area;
%   . non-dimensional duration; see [1]-[4] for details.
%   
%   INPUTS: 
%       beta - beta parameter of the GMW.
%       gamma - gamma parameter of of the GMW.
%
%   OUTPUTS:
%       peak_freq - peak frequecy.
%       energy_freq - energy frequency.   
%       inst_freq - central instanataneous frequency. 
%       K_peak - normalizing constant for amplitude at peak frequncy equal 
%                to 2.
%       K_energy - normalizing constant for unit energy.
%       FF_peak - Fourier-factor using peak frequency. 
%       FF_energy - Fourier-factor using energy frequency.
%       FF_inst - Fourier-factor using central instantaneous frequency.
%       sigma_t - radius (standard deviation) in time.
%       sigma_f - radius (standard deviation) in the frequency domain.
%       area_Heisen - Heisenberg area.
%       duration - non-dimensional duration.
% 
%   Example:
%       beta=3;
%       gamma=3;
%       [w_P,w_E,w_I,K_P,K_E,ff_P,ff_E,ff_I,sigma_T,sigma_w,aH,P] = ...
%         GMWMeasures(beta,gamma)
%
%   See also: AWT

%   References:
%   [1] Aguiar-Conraria, L. and Soares, M.J. (2014)
%         "The continuous wavelet transform: moving beyond the 
%         uni- and bivariate analyis", Journal of Economic Surveys
%         28(2), 344-375 
%   [2] Lilly, J.M. and Olhede, S. C.(2009) 
%       "Higher-Order Properties of  Analytic Wavelets", 
%       IEEE Trans. Signal Process. 57 (1), 146-160.
%   [3] Lilly, J.M. and Olhede, S. C.(2010), "On the Analytic Wavelet 
%         Transform", IEEE Transactions on Information Theory, 56 (8), 
%         4135-4156.
%   [4] Olhede, S.C. and Walden, A.T. (2002), "Generalized Morse Wavelets",
%       IEEE Transactions on Signal Processing, 50(1), 2661-2670
%

% Copyright 2018, L. Aguiar-Conraria and M.J. Soares

%   This software may be used, copied, or redistributed as long as it is 
%   not sold and this copyright notice is reproduced on each copy made. 
%   This routine is provided as is without any expres or implied 
%   warranties whatsoever. 
%   Please acknowledge the use of this software in any publications:
%   Wavelet software was provided by 
%   L. Aguiar-Conraria and M.J. Soares and is available at URL: 
%   https://sites.google.com/site/aguiarconraria/joanasoares-wavelets
%   Please send a copy of such publications to either L. Aguiar-Conraria or
%   M. J. Soares

%   Luís AGUIAR-CONRARIA              Maria Joana SOARES                      
%   Dep. Economics                    Dep. Mathematics and Applications   
%   University of Minho               University of Minho
%   4710-057 Braga                    4710-057 Braga
%   PORTUGAL                          PORTUGAL                          
%   lfaguiar@eeg.uminho.pt            jsoares@math.uminho.pt 

%%%%%%%%%%%%%%%%%%%%%%%    Test on inputs     %%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 2, 
    error('must input beta and gamma')
    
elseif ( beta<=0 || gamma<=0 )
    error('beta and gamma must be positive')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%    Outputs          %%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% Normalizing constant for amplitude at peak frequncy equal to 2 %%%%%
K_peak = 2*(exp(1)*gamma/beta).^(beta/gamma); % Formula (7) in [4] 

%%%%%%%%%%      Normalizing constant for unit energy         %%%%%%%%%%%%%%
r = (2*beta+1)/gamma; % This line was missing in a first version!!
K_energy = (2^((r+1)/2)*sqrt(pi*gamma))/sqrt(exp(gammaln(r)));  

%%%%%%%%%%%%           Energy frequency                       %%%%%%%%%%%%%
energy_freq = exp( gammaln((2*beta+2)/gamma)-gammaln((2*beta+1)/gamma) )*2^(-1/gamma);
                                                     % Formula (11) in [1]
%%%%%%%%%%%%           Peak frequency                       %%%%%%%%%%%%%%%
peak_freq = (beta/gamma)^(1/gamma); % See [4, p.147]

%%%%%%%%%%%%        Instantaneous  frequency                 %%%%%%%%%%%%%%
inst_freq = exp(gammaln((beta+2)/gamma)-gammaln((beta+1)/gamma)); % Formula (44) in [4]

%%%%%%%%%%%%       Fourier factor (see [1, p.12])            %%%%%%%%%%%%%%
FF_energy = (2*pi)/energy_freq; % For energy_freq 

FF_peak = (2*pi)/peak_freq;  % For peak_freq 

FF_inst = (2*pi)/inst_freq;  % For inst_freq 

%%%%%%%%%%%%%%%%%%       Radius in time          %%%%%%%%%%%%%%%%%%%%%%%%%%
sigma_t = sqrt((beta^2*gammatil((2*beta-1)/gamma)+gamma^2*...
            gammatil((2*beta+2*gamma-1)/gamma)-2*beta*gamma*...
            gammatil((2*beta+gamma-1)/gamma))/gammatil((2*beta+1)/gamma));
                        % Formula () in [1]
%%%%%%%%%%%%%%%%%%       Radius in frequency     %%%%%%%%%%%%%%%%%%%%%%%%%%
sigma_f = sqrt(2^(-2/gamma)*(exp(gammaln((2*beta+3)/gamma)-gammaln((2*beta+1)/gamma))-... 
      (exp(gammaln((2*beta+2)/gamma)-gammaln((2*beta+1)/gamma)))^2));
                        
  
area_Heisen = sigma_t*sigma_f; % Formula (13) in [1]
                        
duration = sqrt(beta*gamma);  % Formula (37) in [4] and Appendix D in [4]
                            % See comment in  [4, pg. 151]
                            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

% FUNCTION gammatil %%%%%%%%%%%%%%%%%%%%%%%%%

 function y = gammatil(x)
    y = exp(gammaln(x))./(2.^x);
 end

%%%%%%%%%%%%%%%%%%%%   END OF AUXILIARY FUNCTION gammatil %%%%%%%%%%%%%%%%%












