function Fspx = FourierSpectrum(x,p,q,plot_op,dt,min_period,max_period)
%Estimate of the Fourier spectrum of a time series
%               
%   FSpx = FourierSpectrum(x,p,q,plot_op,dt,min_period,max_period)
%   Estimates the Fourier spectrum of a series x. The spectrum 
%   is estimated by first adjusting an ARMA(p,q) model to x. It can also 
%   plot the spectral density, with frequencies or periods as abcissas.
%   If q=0, it uses the ASToolbox function AROLS.
%   If q>0, it makes use of the function ARMA2SR of the toolbox
%   to fit the ARMA(p,q) model to x. 
%
%   INPUTS: 
%       x - a vector (the series).
%       p,q - non-negative inetgers, the orders of the ARMA(P,Q) model.
%   Optional inputs: 
%       plot_op-  defines plot option:
%           0 -> no plot (default)
%           1 -> plot with frequencies as abcissas.
%           2 -> plot with periods as abcissas.
%      If plot_op=2:
%       You MUST also input:     
%       dt - observation time step
%       and can select:
%       min_period- minimum period for plot.
%       max_period- maximum period for plot.
%   OUTPUT:	    
%       Fspx -the Fourier spectral density of x.      

%   Written by:
%
%   Luís AGUIAR-CONRARIA              Maria Joana SOARES                      
%   Dep. Economics                    Dep. Mathematics and Applications   
%   University of Minho               University of Minho
%   4710-057 Braga                    4710-057 Braga
%   PORTUGAL                          PORTUGAL
%                           
%   lfaguiar@eeg.uminho.pt            jsoares@math.uminho.pt 

%-------------- Test for inputs and default values ---------------------
if nargin < 3
    error('You have to input the series and the parameters p and q');
end
if nargin < 4
    plot_op=0;
end
if (plot_op == 2 && nargin <5)
    error('You must input time step dt')
end

if ( ~( p == round(p)) || p<0 )
   error('p must be a non-negative integer')
end
 if ( ~( q == round(q)) || q<0 )
   error('q must be a non-negative integer')
 end
%------------ Estimate ARMA(p,q) parameters -----------------
if q > 0
    [beta,alfa] = ARMA2SR(x,p,q);
else
        beta = AROLS(x,p);
        alfa = 0;
end

% --------------- Compute Power Spectrum ----------------------------- 
freqs = 0:pi/1000:pi; % Can be changed 
n_freqs = length(freqs);
Fspx = zeros(1,n_freqs); % Initialize vector of Fourier Spectrum
for iFreq = 1:n_freqs
	Aux1 = 1;
	Aux2 = 1;
	Aux3 = 1;
	Aux4 = 1;
	for k = 1:length(beta)
        Aux1 = Aux1-beta(k)*exp(-freqs(iFreq)*1i*k);
        Aux2 = Aux2-beta(k)*exp(freqs(iFreq)*1i*k);   
	end  
	for k=1:length(alfa)
        Aux3 = Aux3+alfa(k)*exp(-freqs(iFreq)*1i*k);
        Aux4 = Aux4+alfa(k)*exp(freqs(iFreq)*1i*k);
	end
    S =(Aux3*Aux4)/(Aux1*Aux2*2*pi);
    S = real(S);
    Fspx(iFreq)= S;
end
%	Output	
Fspx = Fspx/sum(Fspx); 

% ---------------------------	Plots  ------------------------------------
if plot_op==1
    plot(freqs,Fspx,'LineWidth',1.5);    
elseif plot_op==2
	periods = 1./freqs(2:end);     
    periods = periods*2*pi*dt;
    plot(periods,Fspx(2:end),'LineWidth',1.5);
    if nargin<6
        min_period = min(periods);
    end
    if nargin<7
        max_period = max(periods);
    end
    set(gca,'XLim',[min_period max_period])
end
%-----------------------------------------------------------------------
end %---------------- END of FUNCTION ---------------------------------



