function [WCO,WCross,periods,coi,pv_WCO,WGain]=...
           AWCOG(x,y,dt,dj,low_period,up_period,pad,mother,beta,gamma,...
                 wt_size,ws_size,n_sur,p,q)                    
%AWCOG Wavelet coherency and wavelet gain of two series
%
% [WCO,WCross,periods,coi,pv_WCO,WGain]=...       
%   AWCOG(x,y,dt,dj,low_period,up_period,pad,mother,beta,gamma,...
%               wt_size,ws_size,n_sur,p,q) 
%
%   Computes the (matrix of complex) wavelet coherency (WCO) of two series 
%   x and y, using (analytic) wavelet transforms with a Morlet or a 
%   Generalized Morse Wavelet.  Also computes the (matrix of normalized)
%   cross-wavelet transform (WCross) and  the (matrix of complex) wavelet 
%   gain (WGain).
%   It  can also compute the  p-values for the coherency (pv_WCO). 
%   These are computed by Monte-Carlo simulation with n_sur surrogate series.
%   The surrogates are obtained by fitting an ARMA(p,q) model to the series
%   and constructing new samples by drawing errors from a Gaussian distribution.
%   These surrogates are constructed with the use of the function SurrogateSeries
%   contained in  the folder Auxiliary.
%   It has subfunctions: WaveletCoherGain, WaveletTransform and Hamming Window
%
%-------------------------------------------------------------------------
%         MAIN DIFFERENCES from corresponding AWCO (2014 version)
%-------------------------------------------------------------------------
%   1. Also computes wavelet gain
%   2. Allows other choices of boundary conditions (not only pad with zero,
%       but also reflexive or constant) when computing the wavelet
%       transforms
%   3. Series are not normalized (just subtract mean)
%   4. Cross is normalized (to avoid bias); see [7]
%   5. Smoothed-cross and scales are not output 
%   6. Does not make use of AWT (it is self-contained)
%   7. Default values may be chosen with 0 (zero) or [] (empty) as input
%   8. Hamming windows are used for smoothing (no choice of window and no 
%      need of Signal Processing  toolbox
%   9. There is no need of Econometrics toolbox (even if we want to use an 
%       ARMA(p,q),q>0, as null)

%--------------------------------------------------------------------------
%   INPUTS:
%       x, y - two vectors of the same length (series).
%
%   Optional inputs: 
%       dt - sampling rate 
%            0 or [] -> Default: dt = 1
%       dj - frequency resolution (1/dj = number of voices per octave)
%            0 or [] -> Default: dj = 1/4
%       low_period - lower period of the decomposition 
%            0 or [] -> Default: low_period = 2*dt
%       up_period - upper period of the decomposition 
%             0 or [] -> Default: up_period = length(x)*dt      
%       pad - the type of pad to the series (boundary conditions):
%             0 or 'zero' or []-> zero-pading (default)
%             1 or 'reflexive' -> reflection at the boundaries
%             2 or 'constant' -> constant (repeat first and last values)
%       mother - the type of mother wavelet function: 
%               0 or [] or 'Morlet'  -> Morlet wavelet (default)
%               1 or  'GMW'   -> Generalized Morse Wavelet    
%       beta - the beta parameter for the GMW or the Morlet (w0) parameter:
%             0 or []-> Default: 6.0 for 'Morlet'; 3.0 for 'GMW'                                                      
%       gamma - the gamma parameter of the GMW 
%             0 or []->  Default: 3.0           
%       wt_size - to determine the size of the Hamming window used for smoothing
%                 in the time direction; sizevaries with scale s and
%                 is given by:  wt_size*s/dt (with a minimum value of 5)
%             0 or [] -> Default: wt_size = 1;
%       ws_size - to determine the size of the Hamming window window used for smoothing
%                 in the scale direction; size depends on dj and is 
%                 given by: ws_size/(2*dj) (with a minimum value of 5)
%            0 or []-> Default: ws_size  = 1;          
%       n_sur - integer, number of surrogate series, if we want to compute
%               p-values for the Wavelet Coherency
%            0 or []-> Default: no computation               
%       p, q -  non-negative integers, orders of the ARMA(p,q) model 
%              0 or [] -> Default: p=0; q=0   
%   OUTPUTS:
%       WCO - (Complex) wavelet coherency (matrix) of x and y
%       WCross - Cross-wavelet transform (matrix)of x and y (normalized, 
%               to avoid bias) 
%       periods - The vector of Fourier periods (in time units)
%                 that correspond to the scales used
%    Optional outputs
%       coi - The "cone-of-influence", which is a vector with the same length
%             as  x  that contains the limit of the region where the 
%             wavelet transforms are influenced by edge effects
%       pv_WCO - p-values for wavelet coherency
%       WGain - (Matrix of complex) gain of x over y
%
%  Example:
%      see Example2_Coherency and Example4_Gains in the folder Examples.
%                           
%  See also: MeanPHASE and MeanGAIN
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

%   References:
%   
%   [1] Aguiar-Conraria, L. and Soares, M.J. (2014)
%       "The Continuous Wavelet Transform: Moving beyond uni- and bivariate
%       analysis", Journal of Economic Surveys, 28(2), 344-375.
%   [2] Aguiar-Conraria, L., Martins, M.M. and Soares, M.J. (2017), 
%      "Estimating the Taylor Rule in the Time-Frequency Domain",
%        NIPE - WP 04/2018 
%   [3] Cazelles, B., Chavez, M., de Magny, G. C., Guégan, J.-F. and 
%       Hales, S. (2007), "Time-Dependent Spectral Analysis of
%       Epidemiological Time-Series with Wavelets", Journal of the Royal 
%       Society Interface, 4, 625—36. 
%   [4] Lilly, J.M. and Olhede, S. C.(2009), "Higher-Order Properties of 
%       Analytic Wavelets", IEEE Trans. Signal Process. 57 (1), 146-160.
%   [5] Lilly, J.M. and Olhede, S. C.(2010), "On the Analytic Wavelet 
%       Transform", IEEE Transactions on Information Theory, 56 (8), 
%       4135-4156.
%   [6] Torrence, C. and Compo, T.C., "A Prectical Guide to Wavelet 
%       Analysis" (1998), Bulletin of the American Meteorological Society,      
%       79, 605—618.
%   [7] Veleda, D., Montagne, R. and Araujo,M. (2012)  "Cross-wavelet Bias
%         Corrected by Normalizing Scales", Journal of Atmospheric and 
%         Oceanic Technology, 29, 1401--1408.
% 

% --------------- Tests on inputs and default values ----------------------
% Tests on inputs
 if (nargin <2)
    error('You must input two series');
 % Check that x and y are vectors
 elseif (~(isvector(x)) || ~(isvector(y))) 
    error('x and y must be vectors');
 end
 % Check that vectors have the same length
 lx = length(x);
 ly = length(y);
 if lx ~= ly
    error('Series are not of the same length');
 end
 % Default parameters 
 if (nargin<3 || isempty(dt) || dt==0),  dt = 1; end
 if (nargin<4 || isempty(dj) || dj==0),  dj = 1/4; end
 if (nargin<5 || isempty(low_period) || low_period==0), low_period = 2*dt; end
 if (nargin<6 || isempty(up_period)|| up_period==0), up_period = lx*dt; end
 if (nargin<7 || isempty(pad)), pad = 0; end                                                                                           % 
 if (nargin<8 || isempty(mother)), mother = 0; end
 if (nargin<9 || isempty(beta)), beta = 0; end
 if (nargin<10 || isempty(gamma) || gamma==0), gamma = 3.0; end
 if (nargin<11 || isempty(wt_size)|| wt_size==0), wt_size = 1; end 
 if (nargin<12 || isempty(ws_size) ||  ws_size==0), ws_size = 1; end 
 if (nargin<13 || isempty(n_sur)), n_sur = 0; end
 if (nargin<14 || isempty(p)), p = 0; end
 if (nargin<15 || isempty(q)), q = 0; end
 
 %------------------------------------------------------------------------
 if ischar(pad)
    if strcmpi(pad,'zero')
        pad = 0;
    elseif strcmpi(pad,'reflexive')
        pad = 1;
    elseif strcmpi(pad,'constant')
        pad = 2;
    else
        disp('Unknown type of padding; we will use zero-padding')
        pad = 0;
    end
 end
if ~(pad == 0||pad == 1||pad == 2)
        disp('Unknown type of padding; we will use zero-padding')
        pad = 0;
end

%--------------------------------------------------------------------------
 if ischar(mother)
    if strcmpi(mother,'Morlet')
        mother = 0;
    elseif strcmpi(mother,'GMW')
        mother = 1;
    else
        disp('Unknown mother wavelet; we will use a Morlet')
        mother = 0;
    end
 end
 if ~(mother == 0||mother == 1)
        disp('Unknown mother wavelet; we will use a Morlet')
        mother = 0;
 end
 
%-------------------------------------------------------------------------
%            This part depends on the mother wavelet 
%-------------------------------------------------------------------------
 if mother == 0 % This is for the Morlet wavelet
    if beta == 0, beta = 6.; end
    energy_f = beta; % Energy frequency (see [1, Eq. (17)])
    sigma_t = 1/sqrt(2); % Time-spread sigma_t (see [1, p. 352])
    K = pi^(-0.25);      % Normalizing constant (see [1, p.352])
    Fourier_factor = (2*pi)/energy_f; % Fourier factor  (see [1, p. 375])
 else % Use function GMWMeasures of this toolbox (in folder Auxiliary)
    if beta == 0, beta = 3.; end
    [~,~,~,~,K,~,Fourier_factor,~,sigma_t]=...
        GMWMeasures(beta,gamma);                                  
                        % Use of function GMWMeasures to obtain                                  
                        % K, Fourier_factor and sigma_t for                                   
                        % Generalized Morse Wavelet
 end
 %------------ End of part dependent on mother wavelet --------------------

 % Determine extra length to pad series (depends on pad type) 
 pot2 = nextpow2(lx);
 total_length = 2^(pot2+2);
 if pad == 0 
    extra_length = total_length-lx;
 else
    extra_length = total_length-3*lx; 
 end

 % Compute vectors of scales and periods 
 s0 = low_period/Fourier_factor; % Convert low_period to minimum scale s0
 if up_period > lx*dt
    disp('up_period is too long; it will be adapted')
    up_period = lx*dt;
 end
 max_scale = up_period/Fourier_factor; % Convert up_period to maximum scale 
 J = fix(log2(max_scale/s0)/dj);  % Index of maximum scale
 scales = s0*2.^((0:J)*dj);       % Vector of scales(see [1, Eq. (B2)])             
 periods = Fourier_factor*scales; % Conversion of scales to periods

 %*********************   MAIN COMPUTATIONS  ******************************
 %            (Uses subfunction WaveletCoherGain) 
 %
 if nargout > 5 % Compute wavelet coherency, (normalized) cross and gain
	[WCO,WCross,WGain]=...
        WaveletCoherGain(x,y,lx,dt,dj,pad,extra_length,mother,beta,gamma,...
                         K,scales,wt_size,ws_size);
 else % Compute only wavelet coherency and (normalized) cross
    [WCO,WCross]=...
        WaveletCoherGain(x,y,lx,dt,dj,pad,extra_length,mother,beta,gamma,...
                         K,scales,wt_size,ws_size);
 end
 if nargout > 3 % Compute coi 
	coiS = Fourier_factor/sigma_t; % see [1, p.37]
	coi = coiS*dt*[1E-5,1:((lx+1)/2-1),fliplr((1:(lx/2-1))),1E-5];
 end
 if nargout > 4 
    %  Compute p-values of Wavelet Coherency 
    %      (if desired; only computed if n_sur>0)                          
        if n_sur > 0
            [n_rows,n_cols] = size(WCO);
            pv_WCO=zeros([n_rows,n_cols]);
            matrix_sur_x = SurrogateSeries(x,n_sur,p,q); % 
            matrix_sur_y = SurrogateSeries(y,n_sur,p,q); % 
                                                        % Use of function
                                                        % SurrogateSeries
                                                        % (given in the folder
                                                        %  Auxiliary)      
            for i_sur = 1:n_sur      
                x_sur = matrix_sur_x(:,i_sur);
                y_sur = matrix_sur_y(:,i_sur);
                WCO_sur =...
                    WaveletCoherGain(x_sur,y_sur,lx,dt,dj,pad,extra_length,...
                                     mother,beta,gamma,K,scales,...
                                     wt_size,ws_size);
                for jj = 1:n_rows
                    for kk = 1:n_cols
                        if abs( WCO_sur(jj,kk) ) >=abs( WCO(jj,kk) )
                            pv_WCO(jj,kk) = pv_WCO(jj,kk)+1;
                        end
                    end
                end
            end
            pv_WCO = pv_WCO/n_sur; 
        else
            warning('n_sur=0; pv_WCO not computed')
            pv_WCO = [];
        end
 end
  
end  
%----------------  END OF FUNCTION AWCOG -------------------------------

%--------------------------------------------------------------------------
%                 SUBFUNCTION WaveletCoherGain
%------------------------------------------------------------------------%           
function [WCO,WCross,WGain] = ...
    WaveletCoherGain(x,y,lx,dt,dj,pad,extra_length,mother,beta,gamma,...
                     K,scales,wt_size,ws_size)

 % Compute wavelet transforms of x and y  
 WTx = WaveletTransform(x,lx,dt,pad,extra_length,mother,beta,gamma,K,scales);
 WTy = WaveletTransform(y,lx,dt,pad,extra_length,mother,beta,gamma,K,scales); 
 
 % Compute Cross-Wavelet Transform 
 WCross = WTx.*conj(WTy); %Formula (23) in [1]
 
 % Smoothing process (to compute coherency and gain)
 [n_scales,n_times]= size(WCross);
 smoothCross = zeros(n_scales,n_times);
 Wxx=abs(WTx).^2; 
 Wyy=abs(WTy).^2;
 % Smoothing in scale direction 
    ws_size = fix(ws_size/(2*dj)); %
    if ws_size < 5
        ws_size = 5; % Minimum size of Hamming window used
    end
    wind = HammingWindow(ws_size);
    for iTime = 1:n_times
        Wxx(:,iTime) = conv(Wxx(:,iTime),wind,'same');
        Wyy(:,iTime) = conv(Wyy(:,iTime),wind,'same');
        smoothCross(:,iTime) = conv(WCross(:,iTime),wind,'same');
    end   
 % Smoothing in time direction 
   for iScale = 1:n_scales                                         
        wTSS = fix(scales(iScale)*wt_size/dt); % Size of window is
                                               % adapted to scale
        if wTSS < 5
            wTSS = 5; % Minimum size of Hamming window used
        end
        wind = HammingWindow(wTSS);
        Wxx(iScale,:) = conv(Wxx(iScale,:),wind,'same');
        Wyy(iScale,:) = conv(Wyy(iScale,:),wind,'same');       
        smoothCross(iScale,:) = conv(smoothCross(iScale,:),wind,'same');
   end
% Compute (Complex) Wavelet Coherency and (Complex) Gain 
 WCO = smoothCross./sqrt(Wxx.* Wyy); % Formula (25) of [1]
 WGain = smoothCross./Wyy; % See [2,p.8]
end
%---------------- END OF FUNCTION WaveletCoherGain---------------------

%-------------------------------------------------------------------------
%                  SUBFUNCTION WaveletTransform                          
%-----------------------------------------------------------------------
function WT = ...
        WaveletTransform(x,lx,dt,pad,extra_length,mother,beta,gamma,K,scales)                                        
                                
 x = x(:).';    % Make x a row vector
 x = x-mean(x); % HERE WE DO NOT normalize series (just subtract mean)
 
 % Padding (boundary conditions)
 half_extra_length = floor(extra_length/2);
 if pad == 0 % Zero padding
	if rem(extra_length,2)==0
        x_extended = [zeros(1,half_extra_length),x,...
                     zeros(1,half_extra_length)];
	else
        x_extended = [zeros(1,half_extra_length),x,...
                     zeros(1,half_extra_length+1)];
	end
        
 elseif pad == 1 % Reflexive padding
    if rem(extra_length,2)==0
        x_extended = [zeros(1,half_extra_length),fliplr(x),x,...
                     fliplr(x),zeros(1,half_extra_length)];
    else
        x_extended = [zeros(1,half_extra_length),fliplr(x),x,...
                     fliplr(x),zeros(1,half_extra_length+1)];
    end
 
 else % Constant padding
    x_initial = x(1);
    x_final = x(end);     
    if rem(extra_length,2)==0
        x_extended = [zeros(1,half_extra_length),repmat(x_initial,1,lx),x,...
                     repmat(x_final,1,lx),zeros(1,half_extra_length)];
    else       
        x_extended = [zeros(1,half_extra_length),repmat(x_initial,1,lx),x,...
            repmat(x_final,1,lx),zeros(1,half_extra_length+1)];
    end
 end
 % End of padding 

 N = length(x_extended); % New data length after extending series

 % Compute angular frequencies (see [1, p.372])                          
 wk = 1:fix(N/2);
 wk = wk*((2*pi)/(N*dt));
 wk = [0., wk, -wk(fix((N-1)/2):-1:1)];

 % Compute Fast Fourier Transform of x_extended 
 ftx = fft(x_extended);

 % Compute Wavelet Transform 
 n_scales = length(scales);   % Number of scales
 WT = zeros(n_scales, N);     % Matrix to accomodate the Wavelet Transform
 WT = WT + 1i*WT;             % Make it complex

 for i_scale = 1:n_scales % Do the computation for each scale 
                          % (at all times simultaneously)
    scaleP = scales(i_scale); % Particular scale 
    norm = K*sqrt(scaleP/dt);
    if mother == 0 % This is for Morlet
        expnt = -(scaleP*wk-beta).^2/2.*(wk > 0.);
        daughter = norm*sqrt(2*pi)*exp(expnt);
        daughter = daughter.*(wk > 0.); 
        WT(i_scale,:) = ifft(ftx.*daughter); % (see [1, pp. 372-373])
    else  % This is for GMW
        expnt = -(scaleP*wk).^gamma.*(wk > 0.); % 
        daughter = norm*(scaleP*wk).^beta.*exp(expnt);
	    daughter = daughter.*(wk > 0.);  %(Formula (18) in [1]) 
        WT(i_scale,:) = ifft(ftx.*daughter);  % (see [1,pp. 372-373])   
    end
 end
 % Truncate WT to correct size
 % (delete part corresponding to zero-padding)
 if  pad==0
    WT = WT(:,half_extra_length+1:half_extra_length+lx);
 else
    WT = WT(:,half_extra_length+lx+1:half_extra_length+2*lx);
 end
end
%--------------  END OF FUNCTION WaveletTransform ---------------------

%-------------------------------------------------------------------------
%               SUBFUNCTION HammingWindow               
%------------------------------------------------------------------------
function w = HammingWindow(L)
% We only use windows of miminum size L=5 (no need to consider special case)
w = .54 - .46*cos(2*pi*(0:L-1)'/(L-1));
end
%-----------   END  of subfunction HammingWindow -----------------------
