function [WT,periods,coi,WPS,pv_WPS] =...
    AWT(x,dt,dj,low_period,up_period,pad,mother,beta,gamma,...
        sig_type,n_sur,p,q)
%AWT Analytic Wavelet Transform of a series 
% 
% [WT,periods,coi,WPS,pv_WPS] = ...
%      AWT(x,dt,dj,low_period,up_period,pad,mother,beta,gamma,...
%          sig_type,n_sur,p,q,sur_mode)
% 
%    Computes the (matrix of the analytic) Wavelet Transform (WT) of a given 
%    series (vector) x, using a Morlet or a Generalized Morse Wavelet. 
%    Decomposition is done between two periods (low_period,up_period). 
%    It also computes the (matrix of normalized) Wavelet Power Spectrum (WPS) 
%    and its p-values (pv_WPS); these p-values can be computed in three ways:
%     -- with the theoretical distribution, with AR(0) as null 
%                          or
%     -- with the theoretical distribution, with AR(1) as null 
%                          or
%     -- with Monte Carlo simulations; 
%        in this case, n_sur surrogates are  obtained by fitting an ARMA(p,q) model to the series
%        and constructing new samples by drawing errors from a Gaussian distribution.
%        These surrogates are constructed with the use of the function SurrogateSeries
%        contained in  the folder Auxiliary.
%    Has a subfunction: WaveletTransform
% 
% NOTE: The computation of pv_WPS using the theoretical distributions
%       is much faster than using Monte-Carlo simulations. 
%       However, it SHOULD ONLY BE USED WITH MORLET WAVELETS.
%
% --------   MAIN DIFFERENCES from 2014 version-------------------------
%    1. Allows other choices of boundary conditions (not only pad with zero,
%       but also reflexive or constant)
%    2. Power is normalized (to avoid bias); see [5] 
%    3. Vector of scales is not output (only periods)
%    4. Avoids repeating some computations when computing p-values by 
%       Monte-Carlo (computed outside function WaveletTransform)
%    5. Default values can be chosen with 0 (zero) or [] (empty) as input
% -------------------------------------------------------------------------
%
%    INPUTS:
%        x - vector (time-series).
% 
%    Optional inputs: 
%        dt - sampling rate 
%             0  or []-> Default: dt = 1
%        dj - frequency resolution (1/dj = number of voices per octave)
%             0  or [] -> Default: dj = 1/4
%        low_period - lower period of the decomposition 
%             0 or [] -> Default: low_period = 2*dt
%        up_period - upper period of the decomposition 
%              0 or []-> Default: up_period = length(x)*dt      
%        pad - the type of pad to the series (boundary conditions):
%              0 or 'zero' or [] -> zero-pading (default)
%              1 or 'reflexive' -> reflection at the boundaries
%              2 or 'constant' -> constant (repeat first and last values)
%       mother - the type of mother wavelet function: 
%               0 or 'Morlet' -> Morlet wavelet (default)
%               1 or  'GMW'   -> Generalized Morse Wavelet    
%       beta - the beta parameter for the GMW or the Morlet (omega_0) parameter:
%             0  or []-> Default: 6.0 for 'Morlet'; 3.0 for 'GMW'                                                      
%       gamma - the gamma parameter of the GMW 
%             0 or [] ->  Default: 3.0 
%       sig_type - how to compute levels of significance for wavelet power spectrum
%             0 or [] or 'NONE'-> no computation (default)
%             1 or 'AR0' -> theoretical distribution, AR(0) as null
%             2 or 'AR1' -> theoretical distribution, AR(1) as null
%             3 or 'MCS'->  Monte Carlo simulation
%       n_sur - number of surrogates for MC simulation
%              0 or [] -> Default: 0 (no simulation);
%       p,q -  non-negative integers, orders of the ARMA(p,q) model 
%               0 or [] ->  Default: p=0, q=0
%  OUTPUTS:
%       WT - (Complex) Wavelet Transform Matrix (num. rows = num. scales used; 
%                                                num. columns = length(x))
%       periods - the vector of Fourier periods (in time units) that
%                 correspond to the used scales
%  Optional outputs:
%       coi - the "cone-of-influence", which is a vector of the same 
%           size as x that contains the limit of the region where 
%           the wavelet transform is influenced by edge effects
%       WPS -  (Normalized) Wavelet Power Spectrum
%       pv_WPS - p-values for Wavelet Power Spectrum
%
%   Example: 
%   see Example1_Power in the folder Examples.
% 
%   See also: MeanPHASE.
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
%   [1] Aguiar-Conraria, L. and Soares, M.J. (2014)
%         "The continuous wavelet transform: moving beyond the 
%         uni- and bivariate analyis", Journal of Economic Surveys
%         28(2), 344-375.
%   [2] Lilly, J.M. and Olhede, S. C.(2009) 
%       "Higher-Order Properties of  Analytic Wavelets", 
%       IEEE Trans. Signal Process. 57 (1), 146-160.
%   [3] Lilly, J.M. and Olhede, S. C.(2010), "On the Analytic Wavelet 
%         Transform", IEEE Transactions on Information Theory, 56 (8), 
%         4135-4156.
%   [4] Torrence, C. and Compo, T.C., "A Practical Guide to Wavelet 
% 	      Analysis" (1998), Bulletin of the American Meteorological Society, 
% 	      79, 605-618.
%   [5] Liu, Y., Liang, X. S. and Weisberg, R. H., "Rectification of the 
%         bias problem in the wavelet power  spectrum" (2007), Journal. Atmos.
%         Ocean. Tech, 24, 2093-2102.
% 

% ------------- Tests on inputs and default values----------------------
% Tests on inputs
if nargin < 1, 
    error('must input series')
elseif ( ~ (isvector(x)) ) % Check that x is a vector
    error('x  must be a vector');
end
lx = length(x);
% Default parameters 
if (nargin<2 || isempty(dt) || dt==0 ),  dt = 1; end
if (nargin<3 || isempty(dj) || dj==0 ),  dj = 1./4; end
if (nargin<4 || isempty(low_period) || low_period==0), low_period = 2*dt; end
if (nargin<5 || isempty(up_period)|| up_period==0), up_period = lx*dt; end
if (nargin<6 || isempty(pad)), pad = 0; end                                                                                           % 
if (nargin<7 || isempty(mother)), mother = 0; end
if (nargin<8 || isempty(beta)), beta = 0; end
if (nargin<9 || isempty(gamma) || gamma==0), gamma = 3.0; end
if (nargin<10 || isempty(sig_type)), sig_type = 0; end
if (nargin<11 || isempty(n_sur)), n_sur = 0; end
if (nargin<12 || isempty(p)), p = 0; end
if (nargin<13 || isempty(q)), q = 0; end

%-----------------------------------------------------------------------
if ischar(pad)
    if strcmpi(pad,'zero')
        pad = 0;
    elseif strcmpi(pad,'reflexive')
        pad = 1;
    elseif strcmpi(pad,'constant')
        pad = 2;
    else
        disp('Unknown type of pad; we will use zero-padding')
        pad = 0;
    end
end
if ~(pad == 0||pad == 1||pad == 2)
    disp('Unknown type of pad; we will use zero-padding')
    pad = 0;
end
%-------------------------------------------------------------------------
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
%--------------------------------------------------------------------------
if ischar(sig_type)
    if strcmpi(sig_type,'NONE')
        sig_type= 0;
    elseif strcmpi(sig_type,'AR0')
        sig_type = 1;
    elseif strcmpi(sig_type,'AR1')
        sig_type = 2;
    elseif strcmpi(sig_type,'MCS')
        sig_type = 3;
    else
        disp('Unknown sig_type; we will not compute significance levels')
        sig_type = 0; 
    end
end
if ~(sig_type==0||sig_type==1||sig_type==2||sig_type==3)
        disp('Unknown sig_type; we will not compute significance levels')
        sig_type=0;
end

%------------------------------------------------------------------------ 
%               This part depends on the mother wavelet 
%------------------------------------------------------------------------
if mother==0 % Morlet
    if beta==0, beta = 6.; end
    energy_f = beta; % Energy frequency (Formula (17) in [1])
    sigma_t = 1/sqrt(2); % Time-spread sigma_t (see [1, p.352])
    K = pi^(-0.25);     %  (see [1, p.9])
    Fourier_factor = (2*pi)/energy_f; % Fourier factor 
                                      % (see [1, p.375])
else % Use of function GMWMeasures (in folder Auxiliary)to obtain                                   
     % K, Fourier_factor and sigma_t for                                   
     % Generalized Morse Wavelet
    if beta == 0, beta = 3.; end
    [~,~,~,~,K,~,Fourier_factor,~,sigma_t]=...
        GMWMeasures(beta,gamma);  
end
%----------  End of part dependent on mother wavelet -----------------------

% Determine extra length (depending on pad type) 
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
up_scale = up_period/Fourier_factor;  % Convert up_period to maximum scale 
J = fix(log2(up_scale/s0)/dj); % Index of maximum scale
scales = s0*2.^((0:J)*dj);     % Vector of scales (Formula (B2) in [1])             
periods = Fourier_factor*scales; % Conversion of scales to periods

% Compute Wavelet Transform (using subfunction WaveletTransform) 
WT = WaveletTransform(x,lx,dt,pad,extra_length,mother,beta,gamma,K,scales);

if nargout > 2 % Compute coi (if desired)
	coiS = Fourier_factor/sigma_t; % see [1, p.375]
	coi = coiS*dt*[1E-5,1:((lx+1)/2-1),fliplr((1:(lx/2-1))),1E-5]; 
end       
if nargout > 3 % Compute power 
    WPS = abs(WT).^2;  % Formula (6) in [1] 
    if nargout > 4 % Compute pv_WPS (theoretical distribution is for 
        %                            no normalized power)
        % ----- Compute p-values of Wavelet Power Spectrum  -------
        %        (only computed if sig_type > 0 and n_sur >0)                          
        % - if sig_type = 1 or sig_type = 2 (analytical distributions); 
        %   uses function ChiSquareCDF (in folder Auxiliary)
        % - if sig > 2, uses function SurrogateSeries (in folder Auxiliary)
        if sig_type == 0
            pv_WPS=[];    % No p-values for WPS are computed 
    
        elseif sig_type == 1  % Analytical distribution with AR(0) as null;  
                              % see [4], formula (18) with Pk=1 
                              % (formula (16) in [4], alpha=0)
                teor = 2*WPS;
                pv_WPS = 1-ChiSquareCDF(teor,2);   
   
        elseif sig_type == 2 % Analytical distribution with AR(1) as null; 
                [~,n_cols] = size(WPS);
                alpha = AROLS(x,1); % lag-1 coefficient in AR(1) model
                freqs = dt./periods;
                Pk = (1-alpha^2)./(1-2*alpha*cos(freqs*2*pi)+alpha^2); 
                                            % formula (16) in [4]
                PkM = repmat(Pk',1,n_cols); % Expand Pk into a matrix 
                teor = 2*WPS./PkM;
                pv_WPS = 1-ChiSquareCDF(teor,2);  
        else  % MC simulations  
            if n_sur==0
                warning('Should indicate a positive number of surrogates; pv_WPS not computed')
                pv_WPS=[];
            else
                [n_rows,n_cols] = size(WT);
                pv_WPS = zeros(n_rows,n_cols);
                matrix_sur_x = SurrogateSeries(x,n_sur,p,q); 
                            % Use of function SurrogateSeries
                            %  (given in folder Auxiliary of the toolbox)   
                for iSur = 1:n_sur
                    sur_x = matrix_sur_x(:,iSur);
                    wave_sur =...
                    WaveletTransform(sur_x,lx,dt,pad,extra_length,mother,...
                                     beta,gamma,K,scales);                                   
                    power_sur = (abs(wave_sur)).^2;
                    for jj = 1:n_rows
                        for kk = 1:n_cols
                            if power_sur(jj,kk) >= WPS(jj,kk)
                                pv_WPS(jj,kk) = pv_WPS(jj,kk)+1;
                            end
                        end
                    end
                end
                    pv_WPS = pv_WPS/n_sur;            
            end
        end
    end
        % Normalize power
	nC = size(WPS,2);
	WPS = WPS./repmat(scales.',1,nC);
end
end %---------------- END of function AWT -------------------------------
        
%------------------------------------------------------------------------
%                 SUBFUNCTION WaveletTransform                          
%-----------------------------------------------------------------------

function WT = ...
        WaveletTransform(x,lx,dt,pad,extra_length,mother,beta,gamma,K,scales)                                        
                                
x = x(:).';    % Make x a row vector
x = (x-mean(x))/std(x); % Normalize series
%--------   Padding (boundary conditions) --------------------
half_extra_length=floor(extra_length/2);
if pad == 0 % Zero-padding
	if rem(extra_length,2)==0
        x_extended = [zeros(1,half_extra_length),x,...
                     zeros(1,half_extra_length)];
    else
        x_extended = [zeros(1,half_extra_length),x,...
                     zeros(1,half_extra_length+1)];
	end
        
elseif pad==1 % Reflexive
    if rem(extra_length,2)==0
        x_extended = [zeros(1,half_extra_length),fliplr(x),x,...
                     fliplr(x),zeros(1,half_extra_length)];
    else
        x_extended = [zeros(1,half_extra_length),fliplr(x),x,...
                     fliplr(x),zeros(1,half_extra_length+1)];
    end
else % Constant padding
    xinitial = x(1);
    xfinal = x(end);     
    if rem(extra_length,2)==0
        x_extended = [zeros(1,half_extra_length),repmat(xinitial,1,lx),x,...
                     repmat(xfinal,1,lx),zeros(1,half_extra_length)];
    else       
        x_extended = [zeros(1,half_extra_length),repmat(xinitial,1,lx),x,...
            repmat(xfinal,1,lx),zeros(1,half_extra_length+1)];
    end
end
%*-------------------End of padding  --------------------------------------

N = length(x_extended); % New data length after extending series

% Compute angular frequencies (see [1,p. 372])
wk = 1:fix(N/2);
wk = wk*((2*pi)/(N*dt));
wk = [0., wk, -wk(fix((N-1)/2):-1:1)];

% Compute Fast Fourier Transform of x_extended 
ftx = fft(x_extended);

% Compute Wavelet Transform
n_scales = length(scales);   % Number of scales
WT = zeros(n_scales, N);   % Matrix to accomodate the Wavelet Transform
WT = WT + 1i*WT;           % Make it complex

for iScales = 1:n_scales % Do the computation for each scale 
                         % (at all times simultaneously)
    scaleP = scales(iScales); % Particular scale 
    norm = K*sqrt(scaleP/dt);  
    if mother==0 % This is for Morlet
        exponent = -(scaleP*wk-beta).^2/2.*(wk > 0.);
        daughter = norm*sqrt(2*pi)*exp(exponent);
        daughter = daughter.*(wk > 0.); 
        WT(iScales,:) = ifft(ftx.*daughter); % see [1, pp. 372-373]
    else  % This is for GMW
        exponent = -(scaleP*wk).^gamma.*(wk > 0.); 
        daughter = norm*(scaleP*wk).^beta.*exp(exponent);
	    daughter = daughter.*(wk > 0.);   % Formula (18) in [1]
        WT(iScales,:) = ifft(ftx.*daughter);  % [1, pp. 372-373]  
    end
end
%-----------     Truncate WT to correct size  ---------------------
%              (delete part corresponding to padding)
if  pad==0
    WT = WT(:,half_extra_length+1:half_extra_length+lx);
else
    WT = WT(:,half_extra_length+lx+1:half_extra_length+2*lx);
end

end %*------------   END of subfunction WaveletTransform --------------

