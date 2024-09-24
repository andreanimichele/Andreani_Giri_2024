function [MCO,PCO,periods,coi,pv_MCO,pv_PCO,PGain] =...
          MPAWCOG(X,dt,dj,low_period,up_period,pad,mother,beta,gamma,...
                  coher_type,index_partial,wt_size,ws_size,n_sur,p,q)
%MPAWCOG Multiple and/or partial wavelet coherency and partial gain of several series
%                  
%   [MCO,PCO,periods,coi,pv_MCO,pv_PCO,PGain] =...
%         MPAWCOG(X,dt,dj,low_period,up_period,pad,mother,beta,gamma,...
%                 coher_type,index_partial,wt_size,ws_size,n_sur,p,q)
%
%   Computes  wavelet multiple coherency matrix (MCO) and/or (complex)
%   wavelet partial coherency matrix (PCO) and (complex) wavelet partial 
%   gain matrix (PGain) of series given as columns of a matrix X. These 
%   are computed using analytic wavelet transforms with a Morlet or a 
%   Generalized Morse Wavelet. 
%   One can also compute the p-values for the multiple/partial wavelet 
%   coherencies (pv_MCO,pv_PCO). These are computed by Monte-Carlo simulation
%   with n_sur surrogate series.
%   The surrogates are obtained by fitting an ARMA(p,q) model to the series
%   and constructing new samples by drawing errors from a Gaussian distribution.
%   These surrogates are constructed with the use of function
%   SurrogateSeries contained in folder Auxiliary
%
%   It has subfunctions: MPWaveletCoherency, WaveletCoherGain,
%                        WaveletTransform and HammingWindow.
%
%
%------    MAIN DIFFERENCES from corresponding MPAWCO (2014 version)  -----
%   1. Also computes partial wavelet gain
%   2. Allows other choices of boundary conditions (not only pad with zero,
%      but also reflexive or constant)
%   3. Scales are not output 
%   4. Series are not normalized (just subtract mean)
%   5. Does not make use of AWCO(it is self-contained)
%   6. Default values may be chosen with 0 (zero) or [] (empty) as input
%   7. Hamming windows  are used for smoothing (no coice of windows and 
%      no need of SignalProcesing toolbox)
%   8. There is no need of Econometrics toolbox, (even if we want to use an 
%     ARMA(p,q),q>0, as null)
%--------------------------------------------------------------------------
%
%   INPUTS:
%       X - matrix with m (m>2) columns (time series) whose partial or multiple
%           coherencies we want to compute.
%           The first column X(:,1) has a special role, e.g. in case of
%           multiple coherency, we want to compute R1.(2...m)
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
%       coher_type - type of coherency we want to compute 
%               0 or [] or 'both' -> both (default)  
%               1 or  'part' -> partial only
%               2 or  'mult' -> multiple only
%       index_partial - index of series for partial coherency and 
%                       partial gain
%               0 or [] -> Default: index_partial = 2     
%       wt_size - to determine size of Hamming window used for smoothing
%                 in the time direction; size used varies with scale s and
%                 is given by:  wt_size*s/dt (with a minimum value of 5)
%             0 or [] -> Default: wt_size = 1;
%       ws_size - to determine size of Hamming window used for smoothing
%                 in the scale direction; size used depends on dj and is 
%                 given by: ws_size/(2*dj) (with a minimum value of 5)
%            0 or []-> Default: ws_size  = 1;                       
%     
%       n_sur - integer, number of surrogate series, if we want to compute
%               p-values for the Multiple/Patial wavelet Coherency/
%            0 or []-> Default: no computation 
%       p, q -  non-negative integers, orders of the ARMA(p,q) model 
%              0 or [] -> Default: p=0; q=0
%  
%   OUTPUTS:
%       MCO - Wavelet Multiple Coherency Matrix
%       PCO - (Complex) Wavelet Partial  Coherency Matrix 
%       periods - the vector of Fourier periods (in time units)
%               that correspond to the used scales
%  Optional Outputs:
%       coi - the "cone-of-influence", which is a vector of length equal 
%             to size(X,1) that contains the limit of the region where the 
%             wavelet transforms are influenced by edge effects
%       pv_MCO - matrix of the p-values for the multiple wavelet coherency
%       pv_PCO - matrix of the p-values for the partial wavelet coherency
%       PGain -  (Matrix of complex) wavelet partial gain 
%
%  Example: 
%  see Example3_MultPartCoher and Example4_Gains in folder Examples.
%
%  See also: AWCOG, MeanPHASE, MeanGAIN.
%
% Copyright 2018, L. Aguiar-Conraria and M.J. Soares

%   This software may be used, copied, or redistributed as long as it is 
%   not sold and this copyright notice is reproduced on each copy made. 
%   This routine is provided as is without any expres or implied 
%   warranties whatsoever. 
%   Please acknowledge the use of this software in any publications:
%   Multiple and partial wavelet coherency software was provided by 
%   L. Aguiar-Conraria and M.J. Soares and is available at URL: 
%   https://sites.google.com/site/aguiarconraria/joanasoares-wavelets
%   Please send a copy of such publications to either L. Aguiar-Conraria or
%   M. J. Soares 
%
%   Luís AGUIAR-CONRARIA              Maria Joana SOARES                      
%   Dep. Economics                    Dep. Mathematics and Applications   
%   University of Minho               University of Minho
%   4710-057 Braga                    4710-057 Braga
%   PORTUGAL                          PORTUGAL
%   lfaguiar@eeg.uminho.pt            jsoares@math.uminho.pt

%   References:
%   [1] Aguiar-Conraria, L. and Soares, M.J. (2014)
%       "The Continuous Wavelet Transform: Moving beyond uni- and bivariate
%       analysis", Journal of Economic Surveys, 28(2), 344-375.
%   [2] Aguiar-Conraria, L., Martins, M.M. and Soares, M.J. (2018), 
%       "Estimating the Taylor Rule in the Time-Frequency Domain",
%        NIPE - WP 04/2018 
%   [3] Lilly, J.M. and Olhede, S. C.(2009), "Higher-Order Properties of 
%       Analytic Wavelets", IEEE Trans. Signal Process. 57 (1), 146-160.
%   [4] Lilly, J.M. and Olhede, S. C.(2010), "On the Analytic Wavelet 
%         Transform", IEEE Transactions on Information Theory, 56 (8), 
%         4135-4156.
%   [5] Torrence, C. and Compo, T.C., "A Prectical Guide to Wavelet 
%       Analysis" (1998), Bulletin of the American Meteorological Society,      
%       79, 605—618."
%

%--------------- Tests on inputs and default values  ---------------------
% Tests on inputs
 if ( nargin < 1 ), error('Must input matrix X'); end
 % Default parameters
 [nr_X,nc_X] = size(X);
 if nc_X < 3
     error('Less than 3 series; use AWCOG instead')
 end
 if (nargin<2 || isempty(dt) || dt == 0), dt = 1; end
 if (nargin<3 || isempty(dj) || dj == 0), dj = 1./4; end
 if (nargin<4 || isempty(low_period) || low_period == 0), low_period = 2*dt; end
 if (nargin<5 || isempty(up_period) || up_period == 0), up_period = nr_X*dt; end
 if (nargin<6 || isempty(pad)), pad = 0; end
 if (nargin<7 || isempty(mother)), mother = 0; end 
 if (nargin<8 || isempty(beta)), beta = 0; end
 if (nargin<9 || isempty(gamma) || gamma == 0), gamma = 3.0; end
 if (nargin<10 || isempty(coher_type) ),  coher_type = 0; end 
 if (nargin<11 || isempty(index_partial) || index_partial == 0), index_partial = 2; end 
 if (nargin<12 || isempty(wt_size)|| wt_size == 0), wt_size = 1; end 
 if (nargin<13 || isempty(ws_size) || ws_size == 0), ws_size = 1; end 
 if (nargin<14 || isempty(n_sur)), n_sur = 0; end
 if (nargin<15 || isempty(p)), p = 0; end
 if (nargin<16|| isempty(q)) , q = 0; end

 %********************

%-------------------------------------------------------------------------
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
 if ~(pad == 0|| pad == 1|| pad == 2)
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
 if ischar(coher_type)
     if strcmpi(coher_type,'both')
         coher_type = 0;
     elseif strcmpi(coher_type,'part')
         coher_type = 1;
     elseif strcmpi(coher_type,'mult')
         coher_type = 2;
     else
          disp('Unknown coherency type; we compute both')
          coher_type = 0;
     end
 end
 if ~(coher_type == 0||coher_type == 1||coher_type == 2)
        disp('Unknown coherency type; we compute both')
        coher_type = 0;
 end
 
 %------------------------------------------------------------------------

 
 %-------------- This part depends on the mother wavelet ----------------
 if mother == 0 %For a Morlet wavelet
    if beta == 0, beta = 6.; end
    energy_f = beta; % Energy frequency (see [1, Eq. (17)])
    sigmaT = 1/sqrt(2); % Time-spread sigma_t (see [1, p. 352])
    K = pi^(-0.25);     % Normalizing constant (see [1, p. 352])
    Fourier_factor = (2*pi)/energy_f; % Fourier factor  (see [1, p. 375])                                      
 else % For a Generalized Morse Wavelet
    if beta==0, beta = 3.; end
    [~,~,~,~,K,~,Fourier_factor,~,sigmaT]=GMWMeasures(beta,gamma);  
                                  % Use of function GMWMeasures to obtain                                  
                                   % K, Fourier-_actor and sigmaT for                                   
                                   % Generalized Morse Wavelet
 end
 %************    End of part dependent on mother wavelet   ***************

 %***** Determination of extra-length to pad series (depends on pad) ******

  pot2 = nextpow2(nr_X);
  total_length = 2^(pot2+2);
  if pad==0
    extra_length = total_length-nr_X;
  else
    extra_length = total_length-3*nr_X; 
  end

  %************  Compute vectors of scales and periods  ************
 s0 = low_period/Fourier_factor; % Convert low_period to minimum scale 
 if up_period > nr_X*dt
    disp('upPeriod is too long; it will be adapted')
    up_period = nr_X*dt;
 end
 up_scale = up_period/Fourier_factor;  % Convert up_period to maximum scale 
 J = fix(log2(up_scale/s0)/dj); % Index of maximum scale
 scales = s0*2.^((0:J)*dj);     % Vector of scales(see [1, Eq. (B2)])             
 periods = Fourier_factor*scales;      % Conversion of scales to periods

 
  %********************** Main computations ******************************
  % ****** Compute Multiple/Partial Wavelet Coherencies and Gain ************
  %                 (Uses subfunction MPWaveletCoherGain
 if nargout > 6     
	[MCO,PCO,PGain] = ...
        MPWaveletCoherGain(X,nr_X,nc_X,dt,dj,pad,extra_length,mother,beta,...
                           gamma,K,scales,coher_type,index_partial,...
                           wt_size,ws_size);
 else
     [MCO,PCO] = ...
     MPWaveletCoherGain(X,nr_X,nc_X,dt,dj,pad,extra_length,mother,beta,...
                        gamma,K,scales,coher_type,index_partial,...
                        wt_size,ws_size); 
 end
  %******** Compute coi (if desired)  ************************************                    
 if nargout > 3 
	coiS = Fourier_factor/sigmaT; % see [1, p.37]
	coi = coiS*dt*[1E-5,1:((nr_X+1)/2-1),fliplr((1:(nr_X/2-1))),1E-5];
 end
    % ******* Compute p-values of WMCO/WPCO (if desired) ****************
    %                      (only computed if n_sur>0)                          
    
 if nargout > 4  
        if n_sur > 0                    
            if  coher_type == 0
                [nrows,ncols] = size(MCO);
                pv_PCO = zeros([nrows,ncols]); % Initialize matrix pv_WPCO
                pv_MCO = zeros([nrows,ncols]);
            elseif coher_type == 1
                [nrows,ncols] = size(PCO);
                pv_PCO = zeros([nrows,ncols]); 
            else
                [nrows,ncols] = size(MCO);
                pv_MCO = zeros([nrows,ncols]);
            end    
            
            cellSur = cell(1,nc_X); % row cell to accomodate the matrices  
                                    % of all the surrogates of each column of X
            X_sur = zeros(nr_X,nc_X);
    
            for jcol = 1:nc_X                
                    cellSur{jcol} = SurrogateSeries(X(:,jcol),n_sur,p,q);
                           % Use of function SurrogateSeries
                           % (given in the folder Auxiliary of
                           %  this toolbox);                         
                           
            end        
            for iSur = 1:n_sur
                for jcol = 1:nc_X
                    X_sur(:,jcol) = cellSur{jcol}(:,iSur);  
                    % Each column of XSur contains a surrogate of the 
                    % corresponding column of X
                end
                [WMCO_sur,WPCO_sur] =...
                    MPWaveletCoherGain(X_sur,nr_X,nc_X,dt,dj,...
                     pad,extra_length,mother,beta,...
                     gamma,K,scales,coher_type,index_partial,...
                     wt_size,ws_size);
                if coher_type == 0
                    for iRow=1:nrows
                        for iCol=1:ncols
                            if WMCO_sur(iRow,iCol) >= MCO(iRow,iCol)
                                pv_MCO(iRow,iCol) = pv_MCO(iRow,iCol)+1;
                            end
                            if abs(WPCO_sur(iRow,iCol)) >= abs(PCO(iRow,iCol))
                                pv_PCO(iRow,iCol) = pv_PCO(iRow,iCol)+1;
                            end
                        end      
                    end
                  
                elseif  coher_type == 1
                    for iRow = 1:nrows
                        for iCol = 1:ncols
                            if abs(WPCO_sur(iRow,iCol)) >= abs(PCO(iRow,iCol))
                                pv_PCO(iRow,iCol) = pv_PCO(iRow,iCol)+1;
                            end
                        end 
                    end
                else
                    for iRow = 1:nrows
                        for iCol = 1:ncols
                            if WMCO_sur(iRow,iCol) >= MCO(iRow,iCol)
                                pv_MCO(iRow,iCol) = pv_MCO(iRow,iCol)+1;
                            end
                        end                   
                    end                     
                end 
            end
                if coher_type == 0
                    pv_MCO=pv_MCO/n_sur;
                    pv_PCO=pv_PCO/n_sur;
                elseif coher_type == 1
                    pv_MCO=[];
                    pv_PCO=pv_PCO/n_sur;                    
                else
                    pv_MCO=pv_MCO/n_sur;
                    pv_PCO=[];
                end               
        else
            pv_MCO=[];
            pv_PCO=[];
        end       
  end
end
%----------------- END OF FUNCTION MPWCOG ---------------------------
%
%------------------------------------------------------------------------
%                SUBFUNCTION MPWaveletCoherGain                         %
%------------------------------------------------------------------------
%     Computation of multiple and partial wavelet coherencies            %
%
function [MCO,PCO,PGain] =...
          MPWaveletCoherGain(X,n_rowsX,n_colsX,dt,dj,pad,extra_length,...
                              mother,beta,gamma,K,scales,coher_type,...
                              index_partial,wt_size,ws_size)                                                            
  
SC = cell(n_colsX);   % Cell array of size (n_colsX X n_colsX) to accomodate 
                       % all the smoothed spectra                    
for iiRow = 1 : n_colsX
    for iiCol = iiRow : n_colsX
            smoothCross =...
               WSmoothCross(X(:,iiRow),X(:,iiCol),n_rowsX,dt,dj,pad,...
                                 extra_length,mother,beta,gamma,K,...
                                 scales,wt_size,ws_size);                                                           
            SC{iiRow,iiCol} = smoothCross;
            if ~(iiRow == iiCol)
                SC{iiCol,iiRow} = conj(smoothCross); % S is an Hermitian matrix
            end                        
    end    
end
SCross11=SC{1,1};
[scale_length,time_length] = size(SCross11);
% Computation of the cell-matrix S11, obtained by deleting 1st row and 1st
% column of CC; this matrix is always needed (for multiple or partial)
SC11 = SC(2:end,2:end); 

%************     CASE:  Multiple Coherency only  *************************
%                     R1^2.(2...m)
if coher_type == 2
    PCO = [];
    if nargout > 2
        error('Cannot compute partial gain with multiple coherency')
    end
    MatNum = zeros(n_colsX);
    MatDen = zeros(n_colsX-1);
    MCO = zeros(scale_length,time_length); % Initialize WMCO
    for nScale = 1:scale_length
        for nTime = 1:time_length
            for iR = 1:n_colsX
                for iC = 1:n_colsX
                    MatAux = SC{iR,iC};
                    MatNum(iR,iC) = MatAux(nScale,nTime);
                end
            end
            for iR = 1:n_colsX-1
                for iC = 1:n_colsX-1
                    MatAux = SC11{iR,iC};
                    MatDen(iR,iC) = MatAux(nScale,nTime);                    
                end
            end
            MCO(nScale,nTime) = ...
            1-det(MatNum)/(det(MatDen)*SCross11(nScale,nTime));
            MCO(nScale,nTime) = sqrt(abs(MCO(nScale,nTime)));     
                              % Formula (A4) of [1]                                                              
        end
    end
        
%*************   CASE:  Partial  Coherency only  **********************
%                     (rho_1j.(2...j-1 j+1  ...m))                    
elseif coher_type == 1
    MCO = [];
    PCO = zeros(scale_length,time_length); % Initialize WPCO
    PGain = zeros(scale_length,time_length); % Initialize WPGain
    MatNum = zeros(n_colsX-1);
    MatDen1 = zeros(n_colsX-1);
    MatDen2 = zeros(n_colsX-1);
    vrows = 1:n_colsX;
    vrows(index_partial) = [];
    SCJ1 = SC(vrows,2:end);
    SCJJ = SC(vrows,vrows);
    for nScale = 1:scale_length
        for nTime = 1:time_length
            for iR = 1:n_colsX-1
                 for iC = 1:n_colsX-1  
                    MatAux = SCJ1{iR,iC};
                    MatNum(iR,iC) = MatAux(nScale,nTime);
                    MatAux = SC11{iR,iC};
                    MatDen1(iR,iC) = MatAux(nScale,nTime);
                    MatAux = SCJJ{iR,iC};
                    MatDen2(iR,iC) = MatAux(nScale,nTime);
                 end
            end
            detMatDen1=det(MatDen1);
            detMatDen2=det(MatDen2);
            numerator = det(MatNum);
            denominatorCO = abs(detMatDen1*detMatDen2);  
            denominatorGain = abs(detMatDen1);
            PCO(nScale,nTime) =...
                    (-1)^index_partial*(numerator/sqrt(denominatorCO)); 
                     % Formula (A5) of [1]                                                    
            PGain(nScale,nTime) = (-1)^index_partial*(numerator/denominatorGain);                        
                     % See [2, p.]
        end
    end
%*************   CASE:  Multiple and Partial  Coherencies  ****************
%                     
else
    MatNum = zeros(n_colsX);
    MatDen = zeros(n_colsX-1);
    MCO = zeros(scale_length,time_length); % Initialize WMCO
    for nScale = 1:scale_length
        for nTime = 1:time_length
            for iR = 1:n_colsX
                for iC = 1:n_colsX
                    MatAux = SC{iR,iC};
                    MatNum(iR,iC) = MatAux(nScale,nTime);
                end
            end
            for iR = 1:n_colsX-1
                for iC = 1:n_colsX-1
                    MatAux = SC11{iR,iC};
                    MatDen(iR,iC) = MatAux(nScale,nTime);
                end
            end
            MCO(nScale,nTime) = ...
            1-det(MatNum)/(det(MatDen)*SCross11(nScale,nTime));
            MCO(nScale,nTime) = sqrt(abs(MCO(nScale,nTime)));
        end
    end
    PCO = zeros(scale_length,time_length);
    PGain = zeros(scale_length,time_length);
    MatNum = zeros(n_colsX-1);
    MatDen1 = zeros(n_colsX-1);
    MatDen2 = zeros(n_colsX-1);
    vrows = 1:n_colsX;
    vrows(index_partial) = [];
    SCJ1 = SC(vrows,2:end);
    SCJJ = SC(vrows,vrows);
	for nScale = 1:scale_length
        for nTime = 1:time_length
            for iR = 1:n_colsX-1
                for iC = 1:n_colsX-1  
                    MatAux = SCJ1{iR,iC};
                    MatNum(iR,iC) = MatAux(nScale,nTime);
                    MatAux=SC11{iR,iC};
                    MatDen1(iR,iC) = MatAux(nScale,nTime);
                    MatAux=SCJJ{iR,iC};
                    MatDen2(iR,iC) = MatAux(nScale,nTime);
                end
            end
            numerator = det(MatNum);
            detMatDen1=det(MatDen1);
            detMatDen2=det(MatDen2);
            denominatorCO = abs(detMatDen1*detMatDen2);
            denominatorGain = abs(detMatDen1);
            PCO(nScale,nTime) =...
                     (-1)^index_partial*(numerator/sqrt(denominatorCO));              
            PGain(nScale,nTime) = (-1)^index_partial*(numerator/denominatorGain);                          
        end
	end
end
end
%---------------- END OF FUNCTION MPWaveletCoherency ------------------

%---------------------------------------------------------------------------
%                 SUBFUNCTION WSmoothCross 
%          Uses subfunctions WaveletTransform and BartlettWindow
%------------------------------------------------------------------------%           
function sCross = WSmoothCross(x,y,lx,dt,dj,pad,extra_length,mother,...
                               beta,gamma,K,scales,wt_size,ws_size)
                                
    
%  Compute Wavelet Transform of x and y   
WTx = WaveletTransform(x,lx,dt,pad,extra_length,mother,beta,gamma,K,scales);
WTy = WaveletTransform(y,lx,dt,pad,extra_length,mother,beta,gamma,K,scales);
                                        
% Compute Cross-Wavelet Transform  
sCross = WTx.*conj(WTy); %  Formula (20) in [1]

% Smoothing Process  
[n_scales,n_times]= size(sCross);
Wxx = abs(WTx).^2; 
Wyy = abs(WTy).^2;
% Smoothing in scale direction
   ws_size = fix(ws_size/(2*dj)); %
    if ws_size < 5
        ws_size = 5; % Minimum size of Hamming window used
    end
    wind = HammingWindow(ws_size);
    for iTime=1:n_times
        Wxx(:,iTime) = conv(Wxx(:,iTime),wind,'same');
        Wyy(:,iTime) = conv(Wyy(:,iTime),wind,'same');
        sCross(:,iTime) = conv(sCross(:,iTime),wind,'same');
    end
    
% Smoothing in time direction
    for iScale=1:n_scales                                         
        wTSS = fix(scales(iScale)*wt_size/dt); % Size of window is
                                               % adapted to scale
        if wTSS < 5
            wTSS = 5; % Minimum size of Hamming window used
        end                                      
        wind = HammingWindow(wTSS);
        Wxx(iScale,:) = conv(Wxx(iScale,:),wind,'same');
        Wyy(iScale,:) = conv(Wyy(iScale,:),wind,'same');       
        sCross(iScale,:) = conv(sCross(iScale,:),wind,'same');
    end
%--------- End of smoothing process ----------------
end
% --------------  END OF FUNCTION WSmoothCross -------------------

%----------------------------------------------------------------------------
%                  SUBFUNCTION WaveletTransform                          
%--------------------------------------------------------------------------

function WT = ...
        WaveletTransform(x,lx,dt,pad,extra_length,mother,beta,gamma,K,scales)                                        
                                
x = x(:).';    % Make x a row vector
x = (x-mean(x)); % We do not normalize series; just subtract mean
%----------    Padding (boundary conditions) ------------
half_extra_length = floor(extra_length/2);
if pad == 0 
	if rem(extra_length,2)==0
        x_extended = [zeros(1,half_extra_length),x,...
                     zeros(1,half_extra_length)];
    else
        x_extended = [zeros(1,half_extra_length),x,...
                     zeros(1,half_extra_length+1)];
	end
        
elseif pad == 1
    if rem(extra_length,2)==0
        x_extended = [zeros(1,half_extra_length),fliplr(x),x,...
                     fliplr(x),zeros(1,half_extra_length)];
    else
        x_extended = [zeros(1,half_extra_length),fliplr(x),x,...
                     fliplr(x),zeros(1,half_extra_length+1)];
    end
else %Constant padding
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
%----------------------- End of padding  -----------------------
N = length(x_extended); % New data length after extending series

% Compute  angular frequencies  (see [1, p.372])
wk = 1:fix(N/2);
wk = wk*((2*pi)/(N*dt));
wk = [0., wk, -wk(fix((N-1)/2):-1:1)];

% Compute Fast Fourier Transform of x_etxended
ftx = fft(x_extended);

% Compute WaveletTransform   
n_scales = length(scales);   % Number of scales
WT = zeros(n_scales, N);     % Matrix to accomodate the Wavelet Transform
WT = WT + 1i*WT;             % Make it complex

for iScales = 1:n_scales % Do the computation for each scale 
                        % (at all times simultaneously)
    scaleP = scales(iScales); % Particular scale 
    norm = K*sqrt(scaleP/dt);
    if mother == 0 % This is for Morlet
        expnt = -(scaleP*wk-beta).^2/2.*(wk > 0.);
        daughter = norm*sqrt(2*pi)*exp(expnt);
        daughter = daughter.*(wk > 0.); 
        WT(iScales,:) = ifft(ftx.*daughter); % see [1, pp. 372-373]
    else  % This is for GMW
        expnt = -(scaleP*wk).^gamma.*(wk > 0.); % 
        daughter = norm*(scaleP*wk).^beta.*exp(expnt);
	    daughter = daughter.*(wk > 0.);  % formula (18) in [1]
        WT(iScales,:) = ifft(ftx.*daughter);  %ee [1, pp. 372-373] 
    end
end
% Truncate WT to correct size (delete part corresponding to zero-padding)
if pad == 0
    WT = WT(:,half_extra_length+1:half_extra_length+lx);
else
    WT = WT(:,half_extra_length+lx+1:half_extra_length+2*lx);
end
end
%------------- End of subfunction WaveletTransform -------------------

%-------------------------------------------------------------------------
%               SUBFUNCTION HammingWindow               
%------------------------------------------------------------------------
function w = HammingWindow(L)

% We only use windows of miminum size 5 (no need to consider a special case)

w = .54 - .46*cos(2*pi*(0:L-1)'/(L-1));
end
%-----------   END  of subfunction HammingWindow -----------------------
