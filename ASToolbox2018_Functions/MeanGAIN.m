function mean_gain = MeanGAIN(WG,periods,low_filter_period,up_filter_period)%MEANGAIN Mean (in a given frequency band) of wavelet (or partial wavelet) gain %       % mean_gain = MeanGAIN(WG,periods,low_filter_period,up_filter_period)%	Computation of mean gain or mean partial wavelet gain (over a selected%   frequency band defined by low_filter_period and up_filter_period).%   Matrix of (complex) wavelet gain or partial wavelet gain has to be %   computed  first, by using AWCOG or MPAWCOG.%	The vector of periods (periods) constructed when using AWCOG/MPAWCOG %   is needed.%%  **********************************************************************%%   INPUTS:%      WG - WGain: matrix of (complex) gain  %                  (Output of AWCOG)%                     or %           WPGain: matrix of (complex) partial wavelet gain %                   (Output of MPAWCOG)%          %      periods - vector of the wavelet periods %                (Output of AWCOG or MPAWCOG)%  Optional inputs:%      low_filter_period - lower value of the period used for filtering %                          By default, low_filter_period = periods(1) %                          %      up_filter_period - upper value of the period used for filtering %                         By default,  up_filter_period = periods(end)%                         %                        %   OUTPUT:      %       mean_gain - mean of wavelet gain or partial wavelet gain %                  (relative to selected periods)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Example:% see Example4_Gains in folder Examples.                           %  See also: AWCOG, MPAWCOG.% Copyright 2018, L. Aguiar-Conraria and M.J. Soares%   This software may be used, copied, or redistributed as long as it is %   not sold and this copyright notice is reproduced on each copy made. %   This routine is provided as is without any expres or implied %   warranties whatsoever. %   Please acknowledge the use of this software in any publications:%   Wavelet gain software was provided by %   L. Aguiar-Conraria and M.J. Soares and is available at URL: %   https://sites.google.com/site/aguiarconraria/joanasoares-wavelets%   Please send a copy of such publications to either L. Aguiar-Conraria or%   M. J. Soares%   Lu�s AGUIAR-CONRARIA              Maria Joana SOARES                      %   Dep. Economics                    Dep. Mathematics and Applications   %   University of Minho               University of Minho%   4710-057 Braga                    4710-057 Braga%   PORTUGAL                          PORTUGAL                          %   lfaguiar@eeg.uminho.pt            jsoares@math.uminho.pt %   References:%   [1] Aguiar-Conraria, L. and Soares, M.J. (2014)%       "The Continuous Wavelet Transform: Moving beyond uni- and bivariate%       analysis", Journal of Economic Surveys, 28(2), 344-375.%   [2] Aguiar-Conraria, L., Martins, M.M. and Soares, M.J., "Estimating %       the Taylor Rule in the Time-Frequency Domain" (2017), NIPE Working %       Papers NIPE-WP ???%   [3] Cazelles, B., Chavez, M., McMichael, A. J., and Hales, S. (2005)%       "Nonstationary Influence of El Ni�o on the Synchronous Dengue %       Epidemics in Thailand", PLoS Med., 2 (4), 313-318%   [4] Cazelles, B., Chavez, M., de Magny, G. C., Gu�gan, J-F and Hales, %       S.(2007), "Time-Dependent Spectral Analysis of Epidemiological %       Time-Series with Wavelets", Journal of the Royal Society Interface, %       4, 625�36. %   [5] Lilly, J.M. and Olhede, S. C.(2009), "Higher-Order Properties of %       Analytic Wavelets", IEEE Trans. Signal Process. 57 (1), 146-160.%   [6] Lilly, J.M. and Olhede, S. C.(2010), "On the Analytic Wavelet %       Transform", IEEE Transactions on Information Theory, 56 (8), %       4135-4156.%   [7] Torrence, C. and Compo, T.C.(1998), "A Prectical Guide to Wavelet %       Analysis", Bulletin of the American Meteorological Society, 79, %       605�618.%%%%%%%%%%%%%%%     Tests on inputs and default values    %%%%%%%%%%%%%%%%% % Test on inputs if nargin < 2    error('Must input WGain or PWGain and periods') end % Default parameters if ( nargin <3 || isempty(low_filter_period)||low_filter_period == 0 )    low_filter_period = periods(1); elseif (low_filter_period < periods(1))        disp('low_filter_period too low; first period used')        low_filter_period=periods(1); end if ( nargin<4 || isempty(up_filter_period) || up_filter_period == 0 )   up_filter_period = periods(end); elseif (up_filter_period > periods(end))    disp('up_filter_period too high; last period used')     up_filter_period = periods(end); end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    ind_sel_periods = ((periods >= low_filter_period)&(periods <= up_filter_period));                                             % Indexes of selected periods  mean_gain = abs(mean(WG(ind_sel_periods,:)));                                                  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                        end % END of FUNCTION