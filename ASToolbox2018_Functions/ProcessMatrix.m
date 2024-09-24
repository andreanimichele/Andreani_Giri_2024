function mat = ProcessMatrix(mat,log_op,norm_op)%Pre-processes columns of a given matrix %%   mat=ProcessMatrix(mat,log_op,norm_op)%   Pre-process columns of matrix mat for later use.%   Pre-processing can include: %       logarithmization %       normalization%                           %   INPUTS:%        mat - a matrix (whose columns we want to process).%%   Optional inputs:%       log_op - defines option for logarithmization:%           if log_op==1 -> take log of data %                          (Not done, by default)%       norm_op - defines option for normalization:%           if norm_op==1 -> normalize data (i.e. subtract mean and divide %                                           by standard deviation)%                                           (Done, by default)                   %   OUTPUT:%       mat - processed matrix.%   Written by L. Aguiar-Conraria and M.J. Soares%%   Lu�s AGUIAR-CONRARIA              Maria Joana SOARES                      %   Dep. Economics                    Dep. Mathematics and Applications   %   University of Minho               University of Minho%   4710-057 Braga                    4710-057 Braga%   PORTUGAL                          PORTUGAL%                           %   lfaguiar@eeg.uminho.pt            jsoares@math.uminho.pt %%%%%%%%%%%%%%%%%%%%    Default values       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%if ( nargin<2 || isempty(log_op) )    log_op=0;  % By default, we do not take log of seriesendif ( nargin<3|| isempty(norm_op) )    norm_op=1; % By default, we normalizeend%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%mat(any(isnan(mat)'),:) = []; % To delete empty spacessmat = size(mat);if  smat(1)==1   % Just one row series    	mat=mat(:);  % Make it a columnendn_rows=size(mat,1);if log_op == 1 % Take log of data    if ~(all(mat(:)>0))        disp('Matrix has negative or zero elements; can not take log')    else	   mat=log(mat);     endendif norm_op==1 % Normalize data	    mean_cols = mean(mat);    mat = mat-repmat(mean_cols,n_rows,1);    std_cols = std(mat);    mat = mat./repmat(std_cols,n_rows,1); endif smat(1)==1    mat = mat(:).';end%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%end