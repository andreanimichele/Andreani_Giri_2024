% Plot and export figure 1 (Global wavelet power spectra)

% Mortgages, House Prices, and Business Cycle Dynamics: A Medium-Run Exploration Using the Continuous Wavelet Transform
% International Review of Economics & Finance, 2024, Volume 94

% Michele Andreani
% Financial Stability Research Office
% Central Bank of Malta

% Federico Giri
% Universitá Politecnica delle Marche 
% Dipartimento di Scienze Economiche e Sociali (DiSES)

% Last edit: 24 September, 2024

clear all;
close all;
clc;

% Commands for ticks' font in LaTeX
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

folderName = 'figures';  % Specify the folder name
mkdir(folderName);
addpath(folderName);
savefigs = 1;            % Scalar: 1 = save figures in .eps
% Add to path Wavelet toolboxes:
% 1 - M. Joana Soares and L. Aguiar-Conraria - ASToolbox
% 2 - Aslak Grinsted - Cross Wavelet and Wavelet Coherence Toolbox
addpath('wavelet-coherence-master','ASToolbox2018_Functions');

% import data
data = xlsread('data_input/dataset.xlsx');

t = data(:,1);                      
% Calculate quarterly growth rates
mf_mort = diff(log(data(:,3)));        % Multifamily Residences
h_price = diff(log(data(:,4)));        % Average Sales Price of Houses

% ------------------
% AWT parameters
% ------------------

dt = 1;
dj = 1/12;
low_period = 2;
up_period = 128;  
pad = 1;
mother = 'Morlet'; % Morlet wavelet
beta = 6.0; % omega parameter for the Morlet
gamma = 0;
sig_type = 'AR0';
alfa = 0.10;


%%% MORTGAGES

%-----  Wavelet Power and corresponding p-values --------
[~,periods_m,coi_m,power_m,pv_power_m] =...
    AWT(mf_mort,dt,dj,low_period,up_period,pad,mother,beta,gamma,sig_type);

%-------   Local maxima of WPS (ridges)  -----------
max_power_m = MatrixMax(power_m,3,.2);

%-----  Global wavelet power spectrum (GWPS) ---------
GWPS_m = mean(power_m,2);

%%% HOUSE PRICES

%-----  Wavelet Power and corresponding p-values --------
[~,periods_h,coi_h,power_h,pv_power_h] =...
    AWT(h_price,dt,dj,low_period,up_period,pad,mother,beta,gamma,sig_type);

%-------   Local maxima of WPS (ridges)  -----------
max_power_h = MatrixMax(power_h,3,.2);

%-----  Global wavelet power spectrum (GWPS) ---------
GWPS_h = mean(power_h,2);

%%% PLOT RESULTS

y_ticks_power_lab = [2 4 8 16 32 64 100];
y_ticks_power = log2(y_ticks_power_lab);

logcoi_m = log2(coi_m);
logperiods_m = log2(periods_m);

logcoi_h = log2(coi_h);
logperiods_h = log2(periods_h);

% Figure for paper
figure(1);
plot(GWPS_m,logperiods_m, 'b', 'LineWidth',1); hold on;
plot(GWPS_h,logperiods_h, 'r','LineWidth',1);
set(gca,'YTick',y_ticks_power,'YTickLabel',y_ticks_power_lab,...
    'YDir','reverse')
set(gca, 'Fontname', 'times','FontSize',16);
set(gca, 'box','off');
ylabel('Period');
legend('Multifamily mortgages','House prices', 'Location', 'southeast');
switch savefigs
    case 1
        print('-depsc2', '-loose', 'figures/figure_1.eps');
end
hold off;







