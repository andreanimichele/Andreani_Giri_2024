% Plot and export subplots for figure 3

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

data = xlsread('data_input/dataset.xlsx');

% Vector of quarterly numeric dates
t = data(2:end,1);        
% Extract the year
years = floor(t);
% Extract the quarter
quarters = (t - years) * 4;
% Map quarters to corresponding months (Q1 -> 01, Q2 -> 04, Q3 -> 07, Q4 -> 10)
months = (quarters * 3) + 1;
% Create the datetime vector
time = datetime(years, months, 1);

gdp = diff(log(data(:,2)));             % REAL GDP (Fred data)

mf_mort = diff(log(data(:,3)));        % Multifamily Residences

h_price = diff(log(data(:,4)));        % Average Sales Price of Houses, Dollars

% ------------------
% AWT parameters
% ------------------

dt = 1;
dj = 1/12;
low_period = 2;
up_period = 128;  
pad = 1;
mother = 'GMW';
beta = 6.0;
gamma = 0;  
sig_type = 'AR0';
alfa = 0.10;

% Window sizes for smoothing (for coherency ang gain)
windTime_size = 2;
windScale_size = 2;

% mc replications
n_sur = 50;

p = 1;
q = 1; % Surrogates based on an ARMA(1,1) model

wt_size = 0; %2
ws_size = 0; %2

% lower and upper periods for Mean Phase computation

lowperiod1 = 6;
upperiod1 = 32;

lowperiod2 = 6;
upperiod2 = 100;

lowperiod3 = 32;
upperiod3 = 100;

%%% Wavelet coherency

[WC1,WXY1,periods1,coi1,pv_WC1, WGain] = ...
 AWCOG(gdp,h_price,dt,dj,low_period,up_period,pad,mother,beta,gamma,wt_size,ws_size,n_sur);

%%% Wavelet coherency phase differences

phaseWC1_1 = MeanPHASE(WC1,periods1,lowperiod1,upperiod1); 

phaseWC1_2 = MeanPHASE(WC1,periods1,lowperiod2,upperiod2); 

phaseWC1_3 = MeanPHASE(WC1,periods1,lowperiod3,upperiod3); 

%  Computation of (mean) gains

gain1 = MeanGAIN(WGain,periods1,lowperiod1,upperiod1); 

gain2 = MeanGAIN(WGain,periods1,lowperiod2,upperiod2);

gain3 = MeanGAIN(WGain,periods1,lowperiod3,upperiod3);


%%% ---------------------------------------------------------------------
%%% PLOT Coherency, phase differences and Gain. 


x_lim = [t(1) t(end)];
x_ticks = 1950:20:2020;
y_lim = [-10,10];
y_ticks = [-8 -4 0 4 8];

logcoi = log2(coi1);
logperiods = log2(periods1);
y_limCO = [min(logperiods), max(logperiods)];

y_ticksCO_lab = [2 4 8 16 32 64 100];
y_ticksCO = log2(y_ticksCO_lab);

y_limPhase = [-pi-0.1 pi+0.1];
y_ticksPhase = -pi: pi/2 : pi;
% y_ticksPhase_lab = {'-\pi','-\pi/2','0','\pi/2','\pi'};
y_ticksPhase_lab = {'$-\pi$','$-\pi/2$','$0$','$\pi/2$','$\pi$'};

perc5 = 0.05;  % Percentil 5
perc10 = 0.10;  % Percentil 10

%
close all;

% Phase difference, GDP vs Mortgages
figure(1);
plot(time,phaseWC1_1,'k--','LineWidth',2);  hold on;
plot(time,phaseWC1_2,'r--','LineWidth',2);  hold on;
plot(time,phaseWC1_3,'k','LineWidth',2);  hold on;
set(gca,'YGrid','on','YLim',y_limPhase,'YTick',y_ticksPhase,...
    'YTickLabel',y_ticksPhase_lab)
set(gca, 'Fontname', 'times','FontSize',20);
set(gca, 'box','off');
xtickangle(45);
switch savefigs
    case 1
        print('-depsc2', '-loose', 'figures/figure_9_left.eps');
end
hold off;

% Wavelet gain, Partial 1
figure(2);
plot(time,gain1,'k--','LineWidth',2); hold on;
plot(time,gain2,'r-.','LineWidth',2); hold on;
plot(time,gain3,'k','LineWidth',2); axis tight;
% recessionplot;
set(gca, 'box','off');
xtickangle(45);
set(gca, 'Fontname', 'times','FontSize',20);
legend('(6,32) quarters', '(6,80) quarters',...
    '(32,80) quarters', 'Location','northwest');
switch savefigs
    case 1
        print('-depsc2', '-loose', 'figures/figure_3ii.eps');
end
hold off;
