% Plot and export figures in appendix.

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


X = [gdp mf_mort h_price];

% ------------------
% PAWC parameters
% ------------------

dt = 1;
dj = 1/12;
low_period = 2;
up_period = 128;
pad = 1;
mother = 'Morlet';
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

lowperiod1 = 6;         % short cycle
upperiod1 = 32;

lowperiod2 = 6;         % total cycle
upperiod2 = 100;

lowperiod3 = 32;        % medium cycle
upperiod3 = 100;

% Compute both multiple and partial coherences!
coher_type = 'both';
index_partial = 2; % GDP vs Mortgages
[WMCO,WPCO_Mort,periods,coi,pvM,pvP_Mort,PWGain_Mort] = ...
    MPAWCOG(X,dt,dj,low_period,up_period,pad,mother,beta,gamma,...
    coher_type,index_partial,windTime_size,...
    windScale_size,n_sur,p,q);

% Computation of  complex partial wavelet coherency and complex partial gain
% of FFR vs Gap (controlling for Inflation)
coher_type = 'part';
index_partial = 3; % GDP vs House Prices

% Here, multiple coherency is not necessary, i.e ~
[~,WPCO_House,~,~,~,pvP_House,PWGain_House] = ...
    MPAWCOG(X,dt,dj,low_period,up_period,pad,mother,beta,gamma,...
    coher_type,index_partial,windTime_size,...
    windScale_size,n_sur,p,q);

% Computation of mean partial phase differences (with 95% CIs)
%  and mean gains (in three different frequency bands)

% ---------------  MORTGAGES ---------------------------------------
[phaseDif1_Mort,low_phaseDif1_Mort,up_phaseDif1_Mort] =...
    MeanPHASE(WPCO_Mort,periods,lowperiod1,upperiod1,alfa);

[phaseDif2_Mort,low_phaseDif2_Mort,up_phaseDif2_Mort] =...
    MeanPHASE(WPCO_Mort,periods,lowperiod2,upperiod2,alfa);

[phaseDif3_Mort,low_phaseDif3_Mort,up_phaseDif3_Mort] = ...
    MeanPHASE(WPCO_Mort,periods,lowperiod3,upperiod3,alfa);

gain1_Mort = MeanGAIN(PWGain_Mort,periods,lowperiod1,upperiod1);

gain2_Mort = MeanGAIN(PWGain_Mort,periods,lowperiod2,upperiod2);

gain3_Mort = MeanGAIN(PWGain_Mort,periods,lowperiod3,upperiod3);

% ---------------  HOUSE PRICES ---------------------------------------
[phaseDif1_House,low_phaseDif1_House,up_phaseDif1_House] =...
    MeanPHASE(WPCO_House,periods,lowperiod1,upperiod1,alfa);

[phaseDif2_House,low_phaseDif2_House,up_phaseDif2_House] =...
    MeanPHASE(WPCO_House,periods,lowperiod2,upperiod2,alfa);

[phaseDif3_House,low_phaseDif3_House,up_phaseDif3_House] = ...
    MeanPHASE(WPCO_House,periods,lowperiod3,upperiod3,alfa);

gain1_House = MeanGAIN(PWGain_House,periods,lowperiod1,upperiod1);

gain2_House = MeanGAIN(PWGain_House,periods,lowperiod2,upperiod2);

gain3_House = MeanGAIN(PWGain_House,periods,lowperiod3,upperiod3);

%%% FIGURES

x_lim = [t(1) t(end)];
x_ticks = 1950:20:2020;
x_ticks_2 = 1950:5:2020;
y_lim = [-10,10];
y_ticks = [-8 -4 0 4 8];

logcoi = log2(coi);
logperiods = log2(periods);
y_limCO = log2([min(periods),max(periods)]);

y_ticksCO_lab = [2 4 8 16 32 64 100];
y_ticksCO = log2(y_ticksCO_lab);

y_limPhase = [-pi-0.1 pi+0.1];
y_ticksPhase = -pi: pi/2 : pi;

y_ticksPhase_lab = {'$-\pi$','$-\pi/2$','$0$','$\pi/2$','$\pi$'};
%y_ticksPhase_lab = {'-\pi','-\pi/2','0','\pi/2','\pi'};

perc5 = 5/100;  % Percentil 5
perc10 = 10/100;  % Percentil 10

%
%close all;

% Phase difference, Partial 1
figure(1);
plot(time,phaseDif1_Mort,'k--','LineWidth',2);  hold on;
plot(time,phaseDif2_Mort,'r--','LineWidth',2);  hold on;
plot(time,phaseDif3_Mort,'k','LineWidth',2);  hold on;
set(gca,'YGrid','on','YLim',y_limPhase,'YTick',y_ticksPhase,...
    'YTickLabel',y_ticksPhase_lab)
set(gca, 'box','off');
xtickangle(45);
set(gca, 'Fontname', 'times','FontSize',20);
switch savefigs
    case 1
        print('-depsc2', '-loose', 'figures/figure_10_left.eps');
end
hold off;

% Phase difference, Partial 2
figure(2);
plot(time,phaseDif1_House,'k--','LineWidth',2);  hold on;
plot(time,phaseDif2_House,'r--','LineWidth',2);  hold on;
plot(time,phaseDif3_House,'k','LineWidth',2);  hold on;
set(gca,'YGrid','on','YLim',y_limPhase,'YTick',y_ticksPhase,...
    'YTickLabel',y_ticksPhase_lab)
set(gca,'FontSize',16);
set(gca, 'box','off');
xtickangle(45);
set(gca, 'Fontname', 'times','FontSize',20);
switch savefigs
    case 1
        print('-depsc2', '-loose', 'figures/figure_11_left.eps');
end
hold off;

% Wavelet gain, Partial 1
figure(3);
plot(time,gain1_Mort,'k--','LineWidth',2); hold on;
plot(time,gain2_Mort,'r--','LineWidth',2); hold on;
plot(time,gain3_Mort,'k','LineWidth',2); axis tight;
% recessionplot;
set(gca, 'box','off');
xtickangle(45);
set(gca, 'Fontname', 'times','FontSize',20);
legend('(6,32) quarters', '(6,80) quarters',...
    '(32,80) quarters', 'Location','best');
switch savefigs
    case 1
        print('-depsc2', '-loose', 'figures/figure_10_right.eps');
end
hold off;

% Wavelet gain, Partial 2
figure(4);
plot(time,gain1_House,'k--','LineWidth',2); hold on;
plot(time,gain2_House,'r--','LineWidth',2); hold on;
plot(time,gain3_House,'k','LineWidth',2); axis tight;
% recessionplot;
legend('(6,32) quarters', '(6,80) quarters',...
    '(32,80) quarters', 'Location','best');
set(gca, 'box','off');
xtickangle(45);
set(gca, 'Fontname', 'times','FontSize',20);
switch savefigs
    case 1
        print('-depsc2', '-loose', 'figures/figure_11_right.eps');
end
hold off;





