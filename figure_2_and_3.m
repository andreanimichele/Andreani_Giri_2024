% Plot and export figures 2 and 3 (Wavelet coherences)

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
savefigs = 1;            % Scalar: 1 = save figures in .eps

addpath('data_input','figures');
% Add to path Wavelet toolboxes:
% 1 - M. Joana Soares and L. Aguiar-Conraria - ASToolbox
% 2 - Aslak Grinsted - Cross Wavelet and Wavelet Coherence Toolbox
addpath('wavelet-coherence-master','ASToolbox2018_Functions');

data = xlsread('data_input/dataset.xlsx');

t = data(:,1);                      % Time (1963:1 - 2019:3)
% Calculate quarterly growth rates
gdp = diff(log(data(:,2)));             % REAL GDP (Fred data)
mf_mort = diff(log(data(:,3)));        % Multifamily Residences
h_price = diff(log(data(:,4)));        % Average Sales Price of Houses, Dollars

X = [gdp mf_mort h_price];
% Years to add as xticklabels
years = 1965:5:2015;
years = years';
% Find indices of integer years in the t vector, for xticks in figures.
indices = find(ismember(t, years));
year_strg = string(years);

% GDP vs. Mortgages
figure(1);
f1 = subplot(1,1,1);
wtc(gdp,mf_mort,'Pad',1,'MaxScale',80);
xticks(indices);
xticklabels(year_strg);
colormap(jet);
set(f1,'FontName','times','FontSize',16);
hold on;
xtickangle(45);
switch savefigs
    case 1
        print('-depsc2', '-loose', 'figures/figure_2.eps');
end
hold off;

% GDP vs. House Prices
figure(2);
f2 = subplot(1,1,1);
wtc(gdp,h_price,'Pad',1,'MaxScale',80);
xticks(indices);
xticklabels(year_strg);
colormap(jet);
set(f2,'FontName','times','FontSize',16);
hold on;
xtickangle(45);
switch savefigs
    case 1
        print('-depsc2', '-loose', 'figures/figure_3.eps');
end
hold off;

