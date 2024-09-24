% Main script file.
% Download data, and plot main output figures 1-5 from the paper 
% "Mortgages, House Prices, and Business Cycle Dynamics: A Medium-Run
% Exploration Using the Continuous Wavelet Transform,"

% Michele Andreani
% Financial Stability Research Office
% Central Bank of Malta

% Federico Giri
% Universitá Politecnica delle Marche 
% Dipartimento di Scienze Economiche e Sociali (DiSES)

% Last edit: 24 September, 2024

clear all
close all
clc

% Download, wrangle and export data sheets for Wavelet coherence estimation
run("data.m");
% Plot and export figure 1 (Global wavelet power spectra)
run("figure_1.m");
% Plot and export figures 2 and 3 (Wavelet coherences)
run("figure_2_and_3.m");
% Plot and export subplots for figure 2
run("figure_2_subplot.m");
% Plot and export subplots for figure 3
run("figure_3_subplot.m");
% Plot and export figures 4 and 5 (Partial Wavelet coherences)
run("figure_4_5.m");
% Plot and export figures in appendix.
run("figure_appendix.m");

% 
% End of script.
% 