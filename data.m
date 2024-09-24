% Download data used for wavelet coherence and partial coherence
% estimation from Fred data. 

% Mortgages, House Prices, and Business Cycle Dynamics: A Medium-Run Exploration Using the Continuous Wavelet Transform
% International Review of Economics & Finance, 2024, Volume 94

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

folderName = 'data_input';  % Specify the folder name

% Check if the folder exists
if ~exist(folderName, 'dir')
    % If the folder does not exist, create it
    mkdir(folderName);
    fprintf('Folder "%s" created successfully.\n', folderName);
else
    fprintf('Folder "%s" already exists.\n', folderName);
end

% First observation: {'1963-01-01'}
t1 = datetime('1963-01-01');
t2 = datetime('2019-10-01');

% 
% FRED DATA
% 

% Real GDP data download
seriesID = 'GDPC1';
csvFile = fullfile(folderName, 'real_gdp.csv');   % Temporary CSV file
xlsxFile = fullfile(folderName, 'real_gdp.xlsx'); % Final Excel file
% Correct URL for downloading the CSV file
url = sprintf('https://fred.stlouisfed.org/graph/fredgraph.csv?id=%s', seriesID);
% Download the CSV file from the FRED website
websave(csvFile, url);
% Read the CSV file into a table
dataTable = readtable(csvFile);
gdp.t = datetime(dataTable.DATE, 'InputFormat', 'yyyy-MM-dd');
gdp.gdp = dataTable.GDPC1;
tb = find(gdp.t == t1);
te = find(gdp.t == t2);
gdp.t = gdp.t(tb:te,:);
gdp.gdp = gdp.gdp(tb:te,:);
% Save the table as an Excel file
writetable(dataTable, xlsxFile);
delete(csvFile);

% All Sectors; Multifamily Residential Mortgages; Asset, Level 
seriesID = 'ASMRMA';
csvFile = fullfile(folderName, 'mortgages.csv');   % Temporary CSV file
xlsxFile = fullfile(folderName, 'mortgages.xlsx'); % Final Excel file
% Correct URL for downloading the CSV file
url = sprintf('https://fred.stlouisfed.org/graph/fredgraph.csv?id=%s', seriesID);
% Download the CSV file from the FRED website
websave(csvFile, url);
% Read the CSV file into a table
dataTable = readtable(csvFile);
mort.t = datetime(dataTable.DATE, 'InputFormat', 'yyyy-MM-dd');
mort.mort = dataTable.ASMRMA;
tb = find(mort.t == t1);
te = find(mort.t == t2);
mort.t = mort.t(tb:te,:);
mort.mort = mort.mort(tb:te,:);
% Save the table as an Excel file
writetable(dataTable, xlsxFile);
delete(csvFile);

% Average Sales Price of Houses Sold for the United States
seriesID = 'ASPUS';
csvFile = fullfile(folderName, 'house_prices.csv');   % Temporary CSV file
xlsxFile = fullfile(folderName, 'house_prices.xlsx'); % Final Excel file
% Correct URL for downloading the CSV file
url = sprintf('https://fred.stlouisfed.org/graph/fredgraph.csv?id=%s', seriesID);
% Download the CSV file from the FRED website
websave(csvFile, url);
% Read the CSV file into a table
dataTable = readtable(csvFile);
housep.t = datetime(dataTable.DATE, 'InputFormat', 'yyyy-MM-dd');
housep.housep = dataTable.ASPUS;
tb = find(housep.t == t1);
te = find(housep.t == t2);
housep.t = housep.t(tb:te,:);
housep.housep = housep.housep(tb:te,:);
% Save the table as an Excel file
writetable(dataTable, xlsxFile);
delete(csvFile);

% Consumer Price Index for All Urban Consumers: All Items in U.S. City Average
seriesID = 'CPIAUCSL';
csvFile = fullfile(folderName, 'cpi.csv');   % Temporary CSV file
xlsxFile = fullfile(folderName, 'cpi.xlsx'); % Final Excel file
% Correct URL for downloading the CSV file
url = sprintf('https://fred.stlouisfed.org/graph/fredgraph.csv?id=%s', seriesID);
% Download the CSV file from the FRED website
websave(csvFile, url);
% Read the CSV file into a table
dataTable = readtable(csvFile);
cpi = table2timetable(dataTable, 'RowTimes', 'DATE');
cpiq = retime(cpi, 'quarterly', 'mean'); % Use 'mean' if you need average values
cpiq = timetable2table(cpiq);
cpiqq.t = datetime(cpiq.DATE, 'InputFormat', 'yyyy-MM-dd');
cpiqq.cpi = cpiq.CPIAUCSL;
tb = find(cpiqq.t == t1);
te = find(cpiqq.t == t2);
cpiqq.t = cpiqq.t(tb:te,:);
cpiqq.cpi = cpiqq.cpi(tb:te,:);
% Save the table as an Excel file
writetable(dataTable, xlsxFile);
delete(csvFile);

% 
% DATA WRANGLING: Define relevant time series
% 

r_gdp           = gdp.gdp;
cpi             = cpiqq.cpi;
mortgages       = mort.mort;
house_p         = housep.housep;

% 
% Create a numeric date vector for quarterly dates
% 
dates           = gdp.t;
% Extract year and quarter
yearArray       = year(dates);
quarterArray    = quarter(dates);
% Convert quarters to decimals
quarterDecimal  = (quarterArray - 1) * 0.25;
% Combine year with decimal quarter
T               = yearArray + quarterDecimal;

% Mortgages and house prices in real terms
r_mortgages     = mortgages./cpi;
r_house         = house_p./cpi;
% Define raw dataset
Y = [T r_gdp r_mortgages r_house];
% Export dataset
writematrix(Y, "data_input/dataset.xlsx");

fprintf('Dataset exported.\n', folderName);
