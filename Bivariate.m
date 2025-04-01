%% Creates bivariate plots for chemistry data
clear all, close all, clc %use AlaskaPCA2022_2023, not EM
[sFile, sPath] = uigetfile('*.xlsx', 'Select Database File');
sFullPath = fullfile(sPath, sFile);

tData = readtable(sFullPath);
%% 
vColumnLabels = tData.Properties.VariableNames(5:end);
vColumnLabelsArray = cellstr(vColumnLabels);
mFullData = table2array(tData(1:end, 5:end)); %first value is the number of samples (each as a row);second value is the numeric values that will be plotted, must be in an order that can be iterated through (DONT INCLUDE STRINGS)

vSampleDates = datetime(convertStringsToChars(string(table2cell(tData(:, 4))))); % just date not datetime
vSampleLocations = string(table2cell(tData(:, 1))); %labels of each sample (KR4, MR3, etc)
vWatershed = string(table2cell(tData(:,3))); % column that could be used to discriminate between watersheds each location is contained in 
%% Define watersheds
mMatData= mFullData(vWatershed == "Matanuska",:);
mMatDates = vSampleDates(vWatershed == "Matanuska");
mMatSites = vSampleLocations(vWatershed == "Matanuska");

mMCData = mFullData(vWatershed == "Moose",:);
mMCDates = vSampleDates(vWatershed == "Moose");
mMCSites = vSampleLocations(vWatershed == "Moose");

mKnikData= mFullData(vWatershed == "Knik",:);
mKnikDates = vSampleDates(vWatershed == "Knik");
mKnikSites = vSampleLocations(vWatershed == "Knik");

mLSData= mFullData(vWatershed == "LS",:);
mLSDates = vSampleDates(vWatershed == "LS");
mLSSites = vSampleLocations(vWatershed == "LS");

mGulkanaData= mFullData(vWatershed == "Gulkana",:);
mGulkanaDates = vSampleDates(vWatershed == "Gulkana");
mGulkanaSites = vSampleLocations(vWatershed == "Gulkana");

mCanwellData = mFullData(vWatershed == "Canwell",:);
mCanwellDates = vSampleDates(vWatershed == "Canwell");
mCanwellSites = vSampleLocations(vWatershed == "Canwell");

mCastnerData= mFullData(vWatershed == "Castner",:);
mCastnerDates = vSampleDates(vWatershed == "Castner");
mCastnerSites = vSampleLocations(vWatershed == "Castner");
%% Identify column number
vElements = tData.Properties.VariableNames;
sPrompt = ('Input an element for a graph');
sTitle = ('Element Input');
sInput =  inputdlg(sPrompt, sTitle); 
dPosition = find(ismember(vColumnLabels, sInput)); % Check if sInput is in the dataset and find its position
if isempty(dPosition)
    fprintf('The input element is not in the dataset.\n');
else
    fprintf('Input element exists in dataset at column number: %d\n', dPosition);
end
%% Bivariate plots
figure('Position', [100, 100, 1000, 800]);
m=8;
t = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
% First subplot (Ca vs Rb)
nexttile(1);   
plot(mKnikData(:,13), mKnikData(:,23), '^', 'DisplayName', 'Knik', 'Color', 'k', 'MarkerFaceColor', [0.35,0.70,0.90],'MarkerSize',m);
hold on;
plot(mMatData(:,13), mMatData(:,23), 's', 'DisplayName', 'Matanuska', 'Color', 'k', 'MarkerFaceColor', [0.95,0.90,0.25],'MarkerSize',m);
plot(mLSData(:,13), mLSData(:,23), 'd', 'DisplayName', 'LS', 'Color', 'k', 'MarkerFaceColor', [0,0.60,0.50],'MarkerSize',m);
plot(mMCData(:,13), mMCData(:,23), '>', 'DisplayName', 'Moose', 'Color', 'k', 'MarkerFaceColor', [0.90,0.60,0],'MarkerSize',m);
plot(mCastnerData(:,13), mCastnerData(:,23), 'v', 'DisplayName', 'Castner', 'Color', 'k', 'MarkerFaceColor', [0,0.45,0.70],'MarkerSize',m);
plot(mCanwellData(:,13), mCanwellData(:,23), 'pentagram', 'DisplayName', 'Canwell', 'Color', 'k', 'MarkerFaceColor', [0.80,0.60,0.70],'MarkerSize',m);
plot(mGulkanaData(:,13), mGulkanaData(:,23), 'o', 'DisplayName', 'Gulkana', 'Color', 'k', 'MarkerFaceColor', [0.80,0.40,0],'MarkerSize',m);
xlabel('Ca (mg/L)');  % Label for x-axis
ylabel('Rb (mg/L)');  % Label for y-axis
set(gca, 'FontSize', 12);
%legend('Location', 'NorthEast');
grid off;
hold off;

% Second subplot (Se vs U)
nexttile(2);
plot(mKnikData(:,48), mKnikData(:,27), '^', 'DisplayName', 'Knik', 'Color', 'k', 'MarkerFaceColor', [0.35,0.70,0.90],'MarkerSize',m);
hold on;
plot(mMatData(:,48), mMatData(:,27), 's', 'DisplayName', 'Matanuska', 'Color', 'k', 'MarkerFaceColor', [0.95,0.90,0.25],'MarkerSize',m);
plot(mLSData(:,48), mLSData(:,27), 'd', 'DisplayName', 'LS', 'Color', 'k', 'MarkerFaceColor', [0,0.60,0.50],'MarkerSize',m);
plot(mMCData(:,48), mMCData(:,27), '>', 'DisplayName', 'Moose', 'Color', 'k', 'MarkerFaceColor', [0.90,0.60,0],'MarkerSize',m);
plot(mCastnerData(:,48), mCastnerData(:,27), 'v', 'DisplayName', 'Castner', 'Color', 'k', 'MarkerFaceColor', [0,0.45,0.70],'MarkerSize',m);
plot(mCanwellData(:,48), mCanwellData(:,27), 'pentagram', 'DisplayName', 'Canwell', 'Color', 'k', 'MarkerFaceColor', [0.80,0.60,0.70],'MarkerSize',m);
plot(mGulkanaData(:,48), mGulkanaData(:,27), 'o', 'DisplayName', 'Gulkana', 'Color', 'k', 'MarkerFaceColor', [0.80,0.40,0],'MarkerSize',m);
xlabel('U (mg/L)');  % Label for x-axis
ylabel('Se (mg/L)');  % Label for y-axis
set(gca, 'FontSize', 12);
legend('Location', 'NorthEast','FontSize',12);
grid off;
hold off;

% Third subplot (Cr vs SO_4)
nexttile(3);
plot(mKnikData(:,16), mKnikData(:,4), '^', 'DisplayName', 'Knik', 'Color', 'k', 'MarkerFaceColor', [0.35,0.70,0.90],'MarkerSize',m);
hold on;
plot(mMatData(:,16), mMatData(:,4), 's', 'DisplayName', 'Matanuska', 'Color', 'k', 'MarkerFaceColor', [0.95,0.90,0.25],'MarkerSize',m);
plot(mLSData(:,16), mLSData(:,4), 'd', 'DisplayName', 'LS', 'Color', 'k', 'MarkerFaceColor', [0,0.60,0.50],'MarkerSize',m);
plot(mMCData(:,16), mMCData(:,4), '>', 'DisplayName', 'Moose', 'Color', 'k', 'MarkerFaceColor', [0.90,0.60,0],'MarkerSize',m);
plot(mCastnerData(:,16), mCastnerData(:,4), 'v', 'DisplayName', 'Castner', 'Color', 'k', 'MarkerFaceColor', [0,0.45,0.70],'MarkerSize',m);
plot(mCanwellData(:,16), mCanwellData(:,4), 'pentagram', 'DisplayName', 'Canwell', 'Color', 'k', 'MarkerFaceColor', [0.80,0.60,0.70],'MarkerSize',m);
plot(mGulkanaData(:,16), mGulkanaData(:,4), 'o', 'DisplayName', 'Gulkana', 'Color', 'k', 'MarkerFaceColor', [0.80,0.40,0],'MarkerSize',m);
xlabel('Cr (mg/L)');  % Label for x-axis
ylabel('SO_4 (mg/L)');  % Label for y-axis (with subscript for 4)
set(gca, 'FontSize', 12);
%legend('Location', 'NorthEast');
grid off;
hold off;

% Fourth subplot (Sr vs U)
nexttile(4);
plot(mKnikData(:,48), mKnikData(:,35), '^', 'DisplayName', 'Knik', 'Color', 'k', 'MarkerFaceColor', [0.35,0.70,0.90],'MarkerSize',m);
hold on;
plot(mMatData(:,48), mMatData(:,35), 's', 'DisplayName', 'Matanuska', 'Color', 'k', 'MarkerFaceColor', [0.95,0.90,0.25],'MarkerSize',m);
plot(mLSData(:,48), mLSData(:,35), 'd', 'DisplayName', 'LS', 'Color', 'k', 'MarkerFaceColor', [0,0.60,0.50],'MarkerSize',m);
plot(mMCData(:,53), mMCData(:,35), '>', 'DisplayName', 'Moose', 'Color', 'k', 'MarkerFaceColor', [0.90,0.60,0],'MarkerSize',m);
plot(mCastnerData(:,53), mCastnerData(:,35), 'v', 'DisplayName', 'Castner', 'Color', 'k', 'MarkerFaceColor', [0,0.45,0.70],'MarkerSize',m);
plot(mCanwellData(:,53), mCanwellData(:,35), 'pentagram', 'DisplayName', 'Canwell', 'Color', 'k', 'MarkerFaceColor', [0.80,0.60,0.70],'MarkerSize',m);
plot(mGulkanaData(:,53), mGulkanaData(:,35), 'o', 'DisplayName', 'Gulkana', 'Color', 'k', 'MarkerFaceColor', [0.80,0.40,0],'MarkerSize',m);
xlabel('U (mg/L)');  % Label for x-axis
ylabel('Cs (mg/L)');  % Label for y-axis
set(gca, 'FontSize', 12);
%legend('Location', 'NorthEast');
grid off;
hold off;

% Save the figure
formatFileName = 'Bivariate_Plots_Combined.svg';  
folderName = 'U:/GoA plots/NewPlots';        
fileName = sprintf(formatFileName);             
fullFilePath = fullfile(folderName, fileName);   
saveas(gcf, fullFilePath, 'svg');                
