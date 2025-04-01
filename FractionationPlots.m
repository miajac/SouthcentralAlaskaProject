%% Fractionation Plots
clear all, close all, clc %use USGSGaugeDataset
[sFile, sPath] = uigetfile('*.xlsx', 'Select Database File');
sFullPath = fullfile(sPath, sFile);

tData = readtable(sFullPath);
%% 
vColumnLabels = tData.Properties.VariableNames(48:end);
vColumnLabelsArray = cellstr(vColumnLabels);
mFullData = table2array(tData(1:end, 48:end)); %first value is the number of samples (each as a row);second value is the numeric values that will be plotted, must be in an order that can be iterated through (DONT INCLUDE STRINGS)
vSampleLocations = string(table2cell(tData(:, 3))); %labels of each sample (KR4, MR3, etc)
%% Create datasets based on site locations 
mKR3Data = mFullData(vSampleLocations == "KR3", :);
mMR4Data = mFullData(vSampleLocations == "MR4", :);
mLS2Data = mFullData(vSampleLocations == "LS2", :);
mMCData = mFullData(vSampleLocations == "MC", :);

%% Define colors, markers, and labels
siteColors = {[0.35,0.70,0.90], [0.95,0.90,0.25], [0.90,0.60,0], [0,0.60,0.50]};
siteMarkers = {'o', 'o', 'o', 'o'};
siteLabels = {'K3', 'M4', 'MC', 'LS2'};
siteData = {mKR3Data, mMR4Data, mMCData, mLS2Data};

%% Plot the data
figure; % 1: Mg/Ca 2:Sr/Ca 3:Ba/Ca 4:U/Ca 5: Si/Ca 6: SO4/Ca 7: SO4/Si
hold on;
for i = 1:length(siteData)
    scatter(siteData{i}(:,1), siteData{i}(:,5), 50, siteColors{i}, 'o', 'filled', ...
    'MarkerEdgeColor', 'k', 'LineWidth', 0.5, 'DisplayName', siteLabels{i});
end
hold off;

xLabelFormatted = strrep(vColumnLabelsArray{1}, '_', '/');
yLabelFormatted = strrep(vColumnLabelsArray{5}, '_', '/');

% Set axis labels with formatted text
xlabel(xLabelFormatted, 'Interpreter', 'none'); 
ylabel(yLabelFormatted, 'Interpreter', 'none');
title('Fractionation Plot');
legend('Location', 'best');
grid on;
folderName = 'U:/GoA plots/NewPlots';
formatFileName = "Si-Ca,Mg-Ca.jpg"; 
fileName = sprintf(formatFileName);
fullFilePath = fullfile(folderName, fileName);
saveas(gcf, fullFilePath, 'jpg');
