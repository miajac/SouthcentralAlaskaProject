%% Discharge plots
clear all, close all, clc %use Calcium
[sFile, sPath] = uigetfile('*.xlsx', 'Select Database File');
sFullPath = fullfile(sPath, sFile);

tData = readtable(sFullPath);
%%
vResults = table2array(tData(:, 20)); 
vSampleLocations = string(table2cell(tData(:, 2)));
vDates = datetime(convertStringsToChars(string(table2cell(tData(:, 3))))); 
%%
vMatData= vResults(vSampleLocations == "15284000",:);
vMatDates = vDates(vSampleLocations == "15284000");

vKnikData= vResults(vSampleLocations == "15281000",:);
vKnikDates = vDates(vSampleLocations == "15281000");

vLSData= vResults(vSampleLocations == "15290000",:);
vLSDates = vDates(vSampleLocations == "15290000");

vMCData= vResults(vSampleLocations == "15283700",:);
vMCDates = vDates(vSampleLocations == "15283700");
%% Plot - log axis 
figure('Position', [100, 100, 2500, 500]);
plot(vKnikDates, vKnikData, '-o', 'DisplayName', 'K3','Color',[0.35,0.70,0.90],'MarkerFaceColor', [0.35,0.70,0.90],'MarkerEdgeColor','k');%'Color','k','MarkerFaceColor',[0.35,0.70,0.90]
hold on;
plot(vMatDates, vMatData, '-o', 'DisplayName', 'M4','Color',[0.95,0.90,0.25],'MarkerFaceColor',[0.95,0.90,0.25],'MarkerEdgeColor','k');
plot(vLSDates, vLSData, '-o', 'DisplayName', 'LS2','Color',[0,0.60,0.50],'MarkerFaceColor',[0,0.60,0.50],'MarkerEdgeColor','k');
plot(vMCDates, vMCData,'-o', 'DisplayName', 'MC','Color',[0.90,0.60,0],'MarkerFaceColor',[0.90,0.60,0],'MarkerEdgeColor','k');
%xlabel('Time');  % Label for x-axis
ylabel('Ca (mg/L)');
%title('LDA without isotopes');
legend('Location', 'NorthEast');
%set(gca, 'YScale', 'log');  % Set y-axis to logarithmic scale
grid off;
hold off;
formatFileName = 'HistoricalUSGS_Ca.svg';
folderName = 'U:/GoA plots/NewPlots';
fileName = sprintf(formatFileName);
fullFilePath = fullfile(folderName, fileName);
saveas(gcf, fullFilePath, 'svg');