%% Discharge plots
clear all, close all, clc %use DailyDischargeUSGS
[sFile, sPath] = uigetfile('*.xlsx', 'Select Database File');
sFullPath = fullfile(sPath, sFile);

tData = readtable(sFullPath);
%%
vMatQ = table2array(tData(:, 4)); 
vMatDates = datetime(convertStringsToChars(string(table2cell(tData(:, 2))))); 
vKnikQ = table2array(tData(:, 8)); 
vKnikDates = datetime(convertStringsToChars(string(table2cell(tData(:, 6))))); 
vLSQ = table2array(tData(:, 16)); 
vLSDates = datetime(convertStringsToChars(string(table2cell(tData(:, 14))))); 
vMCQ = table2array(tData(:, 12)); 
vMCDates = datetime(convertStringsToChars(string(table2cell(tData(:, 10))))); 

%% Plot - log axis 
figure;
plot(vKnikDates, vKnikQ, '-', 'DisplayName', 'K3','Color',[0.35,0.70,0.90]);%'Color','k','MarkerFaceColor',[0.35,0.70,0.90]
hold on;
plot(vMatDates, vMatQ, '-', 'DisplayName', 'M4','Color',[0.95,0.90,0.25]);
plot(vLSDates, vLSQ, '-', 'DisplayName', 'LS2','Color',[0,0.60,0.50]);
plot(vMCDates, vMCQ,'-', 'DisplayName', 'MC','Color',[0.90,0.60,0]);
%xlabel('Time');  % Label for x-axis
ylabel('Discharge (m^{3}/s)');
%title('LDA without isotopes');
legend('Location', 'North');
set(gca, 'YScale', 'log');  % Set y-axis to logarithmic scale
grid off;
hold off;
formatFileName = 'DailyDischargeUSGS.svg';
folderName = 'U:/GoA plots/NewPlots';
fileName = sprintf(formatFileName);
fullFilePath = fullfile(folderName, fileName);
saveas(gcf, fullFilePath, 'svg');