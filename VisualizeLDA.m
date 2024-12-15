%% LDA plot in matlab
clear all, close all, clc
[sFile, sPath] = uigetfile('*.xlsx', 'Select Database File'); %use PCA
sFullPath = fullfile(sPath, sFile);
tData = readtable(sFullPath);
%%
vSampleLocations = string(table2cell(tData(1:end, 7))); %sample locations set
mFullData = table2array(tData(1:end, 1:6)); %data for each column
vLD1 = mFullData(:, 1); %LD1
vLD2 = mFullData(:,2); %LD2
vLD3 = mFullData(:,3); %LD3
vLD4 = mFullData(:,4); %LD4
vLD5 = mFullData(:,5); %LD5
vLD6 = mFullData(:,6); %LD6
%%
mKnikData = mFullData(vSampleLocations == "Knik", :);
mMatData = mFullData(vSampleLocations == "Matanuska", :);
mLSData = mFullData(vSampleLocations == "LS", :);
mMooseData = mFullData(vSampleLocations == "Moose", :);
mCastnerData = mFullData(vSampleLocations == "Castner", :);
mCanwellData = mFullData(vSampleLocations == "Canwell", :);
mGulkanaData = mFullData(vSampleLocations == "Gulkana", :);
%% With different symbols
figure;
plot(mKnikData(:,1), mKnikData(:,2), '^', 'DisplayName', 'Knik','Color','k','MarkerFaceColor',[0.35,0.70,0.90]);
hold on;
plot(mMatData(:,1), mMatData(:,2), 's', 'DisplayName', 'Matanuska','Color','k','MarkerFaceColor',[0.95,0.90,0.25]);
plot(mLSData(:,1), mLSData(:,2), 'd', 'DisplayName', 'LS','Color','k','MarkerFaceColor',[0,0.60,0.50]);
plot(mMooseData(:,1), mMooseData(:,2),'>', 'DisplayName', 'Moose','Color','k','MarkerFaceColor',[0.90,0.60,0]);
plot(mCastnerData(:,1), mCastnerData(:,2),'v', 'DisplayName', 'Castner','Color','k','MarkerFaceColor',[0,0.45,0.70]);
plot(mCanwellData(:,1), mCanwellData(:,2),'pentagram', 'DisplayName', 'Canwell','Color','k','MarkerFaceColor',[0.80,0.60,0.70]);
plot(mGulkanaData(:,1), mGulkanaData(:,2),'o', 'DisplayName', 'Gulkana','Color','k','MarkerFaceColor',[0.80,0.40,0]);
xlabel('LD1');  % Label for x-axis
ylabel('LD2');  % Label for y-axis
%title('LDA without isotopes');
legend('Location', 'SouthEast');
grid off;
hold off;
formatFileName = 'LDA_allwatersheds_without_isotopes_symbols.svg';
folderName = 'U:/GoA plots/NewPlots';
fileName = sprintf(formatFileName);
fullFilePath = fullfile(folderName, fileName);
saveas(gcf, fullFilePath, 'svg');
%% With same symbols
figure;
plot(mKnikData(:,1), mKnikData(:,2), 'o', 'DisplayName', 'Knik','Color','k','MarkerFaceColor',[0.35,0.70,0.90]);
hold on;
plot(mMatData(:,1), mMatData(:,2), 'o', 'DisplayName', 'Matanuska','Color','k','MarkerFaceColor',[0.95,0.90,0.25]);
plot(mLSData(:,1), mLSData(:,2), 'o', 'DisplayName', 'LS','Color','k','MarkerFaceColor',[0,0.60,0.50]);
plot(mMooseData(:,1), mMooseData(:,2),'o', 'DisplayName', 'Moose','Color','k','MarkerFaceColor',[0.90,0.60,0]);
plot(mCastnerData(:,1), mCastnerData(:,2),'o', 'DisplayName', 'Castner','Color','k','MarkerFaceColor',[0,0.45,0.70]);
plot(mCanwellData(:,1), mCanwellData(:,2),'o', 'DisplayName', 'Canwell','Color','k','MarkerFaceColor',[0.80,0.60,0.70]);
plot(mGulkanaData(:,1), mGulkanaData(:,2),'o', 'DisplayName', 'Gulkana','Color','k','MarkerFaceColor',[0.80,0.40,0]);
xlabel('LD1');  % Label for x-axis
ylabel('LD2');  % Label for y-axis
title('LDA without isotopes');
legend('Location', 'eastoutside');
grid off;
hold off;
formatFileName = 'LDA_allwatersheds_without_isotopes.jpg';
folderName = 'U:/GoA plots/NewPlots';
fileName = sprintf(formatFileName);
fullFilePath = fullfile(folderName, fileName);
saveas(gcf, fullFilePath, 'jpg');