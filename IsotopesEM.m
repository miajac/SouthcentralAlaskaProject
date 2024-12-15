%% Plots
clear all, close all, clc %use GOAendmemberJuly
[sFile, sPath] = uigetfile('*.xlsx', 'Select Database File');
sFullPath = fullfile(sPath, sFile);

tData = readtable(sFullPath);
%% 
fullDataset = table2array(tData(1:end, 6:7)); 
vSampleID = string(table2cell(tData(:, 2))); %labels of each sample (KR4, MR3, etc)
vEndmember = string(table2cell(tData(:,4))); % column that defines types of endmembers 
vWatershed = string(table2cell(tData(:,3)));
%% Watershed breakdown
mKnikFull = fullDataset(vWatershed == "Knik",:);
vK_EM = vEndmember(vWatershed == "Knik");
vK_SiteID = vSampleID(vWatershed == "Knik");
mMatFull = fullDataset(vWatershed == "Matanuska",:);
vM_EM = vEndmember(vWatershed == "Matanuska");
vM_SiteID = vSampleID(vWatershed == "Matanuska");
mLSFull = fullDataset(vWatershed == "Little Susitna",:);
vLS_EM = vEndmember(vWatershed == "Little Susitna");
vLS_SiteID = vSampleID(vWatershed == "Little Susitna");
mCastnerFull = fullDataset(vWatershed == "Castner",:);
vCT_EM = vEndmember(vWatershed == "Castner");
vCT_SiteID = vSampleID(vWatershed == "Castner");
mCanwellFull = fullDataset(vWatershed == "Canwell",:);
vCW_EM = vEndmember(vWatershed == "Canwell");
vCW_SiteID = vSampleID(vWatershed == "Canwell");
mGulkanaFull = fullDataset(vWatershed == "Gulkana",:);
vG_EM = vEndmember(vWatershed == "Gulkana");
vG_SiteID = vSampleID(vWatershed == "Gulkana");

DeltaWatersheds = ["Castner", "Canwell", "Gulkana"];
index = ismember(vWatershed, DeltaWatersheds);
mDeltasFull = fullDataset(index, :);
index_EM = ismember(vWatershed,DeltaWatersheds);
vDelta_EM = vEndmember(index_EM,:);
index_SiteID = ismember(vWatershed,DeltaWatersheds);
vDelta_SiteID = vSampleID(index_SiteID,:);

%%
KR2 = mKnikFull(vK_SiteID == "KR2",:);
KR3 = mKnikFull(vK_SiteID == "KR3",:);
KR4 = mKnikFull(vK_SiteID == "KR4",:);
HC = mKnikFull(vK_SiteID == "HC",:);
GC = mKnikFull(vK_SiteID == "Goat Creek",:);
K_MR4 = mKnikFull(vK_SiteID == "Mat",:);
K_Spring = mKnikFull(vK_EM == "Spring",:);

MR1 = mMatFull(vM_SiteID == "MR1",:);
MR2 = mMatFull(vM_SiteID == "MR2",:);
MR3 = mMatFull(vM_SiteID == "MR3",:);
MR4 = mMatFull(vM_SiteID == "MR4",:);
MR5 = mMatFull(vM_SiteID == "MR5",:);
MC = mMatFull(vM_SiteID == "MC",:);
CC = mMatFull(vM_SiteID == "CC",:);
Hick = mMatFull(vM_SiteID == "Hick's Creek",:);
Pino = mMatFull(vM_SiteID == "Pinocle Creek",:);
Pur = mMatFull(vM_SiteID == "Purinton Creek",:);
Chick = mMatFull(vM_SiteID == "Chickaloon River",:);
KR = mMatFull(vM_SiteID == "Kings River",:);
YJ = mMatFull(vM_SiteID == "Yellow Jacket Creek",:);
UM = mMatFull(vM_SiteID == "Upper Matanuska River",:);
M_Spring = mMatFull(vM_EM == "Spring",:);

LS1 = mLSFull(vLS_SiteID == "LS1",:);
LS15 = mLSFull(vLS_SiteID == "LS1.5",:);
LS2 = mLSFull(vLS_SiteID == "LS2",:);
LS3 = mLSFull(vLS_SiteID == "LS3",:);
LS4 = mLSFull(vLS_SiteID == "LS4",:);
AR = mLSFull(vLS_SiteID == "Archangel River",:);
IL = mLSFull(vLS_SiteID == "Ivory Lake",:);
MGS = mLSFull(vLS_SiteID == "Mint Glacier Subglacial",:);
SS = mLSFull(vLS_SiteID == "Snowbird Supraglacial",:);
SP = mLSFull(vLS_SiteID == "Snowbird Periglacial",:);
FC = mLSFull(vLS_SiteID == "Fishhook Creek",:);
LS_Spring = mLSFull(vLS_EM == "Spring",:);

CT1 = mDeltasFull(vDelta_SiteID == "CT1",:);
CT2 = mDeltasFull(vDelta_SiteID == "CT2",:);
CS = mDeltasFull(vDelta_SiteID == "Castner subglacial",:);
G1 = mDeltasFull(vDelta_SiteID == "G1",:);
G2 = mDeltasFull(vDelta_SiteID == "G2",:);
GS = mDeltasFull(vDelta_SiteID == "Gulkana Supraglacial 1",:);
CoC = mDeltasFull(vDelta_SiteID == "College Creek",:);
CW1 = mDeltasFull(vDelta_SiteID == "CW1",:);
CW2 = mDeltasFull(vDelta_SiteID == "CW2",:);

%% combined endmembers
figure;
legendFontSize = 10;
markerSize = 75;
markerSizeD = 40;
markerSizeS = 75;
markerSizeG = 45;
l=2.5; %line width
position1 = [0.07, 0.57, 0.4, 0.4];
position2 = [0.55, 0.57, 0.4, 0.4];
position3 = [0.07, 0.1, 0.4, 0.4];
position4 = [0.55, 0.1, 0.4, 0.4];
subplot('Position', position1); % Knik plot
% Define the range for DO (δ¹⁸O) based on your x-axis range
DO_rangeK = linspace(-24,-19, 100); 
% Compute DD (δ²H) using the meteoric water line equation
DD_meteoricK = 8 * DO_rangeK + 10;
% Plot the meteoric water line
plot(DO_rangeK, DD_meteoricK, 'k--', 'DisplayName', 'GMWL', 'LineWidth', 0.5, 'HandleVisibility', 'off');
hold on
scatter(KR2(:,1), KR2(:,2), 'ko','MarkerFaceColor',[0.992, 0.851, 0.627],'DisplayName', 'KR2','SizeData',markerSize, 'LineWidth', 0.5, 'HandleVisibility', 'on');
scatter(KR3(:,1), KR3(:,2), 'ko','MarkerFaceColor',[0.961, 0.498, 0.090],'DisplayName', 'KR3','SizeData',markerSize, 'LineWidth',0.5, 'HandleVisibility', 'on');
scatter(KR4(:,1), KR4(:,2), 'ko','MarkerFaceColor',[0.651, 0.212, 0.012],'DisplayName', 'KR4','SizeData',markerSize, 'LineWidth', 0.5, 'HandleVisibility', 'on');
scatter(K_Spring(:,1),K_Spring(:,2),'kd','MarkerFaceColor',[0, 0.620, 0.451],'DisplayName', 'Springs','SizeData',markerSizeD);
scatter(HC(:,1), HC(:,2), 'ks','MarkerFaceColor',[1.000, 0.843, 0.000],'DisplayName', 'Tributaries','SizeData',markerSizeS);
scatter(GC(:,1), GC(:,2), 'ks','MarkerFaceColor',[1.000, 0.843, 0.000],'DisplayName', 'Goat Creek','SizeData',markerSizeS, 'HandleVisibility', 'off'); 
scatter(K_MR4(:,1),K_MR4(:,2),'ks','MarkerFaceColor', [1.000, 0.843, 0.000],'DisplayName', 'Matanuska','SizeData',markerSizeS, 'HandleVisibility', 'off');

ylabel('δ²H (‰)'); 
ax = gca;
ax.FontSize = 10; 
title('Knik');
legend('Location', 'SouthEast','FontSize',legendFontSize);
grid off;
hold off;

subplot('Position', position2); % Matanuska plot
markerSize1=markerSize;
DO_rangeM = linspace(-25,-19, 100); 
DD_meteoricM = 8 * DO_rangeM + 10;
plot(DO_rangeM, DD_meteoricM, 'k--', 'DisplayName', 'GMWL', 'LineWidth', 0.5, 'HandleVisibility', 'off');
hold on
scatter(MR1(:, 1), MR1(:, 2), 'ko','MarkerFaceColor',[0.992, 0.851, 0.627],'DisplayName', 'MR1','SizeData',markerSize, 'LineWidth', 0.5, 'HandleVisibility', 'on'); 
scatter(MR2(:,1), MR2(:,2), 'ko','MarkerFaceColor',[0.992, 0.722, 0.388],'DisplayName', 'MR2','SizeData',markerSize, 'LineWidth', 0.5, 'HandleVisibility', 'on');
scatter(MR3(:,1), MR3(:,2), 'ko','MarkerFaceColor',[0.961, 0.498, 0.090],'DisplayName', 'MR3','SizeData',markerSize, 'LineWidth', 0.5, 'HandleVisibility', 'on');
scatter(MR4(:,1), MR4(:,2), 'ko','MarkerFaceColor',[0.902, 0.333, 0.051],'DisplayName', 'MR4','SizeData',markerSize, 'LineWidth', 0.5, 'HandleVisibility', 'on');
scatter(MR5(:,1), MR5(:,2), 'ko','MarkerFaceColor',[0.651, 0.212, 0.012],'DisplayName', 'MR5','SizeData',markerSize, 'LineWidth', 0.5, 'HandleVisibility', 'on');
scatter(M_Spring(:,1),M_Spring(:,2),'kd','MarkerFaceColor',[0, 0.620, 0.451],'DisplayName', 'Springs','SizeData',markerSizeD);
scatter(YJ(:,1), YJ(:,2), 'ks','MarkerFaceColor',[1.000, 0.843, 0.000],'DisplayName', 'Tributaries','SizeData',markerSizeS);
scatter(CC(:,1), CC(:,2), 'ks','MarkerFaceColor',[1.000, 0.843, 0.000],'DisplayName', 'Caribou Creek','SizeData',markerSizeS, 'HandleVisibility', 'off');
scatter(UM(:,1), UM(:,2), 'ks','MarkerFaceColor',[1.000, 0.843, 0.000],'DisplayName', 'Upper Matanuska River','SizeData',markerSizeS, 'HandleVisibility', 'off');
scatter(Pino(:,1), Pino(:,2), 'ks','MarkerFaceColor',[1.000, 0.843, 0.000],'DisplayName', 'Pinochle Creek','SizeData',markerSizeS, 'HandleVisibility', 'off');
scatter(Hick(:,1), Hick(:,2), 'ks','MarkerFaceColor',[1.000, 0.843, 0.000],'DisplayName', 'Hicks Creek','SizeData',markerSizeS, 'HandleVisibility', 'off');
scatter(Pur(:,1), Pur(:,2), 'ks','MarkerFaceColor',[1.000, 0.843, 0.000],'DisplayName', 'Purinton Creek','SizeData',markerSizeS, 'HandleVisibility', 'off');
scatter(Chick(:,1), Chick(:,2), 'ks','MarkerFaceColor',[1.000, 0.843, 0.000],'DisplayName', 'Chickaloon River','SizeData',markerSizeS, 'HandleVisibility', 'off');
scatter(KR(:,1), KR(:,2), 'ks','MarkerFaceColor',[1.000, 0.843, 0.000],'DisplayName', 'Kings River','SizeData',markerSizeS, 'HandleVisibility', 'off');
scatter(MC(:,1), MC(:,2), 'ks','MarkerFaceColor',[1.000, 0.843, 0.000],'DisplayName', 'Moose Creek','SizeData',markerSizeS, 'HandleVisibility', 'off');

ax = gca;
ax.FontSize = 10; 
title('Matanuska');
legend('Location', 'SouthEast','FontSize',legendFontSize);
grid off;
hold off;

subplot('Position', position3); %LS plot
DO_rangeLS = linspace(-21.5,-18, 100); 
DD_meteoricLS = 8 * DO_rangeLS + 10;
plot(DO_rangeLS, DD_meteoricLS, 'k--', 'DisplayName', 'GMWL', 'LineWidth', 0.5, 'HandleVisibility', 'off');
hold on
scatter(LS1(:, 1), LS1(:, 2), 'ko','MarkerFaceColor',[0.992, 0.851, 0.627],'DisplayName', 'LS1','SizeData',markerSize, 'LineWidth', 0.5, 'HandleVisibility', 'on'); 
scatter(LS15(:,1), LS15(:,2), 'ko','MarkerFaceColor',[0.992, 0.722, 0.388],'DisplayName', 'LS1.5','SizeData',markerSize, 'LineWidth', 0.5, 'HandleVisibility', 'on');
scatter(LS2(:,1), LS2(:,2), 'ko','MarkerFaceColor',[0.961, 0.498, 0.090],'DisplayName', 'LS2','SizeData',markerSize, 'LineWidth', 0.5, 'HandleVisibility', 'on');
scatter(LS3(:,1), LS3(:,2), 'ko','MarkerFaceColor',[0.902, 0.333, 0.051],'DisplayName', 'LS3','SizeData',markerSize, 'LineWidth', 0.5, 'HandleVisibility', 'on');
scatter(LS4(:,1), LS4(:,2), 'ko','MarkerFaceColor',[0.651, 0.212, 0.012],'DisplayName', 'LS4','SizeData',markerSize, 'LineWidth', 0.5, 'HandleVisibility', 'on');
scatter(LS_Spring(:,1),LS_Spring(:,2),'kd','MarkerFaceColor',[0, 0.620, 0.451],'DisplayName', 'Springs','SizeData',markerSizeD);
scatter(AR(:,1), AR(:,2), 'ks','MarkerFaceColor',[1.000, 0.843, 0.000],'DisplayName', 'Tributaries','SizeData',markerSizeS);
scatter(FC(:,1), FC(:,2), 'ks','MarkerFaceColor',[1.000, 0.843, 0.000],'DisplayName', 'Fishhook Creek','SizeData',markerSizeS, 'HandleVisibility', 'off');
scatter(IL(:,1), IL(:,2), 'k^','MarkerFaceColor',[0.616, 0.765, 0.902],'DisplayName', 'Glacial','SizeData',markerSizeG); %Glacial Lake
scatter(SS(:,1), SS(:,2), 'k^','MarkerFaceColor',[0.616, 0.765, 0.902],'DisplayName','Supraglacial','SizeData',markerSizeG,'HandleVisibility','off');
scatter(MGS(:,1), MGS(:,2), 'k^','MarkerFaceColor',[0.616, 0.765, 0.902],'DisplayName', 'Subglacial','SizeData',markerSizeG,'HandleVisibility','off');
scatter(SP(:,1), SP(:,2), 'k^','MarkerFaceColor',[0.616, 0.765, 0.902],'DisplayName','Periglacial','SizeData',markerSizeG, 'HandleVisibility','off');

xlabel('δ¹⁸O (‰)');  
ylabel('δ²H (‰)');
ax = gca;
ax.FontSize = 10; 
title('Little Susitna');
legend('Location', 'SouthEast','FontSize',legendFontSize);
grid off;
hold off;

subplot('Position', position4); %Delta's Plot
DO_rangeD = linspace(-23,-20, 100); 
DD_meteoricD = 8 * DO_rangeD + 10;
plot(DO_rangeD, DD_meteoricD, 'k--', 'DisplayName', 'GMWL', 'LineWidth', 0.5, 'HandleVisibility', 'off');
hold on 
scatter(CT1(:, 1), CT1(:, 2), 'ko', 'MarkerFaceColor', [1.000, 0.949, 0.800], 'DisplayName', 'CT1', 'SizeData', markerSize, 'LineWidth', 0.5, 'HandleVisibility', 'on');
scatter(CT2(:, 1), CT2(:, 2), 'ko', 'MarkerFaceColor', [0.961, 0.871, 0.702],'DisplayName', 'CT2', 'SizeData', markerSize, 'LineWidth', 0.5, 'HandleVisibility', 'on');
scatter(G1(:, 1), G1(:, 2), 'ko', 'MarkerFaceColor', [0.992, 0.722, 0.388],'DisplayName', 'G1', 'SizeData', markerSize, 'LineWidth', 0.5, 'HandleVisibility', 'on');
scatter(G2(:, 1), G2(:, 2), 'ko', 'MarkerFaceColor', [0.961, 0.498, 0.090],'DisplayName', 'G2', 'SizeData', markerSize, 'LineWidth', 0.5, 'HandleVisibility', 'on');
scatter(CW1(:, 1), CW1(:, 2), 'ko', 'MarkerFaceColor', [0.800, 0.361, 0.400],'DisplayName', 'CW1', 'SizeData', markerSize, 'LineWidth', 0.5, 'HandleVisibility', 'on');
scatter(CW2(:, 1), CW2(:, 2), 'ko', 'MarkerFaceColor', [0.600, 0.200, 0.251],'DisplayName', 'CW2', 'SizeData', markerSize, 'LineWidth', 0.5, 'HandleVisibility', 'on');
scatter(CoC(:, 1), CoC(:, 2), 'ks', 'MarkerFaceColor', [1.000, 0.843, 0.000],'DisplayName', 'Tributary', 'SizeData', markerSizeS);
scatter(GS(:, 1), GS(:, 2), 'k^', 'MarkerFaceColor', [0.616, 0.765, 0.902], 'DisplayName', 'Glacial', 'SizeData', markerSizeG);
scatter(CS(:, 1), CS(:, 2), 'k^', 'MarkerFaceColor', [0.616, 0.765, 0.902],'DisplayName', 'Subglacial', 'SizeData', markerSizeG,'HandleVisibility','off');
xlabel('δ¹⁸O (‰)');  
ax = gca;
ax.FontSize = 10; 
title('Delta Range');
legend('Location', 'SouthEast','FontSize',legendFontSize);
grid off;
hold off;
set(gcf, 'Position', [100, 100, 1300, 900]);
% Save the combined figure
folderName = 'U:/GoA plots/NewPlots';
fileName = 'Endmember_isotopes_tribcombined.svg';
fullFilePath = fullfile(folderName, fileName);
saveas(gcf, fullFilePath, 'svg');