%% Endmember PCA
clear all, close all, clc
[sFile, sPath] = uigetfile('*.xlsx', 'Select Database File');
sFullPath = fullfile(sPath, sFile);
mC = readtable(sFullPath);
mA = mC(:, 6:end);
fullDataset = table2array(mA);
vWatershed = string(table2cell(mC(:, 3)));
vEndmember = string(table2cell(mC(:,4)));
vSiteID = string(table2cell(mC(:,2)));
%% Define Endmember Categories
mKnikFull = fullDataset(vWatershed == "Knik",:);
vK_EM = vEndmember(vWatershed == "Knik");
vK_SiteID = vSiteID(vWatershed == "Knik");
mMatFull = fullDataset(vWatershed == "Matanuska",:);
vM_EM = vEndmember(vWatershed == "Matanuska");
vM_SiteID = vSiteID(vWatershed == "Matanuska");
mLSFull = fullDataset(vWatershed == "Little Susitna",:);
vLS_EM = vEndmember(vWatershed == "Little Susitna");
vLS_SiteID = vSiteID(vWatershed == "Little Susitna");
mCastnerFull = fullDataset(vWatershed == "Castner",:);
vCT_EM = vEndmember(vWatershed == "Castner");
vCT_SiteID = vSiteID(vWatershed == "Castner");
mCanwellFull = fullDataset(vWatershed == "Canwell",:);
vCW_EM = vEndmember(vWatershed == "Canwell");
vCW_SiteID = vSiteID(vWatershed == "Canwell");
mGulkanaFull = fullDataset(vWatershed == "Gulkana",:);
vG_EM = vEndmember(vWatershed == "Gulkana");
vG_SiteID = vSiteID(vWatershed == "Gulkana");

DeltaWatersheds = ["Castner", "Canwell", "Gulkana"];
index = ismember(vWatershed, DeltaWatersheds);
mDeltasFull = fullDataset(index, :);
index_EM = ismember(vWatershed,DeltaWatersheds);
vDelta_EM = vEndmember(index_EM,:);
index_SiteID = ismember(vWatershed,DeltaWatersheds);
vDelta_SiteID = vSiteID(index_SiteID,:);

%% zscore the data.
[Knik_Z, Knik_mu, Knik_sigma] = zscore(mKnikFull);
[Mat_Z, Mat_mu, Mat_sigma] = zscore(mMatFull);
[LS_Z, LS_mu, LS_sigma] = zscore(mLSFull);
[CW_Z, CW_mu, CW_sigma] = zscore(mCanwellFull);
[CT_Z, CT_mu, CT_sigma] = zscore(mCastnerFull);
[G_Z, G_mu, G_sigma] = zscore(mGulkanaFull);
[Delta_Z,Delta_mu,Delta_sigma] = zscore(mDeltasFull);

%% Perform PCA on Knik Data
[Knik_coeff,Knik_score,Knik_latent,Knik_tsquared,Knik_explained,Knik_mu2] = pca(Knik_Z);
KR1 = Knik_score(vK_SiteID == "KR1",:);
KR2 = Knik_score(vK_SiteID == "KR2",:);
KR3 = Knik_score(vK_SiteID == "KR3",:);
KR4 = Knik_score(vK_SiteID == "KR4",:);
HC = Knik_score(vK_SiteID == "HC",:);
GC = Knik_score(vK_SiteID == "Goat Creek",:);
K_Spring = Knik_score(vK_EM == "Spring",:);
%% plot knik
figure(1)
scatter(KR1(:, 1), KR1(:, 2), 'ko','MarkerFaceColor',[0.80,0.40,0],'DisplayName', 'KR1');
hold on 
scatter(KR2(:,1), KR2(:,2), 'ko','MarkerFaceColor',[0.90,0.60,0],'DisplayName', 'KR2');
scatter(KR3(:,1), KR3(:,2), 'ko','MarkerFaceColor',[0.95,0.90,0.25],'DisplayName', 'KR3');
scatter(KR4(:,1), KR4(:,2), 'ko','MarkerFaceColor',[0,0.60,0.50],'DisplayName', 'KR4');
scatter(HC(:,1), HC(:,2), 'kd','MarkerFaceColor',[0.35,0.70,0.90],'DisplayName', 'Hunter Creek');
scatter(GC(:,1), GC(:,2), 'kd','MarkerFaceColor',[0,0.45,0.70],'DisplayName', 'Goat Creek');
scatter(K_Spring(:,1),K_Spring(:,2),'ks','MarkerFaceColor',[0.80,0.60,0.70],'DisplayName', 'Springs');
xlabel('PC1');  
ylabel('PC2');  
title('PCA Knik Endmember');
legend('Location', 'eastoutside');
grid off;
hold off;
folderName = 'U:/GoA plots/NewPlots';
fileName = 'Knik_Endmember_PCA.svg';
fullFilePath = fullfile(folderName, fileName);
saveas(gcf, fullFilePath, 'svg');

%% PCA on Matanuska Data
[M_coeff,M_score,M_latent,M_tsquared,M_explained,M_mu2] = pca(Mat_Z);
MR1 = M_score(vM_SiteID == "MR1",:);
MR2 = M_score(vM_SiteID == "MR2",:);
MR3 = M_score(vM_SiteID == "MR3",:);
MR4 = M_score(vM_SiteID == "MR4",:);
MR5 = M_score(vM_SiteID == "MR5",:);
MC = M_score(vM_SiteID == "Moose Creek",:);
CC = M_score(vM_SiteID == "Caribou Creek",:);
Hick = M_score(vM_SiteID == "Hick's Creek",:);
Pino = M_score(vM_SiteID == "Pinocle Creek",:);
Pur = M_score(vM_SiteID == "Purinton Creek",:);
Chick = M_score(vM_SiteID == "Chickaloon River",:);
KR = M_score(vM_SiteID == "Kings River",:);
YJ = M_score(vM_SiteID == "Yellow Jacket Creek",:);
UM = M_score(vM_SiteID == "Upper Matanuska River",:);
M_Spring = M_score(vM_EM == "Spring",:);
%% plot mat
figure(1)
scatter(MR1(:, 1), MR1(:, 2), 'ko','MarkerFaceColor','r','DisplayName', 'MR1');
hold on 
scatter(MR2(:,1), MR2(:,2), 'ko','MarkerFaceColor',[0.90,0.60,0],'DisplayName', 'MR2');
scatter(MR3(:,1), MR3(:,2), 'ko','MarkerFaceColor',[0.95,0.90,0.25],'DisplayName', 'MR3');
scatter(MR4(:,1), MR4(:,2), 'ko','MarkerFaceColor',[0,0.60,0.50],'DisplayName', 'MR4');
scatter(MR5(:,1), MR5(:,2), 'ko','MarkerFaceColor',[0.35,0.70,0.90],'DisplayName', 'MR5');
scatter(MC(:,1), MC(:,2), 'kd','MarkerFaceColor','r','DisplayName', 'Moose Creek');
scatter(CC(:,1), CC(:,2), 'kd','MarkerFaceColor',[0.90,0.60,0],'DisplayName', 'Caribou Creek');
scatter(Hick(:,1), Hick(:,2), 'kd','MarkerFaceColor',[0.95,0.90,0.25],'DisplayName', 'Hicks Creek');
scatter(Pino(:,1), Pino(:,2), 'kd','MarkerFaceColor',[0.6,0.9,0],'DisplayName', 'Pinocle Creek');
scatter(Pur(:,1), Pur(:,2), 'kd','MarkerFaceColor',[0,0.60,0.50],'DisplayName', 'Purinton Creek');
scatter(Chick(:,1), Chick(:,2), 'kd','MarkerFaceColor',[0.35,0.70,0.90],'DisplayName', 'Chickaloon River');
scatter(KR(:,1), KR(:,2), 'kd','MarkerFaceColor',[0,0.45,0.70],'DisplayName', 'Kings River');
scatter(YJ(:,1), YJ(:,2), 'kd','MarkerFaceColor',[0.52,0,0.66],'DisplayName', 'Yellow Jacket Creek');
scatter(UM(:,1), UM(:,2), 'kd','MarkerFaceColor',[1,0.75,0.91],'DisplayName', 'Upper Matanuska River');
scatter(M_Spring(:,1),M_Spring(:,2),'ks','MarkerFaceColor',[0.45,0,0],'DisplayName', 'Springs');
xlabel('PC1');  
ylabel('PC2');  
title('PCA Matanuska Endmember');
legend('Location', 'eastoutside');
grid off;
hold off;
folderName = 'U:/GoA plots/NewPlots';
fileName = 'Mat_Endmember_PCA.svg';
fullFilePath = fullfile(folderName, fileName);
saveas(gcf, fullFilePath, 'svg');

%% PCA on LS Data
[LS_coeff,LS_score,LS_latent,LS_tsquared,LS_explained,LS_mu2] = pca(LS_Z);
LS1 = LS_score(vLS_SiteID == "LS1",:);
LS15 = LS_score(vLS_SiteID == "LS1.5",:);
LS2 = LS_score(vLS_SiteID == "LS2",:);
LS3 = LS_score(vLS_SiteID == "LS3",:);
LS4 = LS_score(vLS_SiteID == "LS4",:);
AR = LS_score(vLS_SiteID == "Archangel River",:);
IL = LS_score(vLS_SiteID == "Ivory Lake",:);
MGS = LS_score(vLS_SiteID == "Mint Glacier Subglacial",:);
SS = LS_score(vLS_SiteID == "Snowbird Supraglacial",:);
SP = LS_score(vLS_SiteID == "Snowbird Periglacial",:);
FC = LS_score(vLS_SiteID == "Fishhook Creek",:);
LS_Spring = LS_score(vLS_EM == "Spring",:);
%% plot ls
figure(1)
scatter(LS1(:, 1), LS1(:, 2), 'ko','MarkerFaceColor',[0.8,0.60,0.7],'DisplayName', 'LS1');
hold on 
scatter(LS15(:,1), LS15(:,2), 'ko','MarkerFaceColor',[0.90,0.60,0],'DisplayName', 'LS1.5');
scatter(LS2(:,1), LS2(:,2), 'ko','MarkerFaceColor',[0.95,0.90,0.25],'DisplayName', 'LS2');
scatter(LS3(:,1), LS3(:,2), 'ko','MarkerFaceColor',[0,0.60,0.50],'DisplayName', 'LS3');
scatter(LS4(:,1), LS4(:,2), 'ko','MarkerFaceColor',[0.35,0.70,0.90],'DisplayName', 'LS4');
scatter(AR(:,1), AR(:,2), 'kd','MarkerFaceColor',[0.8,0.60,0.7],'DisplayName', 'Archangel River');
scatter(FC(:,1), FC(:,2), 'kd','MarkerFaceColor',[0.35,0.70,0.90],'DisplayName', 'Fishhook Creek');
scatter(IL(:,1), IL(:,2), 'k^','MarkerFaceColor',[0.8,0.60,0.7],'DisplayName', 'Ivory Glacial Lake');
scatter(MGS(:,1), MGS(:,2), 'k^','MarkerFaceColor',[0.90,0.60,0],'DisplayName', 'Mint Glacier Subglacial');
scatter(SS(:,1), SS(:,2), 'k^','MarkerFaceColor',[0.95,0.90,0.25],'DisplayName', 'Snowbird Supraglacial');
scatter(SP(:,1), SP(:,2), 'k^','MarkerFaceColor',[0,0.60,0.50],'DisplayName', 'Snowbird Periglacial');
scatter(LS_Spring(:,1),LS_Spring(:,2),'ks','MarkerFaceColor',[0.8,0.60,0.7],'DisplayName', 'Springs');
xlabel('PC1');  
ylabel('PC2');  
title('PCA Little Susitna Endmember');
legend('Location', 'eastoutside');
grid off;
hold off;
folderName = 'U:/GoA plots/NewPlots';
fileName = 'LS_Endmember_PCA.svg';
fullFilePath = fullfile(folderName, fileName);
saveas(gcf, fullFilePath, 'svg');

%% PCA Gulkana 
[G_coeff,G_score,G_latent,G_tsquared,G_explained,G_mu2] = pca(G_Z);
G1 = G_score(vG_SiteID == "G1",:);
G2 = G_score(vG_SiteID == "G2",:);
G25 = G_score(vG_SiteID == "G2.5",:);
G3 = G_score(vG_SiteID == "G3",:);
GS = G_score(vG_SiteID == "Gulkana Supraglacial 1",:);
CC = G_score(vG_SiteID == "College Creek",:);
%% plot gulk
figure(1)
scatter(G1(:, 1), G1(:, 2), 'ko','MarkerFaceColor',[0.8,0.60,0.7],'DisplayName', 'G1');
hold on 
scatter(G2(:,1), G2(:,2), 'ko','MarkerFaceColor',[0.90,0.60,0],'DisplayName', 'G2');
scatter(G25(:,1), G25(:,2), 'ko','MarkerFaceColor',[0.95,0.90,0.25],'DisplayName', 'G2.5');
scatter(G3(:,1), G3(:,2), 'ko','MarkerFaceColor',[0,0.60,0.50],'DisplayName', 'G3');
scatter(GS(:,1), GS(:,2), 'k^','MarkerFaceColor',[0.35,0.70,0.90],'DisplayName', 'Gulkana Supraglacial');
scatter(CoC(:,1), CoC(:,2), 'kd','MarkerFaceColor',[0,0.45,0.7],'DisplayName', 'College Creek');
xlabel('PC1');  
ylabel('PC2');  
title('PCA Gulkana Endmember');
legend('Location', 'eastoutside');
grid off;
hold off;
folderName = 'U:/GoA plots/NewPlots';
fileName = 'Gulkana_Endmember_PCA.svg';
fullFilePath = fullfile(folderName, fileName);
saveas(gcf, fullFilePath, 'svg');

%% PCA Castner 
[CT_coeff,CT_score,CT_latent,CT_tsquared,CT_explained,CT_mu2] = pca(CT_Z);
CT1 = CT_score(vCT_SiteID == "CT1",:);
CT2 = CT_score(vCT_SiteID == "CT2",:);
CS = CT_score(vCT_SiteID == "Castner subglacial",:);
%% plots castner
figure(1)
scatter(CT1(:, 1), CT1(:, 2), 'ko','MarkerFaceColor',[0.8,0.60,0.7],'DisplayName', 'CT1');
hold on 
scatter(CT2(:,1), CT2(:,2), 'ko','MarkerFaceColor',[0.90,0.60,0],'DisplayName', 'CT2');
scatter(CS(:,1), CS(:,2), 'k^','MarkerFaceColor',[0,0.60,0.50],'DisplayName', 'Castner Subglacial');
xlabel('PC1');  
ylabel('PC2');  
title('PCA Castner Endmember');
legend('Location', 'eastoutside');
grid off;
hold off;
folderName = 'U:/GoA plots/NewPlots';
fileName = 'Castner_Endmember_PCA.svg';
fullFilePath = fullfile(folderName, fileName);
saveas(gcf, fullFilePath, 'svg');

%% PCA Deltas
[D_coeff,D_score,D_latent,D_tsquared,D_explained,D_mu2] = pca(Delta_Z);
CT1 = D_score(vDelta_SiteID == "CT1",:);
CT2 = D_score(vDelta_SiteID == "CT2",:);
CS = D_score(vDelta_SiteID == "Castner subglacial",:);
G1 = D_score(vDelta_SiteID == "G1",:);
G2 = D_score(vDelta_SiteID == "G2",:);
G25 = D_score(vDelta_SiteID == "G2.5",:);
G3 = D_score(vDelta_SiteID == "G3",:);
GS = D_score(vDelta_SiteID == "Gulkana Supraglacial 1",:);
CoC = D_score(vDelta_SiteID == "College Creek",:);
CW1 = D_score(vDelta_SiteID == "CW1",:);
CW2 = D_score(vDelta_SiteID == "CW2",:);
%% plot deltas
figure(1)
scatter(CT1(:, 1), CT1(:, 2), 'ko','MarkerFaceColor',[0.66,0,0],'DisplayName', 'CT1');
hold on 
scatter(CT2(:,1), CT2(:,2), 'ko','MarkerFaceColor',[0.90,0.60,0],'DisplayName', 'CT2');
scatter(G1(:, 1), G1(:, 2), 'ko','MarkerFaceColor',[0.95,0.90,0.25],'DisplayName', 'G1');
scatter(G2(:,1), G2(:,2), 'ko','MarkerFaceColor',[0,0.60,0.50],'DisplayName', 'G2')
scatter(G25(:,1), G25(:,2), 'ko','MarkerFaceColor',[0.35,0.70,0.90],'DisplayName', 'G2.5');
scatter(G3(:,1), G3(:,2), 'ko','MarkerFaceColor',[0,0.45,0.7],'DisplayName', 'G3');
scatter(CW1(:, 1), CW1(:, 2), 'ko','MarkerFaceColor',[0.8,0.60,0.7],'DisplayName', 'CW1');
scatter(CW2(:, 1), CW2(:, 2), 'ko','MarkerFaceColor',[0.67,0.4,0.8],'DisplayName', 'CW2');
scatter(CS(:,1), CS(:,2), 'k^','MarkerFaceColor',[0.66,0,0],'DisplayName', 'Castner Subglacial');
scatter(GS(:,1), GS(:,2), 'k^','MarkerFaceColor',[0.90,0.60,0],'DisplayName', 'Gulkana Supraglacial');
scatter(CoC(:,1), CoC(:,2), 'kd','MarkerFaceColor',[0.95,0.90,0.25],'DisplayName', 'College Creek');
xlabel('PC1');  
ylabel('PC2');  
title('PCA Deltas Endmember');
legend('Location', 'eastoutside');
grid off;
hold off;
folderName = 'U:/GoA plots/NewPlots';
fileName = 'Deltas_Endmember_PCA.svg';
fullFilePath = fullfile(folderName, fileName);
saveas(gcf, fullFilePath, 'svg');

%% combined plots
figure;
legendFontSize = 7;
markerSize = 30;
markerSizeS = 35;
position1 = [0.1, 0.55, 0.4, 0.4];
position2 = [0.55, 0.55, 0.4, 0.4];
position3 = [0.1, 0.1, 0.4, 0.4];
position4 = [0.55, 0.1, 0.4, 0.4];
subplot('Position', position1); % Knik plot
scatter(KR1(:, 1), KR1(:, 2), 'o','MarkerEdgeColor',[0.80,0.60,0.70],'DisplayName', 'KR1','SizeData',markerSize);
hold on 
scatter(KR2(:,1), KR2(:,2), 'o','MarkerEdgeColor',[0.90,0.60,0],'DisplayName', 'KR2','SizeData',markerSize);
scatter(KR3(:,1), KR3(:,2), 'o','MarkerEdgeColor',[0.95,0.90,0.25],'DisplayName', 'KR3','SizeData',markerSize);
scatter(KR4(:,1), KR4(:,2), 'o','MarkerEdgeColor',[0,0.60,0.50],'DisplayName', 'KR4','SizeData',markerSize);
scatter(HC(:,1), HC(:,2), 'kd','MarkerFaceColor',[0.35,0.70,0.90],'DisplayName', 'Hunter Creek','SizeData',markerSize);
scatter(GC(:,1), GC(:,2), 'kd','MarkerFaceColor',[0,0.45,0.70],'DisplayName', 'Goat Creek','SizeData',markerSize);
scatter(K_Spring(:,1),K_Spring(:,2),'ks','MarkerFaceColor',[0.80,0.60,0.70],'DisplayName', 'Springs','SizeData',markerSizeS);
%xlabel('PC1');  
ylabel('PC2');  
%title('Knik Watershed');
legend('Location', 'NorthEast','FontSize',legendFontSize);
grid off;
hold off;

subplot('Position', position2); % Matanuska plot
markerSize1=30;
scatter(MR1(:, 1), MR1(:, 2), 'o','MarkerEdgeColor',[0.80,0.60,0.70],'DisplayName', 'MR1','SizeData',markerSize1);
hold on 
scatter(MR2(:,1), MR2(:,2), 'o','MarkerEdgeColor',[0.90,0.60,0],'DisplayName', 'MR2','SizeData',markerSize1);
scatter(MR3(:,1), MR3(:,2), 'o','MarkerEdgeColor',[0.95,0.90,0.25],'DisplayName', 'MR3','SizeData',markerSize1);
scatter(MR4(:,1), MR4(:,2), 'o','MarkerEdgeColor',[0,0.60,0.50],'DisplayName', 'MR4','SizeData',markerSize1);
scatter(MR5(:,1), MR5(:,2), 'o','MarkerEdgeColor',[0.52,0,0.66],'DisplayName', 'MR5','SizeData',markerSize1);
scatter(MC(:,1), MC(:,2), 'kd','MarkerFaceColor','r','DisplayName', 'Moose Creek');
scatter(CC(:,1), CC(:,2), 'kd','MarkerFaceColor',[0.90,0.60,0],'DisplayName', 'Caribou Creek','SizeData',markerSize1);
scatter(Hick(:,1), Hick(:,2), 'kd','MarkerFaceColor',[0.95,0.90,0.25],'DisplayName', 'Hicks Creek','SizeData',markerSize1);
scatter(Pino(:,1), Pino(:,2), 'kd','MarkerFaceColor',[0.6,0.9,0],'DisplayName', 'Pinocle Creek','SizeData',markerSize1);
scatter(Pur(:,1), Pur(:,2), 'kd','MarkerFaceColor',[0,0.60,0.50],'DisplayName', 'Purinton Creek','SizeData',markerSize1);
scatter(Chick(:,1), Chick(:,2), 'kd','MarkerFaceColor',[0.35,0.70,0.90],'DisplayName', 'Chickaloon River','SizeData',markerSize1);
scatter(KR(:,1), KR(:,2), 'kd','MarkerFaceColor',[0,0.45,0.70],'DisplayName', 'Kings River','SizeData',markerSize1);
scatter(YJ(:,1), YJ(:,2), 'kd','MarkerFaceColor',[0.52,0,0.66],'DisplayName', 'Yellow Jacket Creek','SizeData',markerSize1);
scatter(UM(:,1), UM(:,2), 'kd','MarkerFaceColor',[1,0.75,0.91],'DisplayName', 'Upper Matanuska River','SizeData',markerSize1);
scatter(M_Spring(:,1),M_Spring(:,2),'ks','MarkerFaceColor',[0.45,0,0],'DisplayName', 'Springs','SizeData',markerSizeS);
%xlabel('PC1');  
%ylabel('PC2');  
%title('Matanuska Watershed');
legend('Location', 'NorthEast','FontSize',legendFontSize);
grid off;
hold off;

subplot('Position', position3); %LS plot
scatter(LS1(:, 1), LS1(:, 2), 'o','MarkerEdgeColor',[0.80,0.60,0.70],'DisplayName', 'LS1');
hold on 
scatter(LS15(:,1), LS15(:,2), 'o','MarkerEdgeColor',[0.90,0.60,0],'DisplayName', 'LS1.5');
scatter(LS2(:,1), LS2(:,2), 'o','MarkerEdgeColor',[0.95,0.90,0.25],'DisplayName', 'LS2');
scatter(LS3(:,1), LS3(:,2), 'o','MarkerEdgeColor',[0,0.60,0.50],'DisplayName', 'LS3');
scatter(LS4(:,1), LS4(:,2), 'o','MarkerEdgeColor',[0.52,0,0.66],'DisplayName', 'LS4');
scatter(AR(:,1), AR(:,2), 'kd','MarkerFaceColor',[0.8,0.60,0.7],'DisplayName', 'Archangel River','SizeData',markerSize);
scatter(FC(:,1), FC(:,2), 'kd','MarkerFaceColor',[0.35,0.70,0.90],'DisplayName', 'Fishhook Creek','SizeData',markerSize);
scatter(IL(:,1), IL(:,2), 'k^','MarkerFaceColor',[0.8,0.60,0.7],'DisplayName', 'Ivory Glacial Lake','SizeData',markerSize);
scatter(MGS(:,1), MGS(:,2), 'k^','MarkerFaceColor',[0.90,0.60,0],'DisplayName', 'Mint Glacier Subglacial','SizeData',markerSize);
scatter(SS(:,1), SS(:,2), 'k^','MarkerFaceColor',[0.95,0.90,0.25],'DisplayName', 'Snowbird Supraglacial','SizeData',markerSize);
scatter(SP(:,1), SP(:,2), 'k^','MarkerFaceColor',[0,0.60,0.50],'DisplayName', 'Snowbird Periglacial','SizeData',markerSize);
scatter(LS_Spring(:,1),LS_Spring(:,2),'ks','MarkerFaceColor',[0.8,0.60,0.7],'DisplayName', 'Springs','SizeData',markerSizeS);
xlabel('PC1');  
ylabel('PC2');  
%title('Little Susitna Watersheds');
legend('Location', 'NorthEast','FontSize',legendFontSize);
grid off;
hold off;

subplot('Position', position4); %Delta's Plot
scatter(CT1(:, 1), CT1(:, 2), 'o','MarkerEdgeColor','r','DisplayName', 'CT1','SizeData',markerSize);
hold on 
scatter(CT2(:,1), CT2(:,2), 'o','MarkerEdgeColor',[0.90,0.60,0],'DisplayName', 'CT2','SizeData',markerSize);
scatter(G1(:, 1), G1(:, 2), 'o','MarkerEdgeColor',[0.95,0.90,0.25],'DisplayName', 'G1','SizeData',markerSize);
scatter(G2(:,1), G2(:,2), 'o','MarkerEdgeColor',[0,0.60,0.50],'DisplayName', 'G2','SizeData',markerSize);
scatter(G25(:,1), G25(:,2), 'o','MarkerEdgeColor',[0.35,0.70,0.90],'DisplayName', 'G2.5','SizeData',markerSize);
scatter(G3(:,1), G3(:,2), 'o','MarkerEdgeColor',[0,0.45,0.7],'DisplayName', 'G3','SizeData',markerSize);
scatter(CW1(:, 1), CW1(:, 2), 'o','MarkerEdgeColor',[0.8,0.60,0.7],'DisplayName', 'CW1','SizeData',markerSize);
scatter(CW2(:, 1), CW2(:, 2), 'o','MarkerEdgeColor',[0.52,0,0.66],'DisplayName', 'CW2','SizeData',markerSize);
scatter(CS(:,1), CS(:,2), 'k^','MarkerFaceColor',[0.80,0.60,0.70],'DisplayName', 'Castner Subglacial','SizeData',markerSize);
scatter(GS(:,1), GS(:,2), 'k^','MarkerFaceColor',[0.90,0.60,0],'DisplayName', 'Gulkana Supraglacial','SizeData',markerSize);
scatter(CoC(:,1), CoC(:,2), 'kd','MarkerFaceColor',[0.95,0.90,0.25],'DisplayName', 'College Creek','SizeData',markerSize);
xlabel('PC1');  
%ylabel('PC2');  
%title('Delta Watersheds');
legend('Location', 'NorthEast','FontSize',legendFontSize);
grid off;
hold off;
set(gcf, 'Position', [100, 100, 1200, 800]);


% Save the combined figure
folderName = 'U:/GoA plots/NewPlots';
fileName = 'Endmember_PCA_open.jpg';
fullFilePath = fullfile(folderName, fileName);
saveas(gcf, fullFilePath, 'jpg');
