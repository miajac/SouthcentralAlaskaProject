%% Endmember July Data
clear all, close all, clc
[sFile, sPath] = uigetfile('*.xlsx', 'Select Database File'); %endmemberJulyAug2023
sFullPath = fullfile(sPath, sFile);
mC = readtable(sFullPath);
mA = mC(:, 6:end);
fullDataset = table2array(mA);
vWatershed = string(table2cell(mC(:, 3)));
vEndmember = string(table2cell(mC(:,4)));
vSiteID = string(table2cell(mC(:,2)));
vLabID = string(table2cell(mC(:,1)));
%% Define Endmember Categories
mKnikFull = fullDataset(vWatershed == "Knik" | vWatershed == "all",:);
vK_EM = vEndmember(vWatershed == "Knik" | vWatershed == "all");
vK_SiteID = vSiteID(vWatershed == "Knik" | vWatershed == "all");
vK_LabID = vLabID(vWatershed == "Knik" | vWatershed == "all");
mMatFull = fullDataset(vWatershed == "Matanuska" | vWatershed == "all",:);
vM_EM = vEndmember(vWatershed == "Matanuska" | vWatershed == "all");
vM_SiteID = vSiteID(vWatershed == "Matanuska" | vWatershed == "all");
vM_LabID = vLabID(vWatershed == "Matanuska" | vWatershed == "all");
mLSFull = fullDataset(vWatershed == "Little Susitna" | vWatershed == "all",:);
vLS_EM = vEndmember(vWatershed == "Little Susitna" | vWatershed == "all");
vLS_SiteID = vSiteID(vWatershed == "Little Susitna" | vWatershed == "all");
vLS_LabID = vLabID(vWatershed == "Little Susitna" | vWatershed == "all");
mCastnerFull = fullDataset(vWatershed == "Castner" | vWatershed == "all",:);
vCT_EM = vEndmember(vWatershed == "Castner" | vWatershed == "all");
vCT_SiteID = vSiteID(vWatershed == "Castner" | vWatershed == "all");
mCanwellFull = fullDataset(vWatershed == "Canwell" | vWatershed == "all",:);
vCW_EM = vEndmember(vWatershed == "Canwell" | vWatershed == "all");
vCW_SiteID = vSiteID(vWatershed == "Canwell" | vWatershed == "all");
mGulkanaFull = fullDataset(vWatershed == "Gulkana" | vWatershed == "all",:);
vG_EM = vEndmember(vWatershed == "Gulkana" | vWatershed == "all");
vG_SiteID = vSiteID(vWatershed == "Gulkana" | vWatershed == "all");

DeltaWatersheds = ["Castner", "Canwell", "Gulkana", "all"];
index = ismember(vWatershed, DeltaWatersheds);
mDeltasFull = fullDataset(index, :);
index_EM = ismember(vWatershed,DeltaWatersheds);
vDelta_EM = vEndmember(index_EM,:);
index_SiteID = ismember(vWatershed,DeltaWatersheds);
vDelta_SiteID = vSiteID(index_SiteID,:);
vDelta_LabID = vLabID(index_SiteID,:);
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
KR2 = Knik_score(vK_SiteID == "KR2",:);
KR3 = Knik_score(vK_SiteID == "KR3",:);
KR4 = Knik_score(vK_SiteID == "KR4",:);
K_Spring = Knik_score(vK_EM == "Spring",:);
K_Trib = Knik_score(vK_EM == 'Tributary',:);
K_G = Knik_score(vK_EM == 'Glacier',:);
K_Precip = Knik_score(vK_EM == 'precip',:);
%comment out if needed
%K_GSupra = Knik_score(vK_SiteID == "Gulkana Supraglacial 1",:);
%K_Den = Knik_score(vK_SiteID == "Denali NP",:);
%K_Kat = Knik_score(vK_SiteID == "Katmai NP",:);
%K_SSupra = Knik_score(vK_SiteID == "Snowbird Supraglacial",:);
%K_SPeri = Knik_score(vK_SiteID == "Snowbird Periglacial",:);
%HC = Knik_score(vK_SiteID == "HC",:);
%GC = Knik_score(vK_SiteID == "Goat Creek",:);
%K_MR4 = Knik_score(vK_SiteID == "Mat",:);
%KR1 = Knik_score(vK_SiteID == "KR1",:);
%% PCA on Matanuska Data
[M_coeff,M_score,M_latent,M_tsquared,M_explained,M_mu2] = pca(Mat_Z);
MR1 = M_score(vM_SiteID == "MR1",:);
MR2 = M_score(vM_SiteID == "MR2",:);
MR3 = M_score(vM_SiteID == "MR3",:);
MR4 = M_score(vM_SiteID == "MR4",:);
MR5 = M_score(vM_SiteID == "MR5",:);
M_Spring = M_score(vM_EM == "Spring",:);
M_Trib = M_score(vM_EM == "Tributary",:);
M_G = M_score(vM_EM == "Glacier",:);
M_Precip = M_score(vM_EM == "precip",:);
%comment out if needed
% M_GSupra = M_score(vM_SiteID == "Gulkana Supraglacial 1",:);
% M_Den = M_score(vM_SiteID == "Denali NP",:);
% M_Kat = M_score(vM_SiteID == "Katmai NP",:);
% M_SSupra = M_score(vM_SiteID == "Snowbird Supraglacial",:);
% M_SPeri = M_score(vM_EM == "Snowbird Periglacial",:);
% MC = M_score(vM_SiteID == "MC",:);
% CC = M_score(vM_SiteID == "CC",:);
% Hick = M_score(vM_SiteID == "Hick's Creek",:);
% Pino = M_score(vM_SiteID == "Pinocle Creek",:);
% Pur = M_score(vM_SiteID == "Purinton Creek",:);
% Chick = M_score(vM_SiteID == "Chickaloon River",:);
% KR = M_score(vM_SiteID == "Kings River",:);
% YJ = M_score(vM_SiteID == "Yellow Jacket Creek",:);
% UM = M_score(vM_SiteID == "Upper Matanuska River",:);
%% PCA on LS Data
[LS_coeff,LS_score,LS_latent,LS_tsquared,LS_explained,LS_mu2] = pca(LS_Z);
LS1 = LS_score(vLS_SiteID == "LS1",:);
LS15 = LS_score(vLS_SiteID == "LS1.5",:);
LS2 = LS_score(vLS_SiteID == "LS2",:);
LS3 = LS_score(vLS_SiteID == "LS3",:);
LS4 = LS_score(vLS_SiteID == "LS4",:);
LS_Spring = LS_score(vLS_EM == "Spring",:);
LS_Trib = LS_score(vLS_EM == "Tributary",:);
LS_G = LS_score(vLS_EM == "Glacier",:);
LS_P = LS_score(vLS_EM == "precip",:);
%comment out if needed
% LS_GSupra = LS_score(vLS_SiteID == "Gulkana Supraglacial 1",:);
% LS_Den = LS_score(vLS_SiteID == "Denali NP",:);
% LS_Kat = LS_score(vLS_SiteID == "Katmai NP",:);
% AR = LS_score(vLS_SiteID == "Archangel River",:);
% IL = LS_score(vLS_SiteID == "Ivory Lake",:);
% MGS = LS_score(vLS_SiteID == "Mint Glacier Subglacial",:);
% SS = LS_score(vLS_SiteID == "Snowbird Supraglacial",:);
% SP = LS_score(vLS_SiteID == "Snowbird Periglacial",:);
% FC = LS_score(vLS_SiteID == "Fishhook Creek",:);
%% PCA Gulkana 
[G_coeff,G_score,G_latent,G_tsquared,G_explained,G_mu2] = pca(G_Z);
G1 = G_score(vG_SiteID == "G1",:);
G2 = G_score(vG_SiteID == "G2",:);
GS = G_score(vG_SiteID == "Gulkana Supraglacial 1",:);
CC = G_score(vG_SiteID == "College Creek",:);
%% PCA Castner 
[CT_coeff,CT_score,CT_latent,CT_tsquared,CT_explained,CT_mu2] = pca(CT_Z);
CT1 = CT_score(vCT_SiteID == "CT1",:);
CT2 = CT_score(vCT_SiteID == "CT2",:);
CS = CT_score(vCT_SiteID == "Castner subglacial",:);
%% PCA Deltas
[D_coeff,D_score,D_latent,D_tsquared,D_explained,D_mu2] = pca(Delta_Z);
CT1 = D_score(vDelta_SiteID == "CT1",:);
CT2 = D_score(vDelta_SiteID == "CT2",:);
G1 = D_score(vDelta_SiteID == "G1",:);
G2 = D_score(vDelta_SiteID == "G2",:);
CW1 = D_score(vDelta_SiteID == "CW1",:);
CW2 = D_score(vDelta_SiteID == "CW2",:);
D_Trib = D_score(vDelta_EM == "Tributary",:);
D_Spring = D_score(vDelta_EM == "Spring",:);
D_G = D_score(vDelta_EM == "Glacier",:);
D_P = D_score(vDelta_EM == "precip",:);
%comment out if needed
% D_Den = D_score(vDelta_SiteID == "Denali NP",:);
% D_Kat = D_score(vDelta_SiteID == "Katmai NP",:);
% D_SSupra = D_score(vDelta_SiteID == "Snowbird Supraglacial",:);
% D_SPeri = D_score(vDelta_SiteID == "Snowbird Periglacial",:);
% CS = D_score(vDelta_SiteID == "Castner subglacial",:);
% GS = D_score(vDelta_SiteID == "Gulkana Supraglacial 1",:);
% CoC = D_score(vDelta_SiteID == "College Creek",:);
%% combined plots
figure;
legendFontSize = 10;
legendFontSize1 = 8;
markerSize = 75;
markerSizeD = 40;
markerSizeS = 75;
markerSizeG = 45;
l=2.5; %line width
position1 = [0.07, 0.57, 0.4, 0.4];  % Slightly adjust positions to reduce space
position2 = [0.55, 0.57, 0.4, 0.4];
position3 = [0.07, 0.1, 0.4, 0.4];
position4 = [0.55, 0.1, 0.4, 0.4];
subplot('Position', position1); % Knik plot
%scatter(KR1(:, 1), KR1(:, 2), 'o','MarkerEdgeColor',[0.80,0.60,0.70],'DisplayName', 'KR1','SizeData',markerSize);
%hold on 
scatter(KR2(:,1), KR2(:,2), 'ko','MarkerFaceColor',[0.992, 0.851, 0.627],'DisplayName', 'KR2','SizeData',markerSize, 'LineWidth', 0.5, 'HandleVisibility', 'on');
hold on
scatter(KR3(:,1), KR3(:,2), 'ko','MarkerFaceColor',[0.961, 0.498, 0.090],'DisplayName', 'KR3','SizeData',markerSize, 'LineWidth',0.5, 'HandleVisibility', 'on');
scatter(KR4(:,1), KR4(:,2), 'ko','MarkerFaceColor',[0.651, 0.212, 0.012],'DisplayName', 'KR4','SizeData',markerSize, 'LineWidth', 0.5, 'HandleVisibility', 'on');
scatter(K_Spring(:,1),K_Spring(:,2),'kd','MarkerFaceColor',[0, 0.620, 0.451],'DisplayName', 'Springs','SizeData',markerSizeD);
scatter(K_Trib(:,1),K_Trib(:,2),'ks','MarkerFaceColor',[1.000, 0.843, 0.000],'DisplayName', 'Tributaries','SizeData',markerSizeS);
scatter(K_G(:,1), K_G(:,2), 'k^','MarkerFaceColor',[0.616, 0.765, 0.902],'DisplayName','Glacial','SizeData',markerSizeG);
scatter(K_Precip(:,1), K_Precip(:,2), 'k*','DisplayName','Precipitation','SizeData',markerSizeG, 'LineWidth',1);

% scatter(HC(:,1), HC(:,2), 'ks','MarkerFaceColor',[1.000, 0.843, 0.000],'DisplayName', 'Tributaries','SizeData',markerSizeS);
% scatter(GC(:,1), GC(:,2), 'ks','MarkerFaceColor',[1.000, 0.843, 0.000],'DisplayName', 'Goat Creek','SizeData',markerSizeS, 'HandleVisibility', 'off'); 
% scatter(K_MR4(:,1),K_MR4(:,2),'ks','MarkerFaceColor', [1.000, 0.843, 0.000],'DisplayName', 'Matanuska','SizeData',markerSizeS, 'HandleVisibility', 'off');
% scatter(K_GSupra(:,1), K_GSupra(:,2), 'k^','MarkerFaceColor',[0.616, 0.765, 0.902],'DisplayName','Glacial','SizeData',markerSizeG,'HandleVisibility','on');
% scatter(K_SSupra(:,1), K_SSupra(:,2), 'k^','MarkerFaceColor',[0.616, 0.765, 0.902],'DisplayName','Supraglacial','SizeData',markerSizeG, 'HandleVisibility','off');
% scatter(K_SPeri(:,1), K_SPeri(:,2), 'k^','MarkerFaceColor',[0.616, 0.765, 0.902],'DisplayName','Periglacial','SizeData',markerSizeG, 'HandleVisibility','off');
% scatter(K_Den(:,1), K_Den(:,2), 'k*','DisplayName','Precipitation','SizeData',markerSizeG, 'LineWidth',1, 'HandleVisibility','on');
% scatter(K_Kat(:,1), K_Kat(:,2), 'k*','DisplayName','Precipitation','SizeData',markerSizeG, 'LineWidth',1, 'HandleVisibility','off');

%xlabel('PC1');  
ylabel('PC2'); 
ax = gca;
ax.FontSize = 10; 
title('Knik', 'FontSize',14);
legend('Location', 'NorthWest','FontSize',legendFontSize1);
grid off;
hold off;

subplot('Position', position2); % Matanuska plot
markerSize1=75;
scatter(MR1(:, 1), MR1(:, 2), 'ko','MarkerFaceColor',[0.992, 0.851, 0.627],'DisplayName', 'MR1','SizeData',markerSize, 'LineWidth', 0.5, 'HandleVisibility', 'on'); 
hold on
scatter(MR2(:,1), MR2(:,2), 'ko','MarkerFaceColor',[0.992, 0.722, 0.388],'DisplayName', 'MR2','SizeData',markerSize, 'LineWidth', 0.5, 'HandleVisibility', 'on');
scatter(MR3(:,1), MR3(:,2), 'ko','MarkerFaceColor',[0.961, 0.498, 0.090],'DisplayName', 'MR3','SizeData',markerSize, 'LineWidth', 0.5, 'HandleVisibility', 'on');
scatter(MR4(:,1), MR4(:,2), 'ko','MarkerFaceColor',[0.902, 0.333, 0.051],'DisplayName', 'MR4','SizeData',markerSize, 'LineWidth', 0.5, 'HandleVisibility', 'on');
scatter(MR5(:,1), MR5(:,2), 'ko','MarkerFaceColor',[0.651, 0.212, 0.012],'DisplayName', 'MR5','SizeData',markerSize, 'LineWidth', 0.5, 'HandleVisibility', 'on');
scatter(M_Spring(:,1),M_Spring(:,2),'kd','MarkerFaceColor',[0, 0.620, 0.451],'DisplayName', 'Springs','SizeData',markerSizeD);
scatter(M_Trib(:,1), M_Trib(:,2), 'ks','MarkerFaceColor',[1.000, 0.843, 0.000],'DisplayName', 'Tributaries','SizeData',markerSizeS);
scatter(M_G(:,1),M_G(:,2), 'k^','MarkerFaceColor',[0.616, 0.765, 0.902],'DisplayName','Glacial','SizeData',markerSizeG);
scatter(M_Precip(:,1),M_Precip(:,2), 'k*','DisplayName','Precipitation','SizeData',markerSizeG,'LineWidth',1);

% scatter(YJ(:,1), YJ(:,2), 'ks','MarkerFaceColor',[1.000, 0.843, 0.000],'DisplayName', 'Tributaries','SizeData',markerSizeS);
% scatter(CC(:,1), CC(:,2), 'ks','MarkerFaceColor',[1.000, 0.843, 0.000],'DisplayName', 'Caribou Creek','SizeData',markerSizeS, 'HandleVisibility', 'off');
% scatter(UM(:,1), UM(:,2), 'ks','MarkerFaceColor',[1.000, 0.843, 0.000],'DisplayName', 'Upper Matanuska River','SizeData',markerSizeS, 'HandleVisibility', 'off');
% scatter(Pino(:,1), Pino(:,2), 'ks','MarkerFaceColor',[1.000, 0.843, 0.000],'DisplayName', 'Pinochle Creek','SizeData',markerSizeS, 'HandleVisibility', 'off');
% scatter(Hick(:,1), Hick(:,2), 'ks','MarkerFaceColor',[1.000, 0.843, 0.000],'DisplayName', 'Hicks Creek','SizeData',markerSizeS, 'HandleVisibility', 'off');
% scatter(Pur(:,1), Pur(:,2), 'ks','MarkerFaceColor',[1.000, 0.843, 0.000],'DisplayName', 'Purinton Creek','SizeData',markerSizeS, 'HandleVisibility', 'off');
% scatter(Chick(:,1), Chick(:,2), 'ks','MarkerFaceColor',[1.000, 0.843, 0.000],'DisplayName', 'Chickaloon River','SizeData',markerSizeS, 'HandleVisibility', 'off');
% scatter(KR(:,1), KR(:,2), 'ks','MarkerFaceColor',[1.000, 0.843, 0.000],'DisplayName', 'Kings River','SizeData',markerSizeS, 'HandleVisibility', 'off');
% scatter(MC(:,1), MC(:,2), 'ks','MarkerFaceColor',[1.000, 0.843, 0.000],'DisplayName', 'Moose Creek','SizeData',markerSizeS, 'HandleVisibility', 'off');
% scatter(M_GSupra(:,1),M_GSupra(:,2), 'k^','MarkerFaceColor',[0.616, 0.765, 0.902],'DisplayName','Glacial','SizeData',markerSizeG,'HandleVisibility','on');
% scatter(M_SSupra(:,1),M_SSupra(:,2), 'k^','MarkerFaceColor',[0.616, 0.765, 0.902],'DisplayName','Supraglacial','SizeData',markerSizeG,'HandleVisibility','off');
% scatter(M_SPeri(:,1),M_SPeri(:,2), 'k^','MarkerFaceColor',[0.616, 0.765, 0.902],'DisplayName','Periglacial','SizeData',markerSizeG,'HandleVisibility','off');
% scatter(M_Den(:,1),M_Den(:,2), 'k*','DisplayName','Precipitation','SizeData',markerSizeG,'LineWidth',1,'HandleVisibility','on');
% scatter(M_Kat(:,1),M_Kat(:,2), 'k*','DisplayName','Precipitation','SizeData',markerSizeG,'LineWidth',1,'HandleVisibility','off');

%xlabel('PC1');  
%ylabel('PC2'); 
ax = gca;
ax.FontSize = 10; 
title('Matanuska', 'FontSize',14);
legendFontSize1 = 8;
legend('Location', 'SouthEast','FontSize',legendFontSize1);
grid off;
hold off;

subplot('Position', position3); %LS plot

scatter(LS1(:, 1), LS1(:, 2), 'ko','MarkerFaceColor',[0.992, 0.851, 0.627],'DisplayName', 'LS1','SizeData',markerSize, 'LineWidth', 0.5, 'HandleVisibility', 'on'); 
hold on
scatter(LS15(:,1), LS15(:,2), 'ko','MarkerFaceColor',[0.992, 0.722, 0.388],'DisplayName', 'LS1.5','SizeData',markerSize, 'LineWidth', 0.5, 'HandleVisibility', 'on');
scatter(LS2(:,1), LS2(:,2), 'ko','MarkerFaceColor',[0.961, 0.498, 0.090],'DisplayName', 'LS2','SizeData',markerSize, 'LineWidth', 0.5, 'HandleVisibility', 'on');
scatter(LS3(:,1), LS3(:,2), 'ko','MarkerFaceColor',[0.902, 0.333, 0.051],'DisplayName', 'LS3','SizeData',markerSize, 'LineWidth', 0.5, 'HandleVisibility', 'on');
scatter(LS4(:,1), LS4(:,2), 'ko','MarkerFaceColor',[0.651, 0.212, 0.012],'DisplayName', 'LS4','SizeData',markerSize, 'LineWidth', 0.5, 'HandleVisibility', 'on');
scatter(LS_Spring(:,1),LS_Spring(:,2),'kd','MarkerFaceColor',[0, 0.620, 0.451],'DisplayName', 'Springs','SizeData',markerSizeD);
scatter(LS_Trib(:,1), LS_Trib(:,2), 'ks','MarkerFaceColor',[1.000, 0.843, 0.000],'DisplayName', 'Tributaries','SizeData',markerSizeS);
scatter(LS_G(:,1), LS_G(:,2), 'k^','MarkerFaceColor',[0.616, 0.765, 0.902],'DisplayName', 'Glacial','SizeData',markerSizeG); 
scatter(LS_P(:,1),LS_P(:,2),'k*','DisplayName','Precipitation','SizeData',markerSizeG);

% scatter(AR(:,1), AR(:,2), 'ks','MarkerFaceColor',[1.000, 0.843, 0.000],'DisplayName', 'Tributaries','SizeData',markerSizeS);
% scatter(FC(:,1), FC(:,2), 'ks','MarkerFaceColor',[1.000, 0.843, 0.000],'DisplayName', 'Fishhook Creek','SizeData',markerSizeS, 'HandleVisibility', 'off');
% scatter(IL(:,1), IL(:,2), 'k^','MarkerFaceColor',[0.616, 0.765, 0.902],'DisplayName', 'Glacial','SizeData',markerSizeG); %Glacial Lake
% scatter(SS(:,1), SS(:,2), 'k^','MarkerFaceColor',[0.616, 0.765, 0.902],'DisplayName','Supraglacial','SizeData',markerSizeG,'HandleVisibility','off');
% scatter(MGS(:,1), MGS(:,2), 'k^','MarkerFaceColor',[0.616, 0.765, 0.902],'DisplayName', 'Subglacial','SizeData',markerSizeG,'HandleVisibility','off');
% scatter(SP(:,1), SP(:,2), 'k^','MarkerFaceColor',[0.616, 0.765, 0.902],'DisplayName','Periglacial','SizeData',markerSizeG, 'HandleVisibility','off');
% scatter(LS_GSupra(:,1),LS_GSupra(:,2),'k^','MarkerFaceColor',[0.616, 0.765, 0.902],'DisplayName','Supraglacial','SizeData',markerSizeG,'HandleVisibility','off');
% scatter(LS_Den(:,1),LS_Den(:,2),'k*','DisplayName','Precipitation','SizeData',markerSizeG,'HandleVisibility','on');
% scatter(LS_Kat(:,1),LS_Kat(:,2),'k*','DisplayName','Precipitation','SizeData',markerSizeG,'HandleVisibility','off');
xlabel('PC1');
ylabel('PC2');
ax = gca;
ax.FontSize = 10; 
title('Little Susitna', 'FontSize',14);
legend('Location', 'NorthWest','FontSize',legendFontSize1);
grid off;
hold off;

subplot('Position', position4); % Delta's Plot  
scatter(CT1(:, 1), CT1(:, 2), 'ko', 'MarkerFaceColor', [1.000, 0.949, 0.800], 'DisplayName', 'CT1', 'SizeData', markerSize, 'LineWidth', 0.5, 'HandleVisibility', 'on');
hold on
scatter(CT2(:, 1), CT2(:, 2), 'ko', 'MarkerFaceColor', [0.961, 0.871, 0.702],'DisplayName', 'CT2', 'SizeData', markerSize, 'LineWidth', 0.5, 'HandleVisibility', 'on');
scatter(G1(:, 1), G1(:, 2), 'ko', 'MarkerFaceColor', [0.992, 0.722, 0.388],'DisplayName', 'G1', 'SizeData', markerSize, 'LineWidth', 0.5, 'HandleVisibility', 'on');
scatter(G2(:, 1), G2(:, 2), 'ko', 'MarkerFaceColor', [0.961, 0.498, 0.090],'DisplayName', 'G2', 'SizeData', markerSize, 'LineWidth', 0.5, 'HandleVisibility', 'on');
scatter(CW1(:, 1), CW1(:, 2), 'ko', 'MarkerFaceColor', [0.800, 0.361, 0.400],'DisplayName', 'CW1', 'SizeData', markerSize, 'LineWidth', 0.5, 'HandleVisibility', 'on');
scatter(CW2(:, 1), CW2(:, 2), 'ko', 'MarkerFaceColor', [0.600, 0.200, 0.251],'DisplayName', 'CW2', 'SizeData', markerSize, 'LineWidth', 0.5, 'HandleVisibility', 'on');
scatter(D_Spring(:,1),D_Spring(:,2),'k^','MarkerFaceColor',[0.616, 0.765, 0.902],'DisplayName', 'Glacial','SizeData',markerSizeG, 'HandleVisibility','off');
scatter(D_Trib(:, 1), D_Trib(:, 2), 'ks', 'MarkerFaceColor', [1.000, 0.843, 0.000],'DisplayName', 'Tributary', 'SizeData', markerSizeS);
scatter(D_G(:, 1), D_G(:, 2), 'k^', 'MarkerFaceColor', [0.616, 0.765, 0.902], 'DisplayName', 'Glacial', 'SizeData', markerSizeG);
scatter(D_P(:,1),D_P(:,2),'k*','DisplayName','Precipitation','SizeData',markerSizeG);

% scatter(CoC(:, 1), CoC(:, 2), 'ks', 'MarkerFaceColor', [1.000, 0.843, 0.000],'DisplayName', 'Tributary', 'SizeData', markerSizeS);
% scatter(GS(:, 1), GS(:, 2), 'k^', 'MarkerFaceColor', [0.616, 0.765, 0.902], 'DisplayName', 'Glacial', 'SizeData', markerSizeG);
% scatter(CS(:, 1), CS(:, 2), 'k^', 'MarkerFaceColor', [0.616, 0.765, 0.902],'DisplayName', 'Subglacial', 'SizeData', markerSizeG,'HandleVisibility','off');
% scatter(D_SSupra(:, 1), D_SSupra(:, 2), 'k^', 'MarkerFaceColor', [0.616, 0.765, 0.902], ...
%     'DisplayName', 'Supraglacial', 'SizeData', markerSizeG, 'HandleVisibility', 'off');
% scatter(D_SPeri(:, 1), D_SPeri(:, 2), 'k^', 'MarkerFaceColor', [0.616, 0.765, 0.902], ...
%     'DisplayName', 'Periglacial', 'SizeData', markerSizeG, 'HandleVisibility', 'off');
% scatter(D_Den(:, 1), D_Den(:, 2), 'k*', ...
%     'DisplayName', 'Precipitation', 'SizeData', markerSizeG, 'HandleVisibility', 'on');
% scatter(D_Kat(:, 1), D_Kat(:, 2), 'k*', ...
%     'DisplayName', 'Precipitation', 'SizeData', markerSizeG, 'HandleVisibility', 'off');

xlabel('PC1');  
%ylabel('PC2'); 
ax = gca;
ax.FontSize = 10; 
title('Delta Range', 'FontSize',14);
legend('Location', 'NorthWest','FontSize',legendFontSize1);
grid off;
hold off;
set(gcf, 'Position', [100, 100, 1300, 900]);
% Save the combined figure
folderName = 'U:/GoA plots/NewPlots';
fileName = 'Endmember_PCA_precip1_data.svg';
fullFilePath = fullfile(folderName, fileName);
saveas(gcf, fullFilePath, 'svg');