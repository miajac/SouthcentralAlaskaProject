[sFile, sPath] = uigetfile('*.xlsx', 'Select Database File');
sFullPath = fullfile(sPath, sFile);
mC = readtable(sFullPath);
mA = mC(:, 5:end);
fullDataset = table2array(mA);
mLocations = string(table2cell(mC(:, 3)));
%% zscore the data.
[Z, mu, sigma] = zscore(fullDataset);
%% Perform PCA on data.
[coeff,score,latent,tsquared,explained,mu2] = pca(Z);
mLS = score(mLocations == 'LS',:);
mMatanuska = score(mLocations == 'Matanuska',:);
mKnik = score(mLocations == 'Knik',:);
mMoose = score(mLocations == 'Moose',:);
mCastner = score(mLocations == 'Castner',:);
mCanwell = score(mLocations == 'Canwell',:);
mGulkana = score(mLocations == 'Gulkana',:);

figure(1)
scatter(mLS(:, 1), mLS(:, 2), 'ko','MarkerFaceColor','b','DisplayName', 'LS')
hold on 
scatter(mMatanuska(:,1), mMatanuska(:,2), 'ko','MarkerFaceColor','r','DisplayName', 'Mat');
scatter(mKnik(:,1),mKnik(:,2),'ko','MarkerFaceColor','g','DisplayName', 'Knik');
scatter(mMoose(:,1),mMoose(:,2),'ko','MarkerFaceColor','c','DisplayName', 'MC');
scatter(mCastner(:,1),mCastner(:,2),'ko','MarkerFaceColor','m','DisplayName', 'Castner');
scatter(mCanwell(:,1),mCanwell(:,2),'ko','MarkerFaceColor','k','DisplayName', 'Canwell');
scatter(mGulkana(:,1),mGulkana(:,2),'ko','MarkerFaceColor','y','DisplayName', 'Gulkana');
xlabel('PC1');  % Label for x-axis
ylabel('PC2');  % Label for y-axis
title('PCA REE+Y');
legend('Location', 'eastoutside');
grid off;
hold off;
