%% Plots
clear all, close all, clc %use Spatial
[sFile, sPath] = uigetfile('*.xlsx', 'Select Database File');
sFullPath = fullfile(sPath, sFile);

tData = readtable(sFullPath);
%% 
tData=removevars(tData,'ICP_MS');
vColumnLabels = tData.Properties.VariableNames(20:end);
vColumnLabelsArray = cellstr(vColumnLabels);
mFullData = table2array(tData(1:345, 20:end)); %first value is the number of samples (each as a row);second value is the numeric values that will be plotted, must be in an order that can be iterated through (DONT INCLUDE STRINGS)

vSampleDates = datetime(convertStringsToChars(string(table2cell(tData(:, 13))))); % just date not datetime
vSampleLocations = string(table2cell(tData(:, 2))); %labels of each sample (KR4, MR3, etc)
vWatershed = string(table2cell(tData(:,5))); % column that could be used to discriminate between watersheds each location is contained in
vEndmember = string(table2cell(tData(:,6))); % column that defines types of endmembers 
vSampleYear = table2array(tData(:,15)); %column for year
%%
% Knik River 1 elemental data and sample dates
mKR1Data = mFullData(vSampleLocations == "K1", :);
mKR1Dates = vSampleDates(vSampleLocations == "K1");
mKR1Year = vSampleYear(vSampleLocations == "K1");
mKR1_2022 = mKR1Data(mKR1Year == 2022, :);
mKR1_2023 = mKR1Data(mKR1Year == 2023, :);
mKR1_Date22 = mKR1Dates(mKR1Year == 2022);
mKR1_Date23 = mKR1Dates(mKR1Year == 2023);

% Hunter Creek elemental data and sample dates
mHCData = mFullData(vSampleLocations == "HC", :);
mHCDates = vSampleDates(vSampleLocations == "HC");
mHCYear = vSampleYear(vSampleLocations == "HC");
mHC_2022 = mHCData(mHCYear == 2022, :);
mHC_2023 = mHCData(mHCYear == 2023, :);
mHC_Date22 = mHCDates(mHCYear == 2022);
mHC_Date23 = mHCDates(mHCYear == 2023);

% Knik River 2 elemental data and sample dates
mKR2Data = mFullData(vSampleLocations == "K2", :);
mKR2Dates = vSampleDates(vSampleLocations == "K2");
mKR2Year = vSampleYear(vSampleLocations == "K2");
mKR2_2022 = mKR2Data(mKR2Year == 2022, :);
mKR2_2023 = mKR2Data(mKR2Year == 2023, :);
mKR2_Date22 = mKR2Dates(mKR2Year == 2022);
mKR2_Date23 = mKR2Dates(mKR2Year == 2023);

% Knik River 3 elemental data and sample dates
mKR3Data = mFullData(vSampleLocations == "K3", :);
mKR3Dates = vSampleDates(vSampleLocations == "K3");
mKR3Year = vSampleYear(vSampleLocations == "K3");
mKR3_2022 = mKR3Data(mKR3Year == 2022, :);
mKR3_2023 = mKR3Data(mKR3Year == 2023, :);
mKR3_Date22 = mKR3Dates(mKR3Year == 2022);
mKR3_Date23 = mKR3Dates(mKR3Year == 2023);

% Knik River 4 elemental data and sample dates
mKR4Data = mFullData(vSampleLocations == "K4", :);
mKR4Dates = vSampleDates(vSampleLocations == "K4");
mKR4Year = vSampleYear(vSampleLocations == "K4");
mKR4_2022 = mKR4Data(mKR4Year == 2022, :);
mKR4_2023 = mKR4Data(mKR4Year == 2023, :);
mKR4_Date22 = mKR4Dates(mKR4Year == 2022);
mKR4_Date23 = mKR4Dates(mKR4Year == 2023);

% Caribou Creek elemental data and sample dates
mCCData = mFullData(vSampleLocations == "CC", :);
mCCDates = vSampleDates(vSampleLocations == "CC");
mCCYear = vSampleYear(vSampleLocations == "CC");
mCC_2022 = mCCData(mCCYear == 2022, :);
mCC_2023 = mCCData(mCCYear == 2023, :);
mCC_Date22 = mCCDates(mCCYear == 2022);
mCC_Date23 = mCCDates(mCCYear == 2023);

% Matanuska River 0.5 elemental data and sample dates
mMRHalfData = mFullData(vSampleLocations == "M0.5", :);
mMRHalfDates = vSampleDates(vSampleLocations == "M0.5");
mMRHalfYear = vSampleYear(vSampleLocations == "M0.5");
mMRHalf_2022 = mMRHalfData(mMRHalfYear == 2022, :);
mMRHalf_2023 = mMRHalfData(mMRHalfYear == 2023, :);
mMRHalf_Date22 = mMRHalfDates(mMRHalfYear == 2022);
mMRHalf_Date23 = mMRHalfDates(mMRHalfYear == 2023);

% Matanuska River 1 elemental data and sample dates
mMR1Data = mFullData(vSampleLocations == "M1", :);
mMR1Dates = vSampleDates(vSampleLocations == "M1");
mMR1Year = vSampleYear(vSampleLocations == "M1");
mMR1_2022 = mMR1Data(mMR1Year == 2022, :);
mMR1_2023 = mMR1Data(mMR1Year == 2023, :);
mMR1_Date22 = mMR1Dates(mMR1Year == 2022);
mMR1_Date23 = mMR1Dates(mMR1Year == 2023);

% Matanuska River 2 elemental data and sample dates
mMR2Data = mFullData(vSampleLocations == "M2", :);
mMR2Dates = vSampleDates(vSampleLocations == "M2");
mMR2Year = vSampleYear(vSampleLocations == "M2");
mMR2_2022 = mMR2Data(mMR2Year == 2022, :);
mMR2_2023 = mMR2Data(mMR2Year == 2023, :);
mMR2_Date22 = mMR2Dates(mMR2Year == 2022);
mMR2_Date23 = mMR2Dates(mMR2Year == 2023);

% Matanuska River 3 elemental data and sample dates
mMR3Data = mFullData(vSampleLocations == "M3", :);
mMR3Dates = vSampleDates(vSampleLocations == "M3");
mMR3Year = vSampleYear(vSampleLocations == "M3");
mMR3_2022 = mMR3Data(mMR3Year == 2022, :);
mMR3_2023 = mMR3Data(mMR3Year == 2023, :);
mMR3_Date22 = mMR3Dates(mMR3Year == 2022);
mMR3_Date23 = mMR3Dates(mMR3Year == 2023);

% Matanuska River 4 elemental data and sample dates
mMR4Data = mFullData(vSampleLocations == "M4", :);
mMR4Dates = vSampleDates(vSampleLocations == "M4");
mMR4Year = vSampleYear(vSampleLocations == "M4");
mMR4_2022 = mMR4Data(mMR4Year == 2022, :);
mMR4_2023 = mMR4Data(mMR4Year == 2023, :);
mMR4_Date22 = mMR4Dates(mMR4Year == 2022);
mMR4_Date23 = mMR4Dates(mMR4Year == 2023);

% Matanuska River 5 elemental data and sample dates
mMR5Data = mFullData(vSampleLocations == "M5", :);
mMR5Dates = vSampleDates(vSampleLocations == "M5");
mMR5Year = vSampleYear(vSampleLocations == "M5");
mMR5_2022 = mMR5Data(mMR5Year == 2022, :);
mMR5_2023 = mMR5Data(mMR5Year == 2023, :);
mMR5_Date22 = mMR5Dates(mMR5Year == 2022);
mMR5_Date23 = mMR5Dates(mMR5Year == 2023);

% Little Susitna River 1 elemental data and samples dates
mLS1Data = mFullData(vSampleLocations == "LS1", :);
mLS1Dates = vSampleDates(vSampleLocations == "LS1");
mLS1Year = vSampleYear(vSampleLocations == "LS1");
mLS1_2022 = mLS1Data(mLS1Year == 2022, :);
mLS1_2023 = mLS1Data(mLS1Year == 2023, :);
mLS1_Date22 = mLS1Dates(mLS1Year == 2022);
mLS1_Date23 = mLS1Dates(mLS1Year == 2023);

% Little Susitna River 1.5 elemental data and samples dates
mLS15Data = mFullData(vSampleLocations == "LS1.5", :);
mLS15Dates = vSampleDates(vSampleLocations == "LS1.5");
mLS15Year = vSampleYear(vSampleLocations == "LS1.5");
mLS15_2022 = mLS15Data(mLS15Year == 2022, :);
mLS15_2023 = mLS15Data(mLS15Year == 2023, :);
mLS15_Date22 = mLS15Dates(mLS15Year == 2022);
mLS15_Date23 = mLS15Dates(mLS15Year == 2023);

% Little Susitna River 2 elemental data and samples dates
mLS2Data = mFullData(vSampleLocations == "LS2", :);
mLS2Dates = vSampleDates(vSampleLocations == "LS2");
mLS2Year = vSampleYear(vSampleLocations == "LS2");
mLS2_2022 = mLS2Data(mLS2Year == 2022, :);
mLS2_2023 = mLS2Data(mLS2Year == 2023, :);
mLS2_Date22 = mLS2Dates(mLS2Year == 2022);
mLS2_Date23 = mLS2Dates(mLS2Year == 2023);

% Little Susitna River 3 elemental data and samples dates
mLS3Data = mFullData(vSampleLocations == "LS3", :);
mLS3Dates = vSampleDates(vSampleLocations == "LS3");
mLS3Year = vSampleYear(vSampleLocations == "LS3");
mLS3_2022 = mLS3Data(mLS3Year == 2022, :);
mLS3_2023 = mLS3Data(mLS3Year == 2023, :);
mLS3_Date22 = mLS3Dates(mLS3Year == 2022);
mLS3_Date23 = mLS3Dates(mLS3Year == 2023);

% Little Susitna River 4 elemental data and samples dates
mLS4Data = mFullData(vSampleLocations == "LS4", :);
mLS4Dates = vSampleDates(vSampleLocations == "LS4");
mLS4Year = vSampleYear(vSampleLocations == "LS4");
mLS4_2022 = mLS4Data(mLS4Year == 2022, :);
mLS4_2023 = mLS4Data(mLS4Year == 2023, :);
mLS4_Date22 = mLS4Dates(mLS4Year == 2022);
mLS4_Date23 = mLS4Dates(mLS4Year == 2023);

% Moose Creek elemental data and samples dates
mMCData = mFullData(vSampleLocations == "MC", :);
mMCDates = vSampleDates(vSampleLocations == "MC");
mMCYear = vSampleYear(vSampleLocations == "MC");
mMC_2022 = mMCData(mMCYear == 2022, :);
mMC_2023 = mMCData(mMCYear == 2023, :);
mMC_Date22 = mMCDates(mMCYear == 2022);
mMC_Date23 = mMCDates(mMCYear == 2023);

% Gulkana 1 elemental data and samples dates
mG1Data = mFullData(vSampleLocations == "G1", :);
mG1Dates = vSampleDates(vSampleLocations == "G1");
mG1Year = vSampleYear(vSampleLocations == "G1");
mG1_2022 = mG1Data(mG1Year == 2022, :);
mG1_2023 = mG1Data(mG1Year == 2023, :);
mG1_Date22 = mG1Dates(mG1Year == 2022);
mG1_Date23 = mG1Dates(mG1Year == 2023);

% Gulkana 2 elemental data and samples dates
mG2Data = mFullData(vSampleLocations == "G2", :);
mG2Dates = vSampleDates(vSampleLocations == "G2");
mG2Year = vSampleYear(vSampleLocations == "G2");
mG2_2022 = mG2Data(mG2Year == 2022, :);
mG2_2023 = mG2Data(mG2Year == 2023, :);
mG2_Date22 = mG2Dates(mG2Year == 2022);
mG2_Date23 = mG2Dates(mG2Year == 2023);

% Gulkana 3 elemental data and samples dates
mG3Data = mFullData(vSampleLocations == "G3", :);
mG3Dates = vSampleDates(vSampleLocations == "G3");
mG3Year = vSampleYear(vSampleLocations == "G3");
mG3_2022 = mG3Data(mG3Year == 2022, :);
mG3_2023 = mG3Data(mG3Year == 2023, :);
mG3_Date22 = mG3Dates(mG3Year == 2022);
mG3_Date23 = mG3Dates(mG3Year == 2023);

% Castner 1 elemental data and samples dates
mCT1Data = mFullData(vSampleLocations == "CT1", :);
mCT1Dates = vSampleDates(vSampleLocations == "CT1");
mCT1Year = vSampleYear(vSampleLocations == "CT1");
mCT1_2022 = mCT1Data(mCT1Year == 2022, :);
mCT1_2023 = mCT1Data(mCT1Year == 2023, :);
mCT1_Date22 = mCT1Dates(mCT1Year == 2022);
mCT1_Date23 = mCT1Dates(mCT1Year == 2023);

% Castner 2 elemental data and samples dates
mCT2Data = mFullData(vSampleLocations == "CT2", :);
mCT2Dates = vSampleDates(vSampleLocations == "CT2");
mCT2Year = vSampleYear(vSampleLocations == "CT2");
mCT2_2022 = mCT2Data(mCT2Year == 2022, :);
mCT2_2023 = mCT2Data(mCT2Year == 2023, :);
mCT2_Date22 = mCT2Dates(mCT2Year == 2022);
mCT2_Date23 = mCT2Dates(mCT2Year == 2023);

% Canwell 1 elemental data and samples dates
mCW1Data = mFullData(vSampleLocations == "CW1", :);
mCW1Dates = vSampleDates(vSampleLocations == "CW1");
mCW1Year = vSampleYear(vSampleLocations == "CW1");
mCW1_2022 = mCW1Data(mCW1Year == 2022, :);
mCW1_2023 = mCW1Data(mCW1Year == 2023, :);
mCW1_Date22 = mCW1Dates(mCW1Year == 2022);
mCW1_Date23 = mCW1Dates(mCW1Year == 2023);

% Canwell 2 elemental data and samples dates
mCW2Data = mFullData(vSampleLocations == "CW2", :);
mCW2Dates = vSampleDates(vSampleLocations == "CW2");
mCW2Year = vSampleYear(vSampleLocations == "CW2");
mCW2_2022 = mCW2Data(mCW2Year == 2022, :);
mCW2_2023 = mCW2Data(mCW2Year == 2023, :);
mCW2_Date22 = mCW2Dates(mCW2Year == 2022);
mCW2_Date23 = mCW2Dates(mCW2Year == 2023);

%% Define watersheds
vWatershedLocations = string(table2cell(tData(:, 4))); 
mMatData = mFullData(vWatershedLocations == "Matanuska",:);
mMatData1= mFullData(vWatershed == "Matanuska",:);
mMatDates = vSampleDates(vWatershedLocations == "Matanuska");
mMatSites = vSampleLocations(vWatershedLocations == "Matanuska");
mMatYear = vSampleYear(vWatershedLocations == "Matanuska");
vMatEM = vEndmember(vWatershed == "Matanuska");
mMatSpring = mMatData1(vMatEM == "Spring",:);
mMatTrib = mMatData1(vMatEM == "Tributary",:);
mMatMain = mMatData1(vMatEM == "",:);
mMat_2022 = mMatData(mMatYear == 2022, :);
mMat_2023 = mMatData(mMatYear == 2023, :);

mMCData = mFullData(vWatershedLocations == "MC",:);
mMCDates = vSampleDates(vWatershedLocations == "MC");
mMCSites = vSampleLocations(vWatershedLocations == "MC");
mMCYear = vSampleYear(vWatershedLocations == "MC");
mMC_2022 = mMatData(mMatYear == 2022, :);
mMC_2023 = mMatData(mMatYear == 2023, :);
vSampleDates(mMatYear == 2)

mKnikData = mFullData(vWatershedLocations == "Knik",:);
mKnikData1= mFullData(vWatershed == "Knik",:);
mKnikDates = vSampleDates(vWatershedLocations == "Knik");
mKnikSites = vSampleLocations(vWatershedLocations == "Knik");
mKnikYear = vSampleYear(vWatershedLocations == "Knik");
vKnikEM = vEndmember(vWatershed == "Knik");
mKnikSpring = mKnikData1(vKnikEM == "Spring",:);
mKnikTrib = mKnikData1(vKnikEM == "Tributary",:); 
mKnikMain = mKnikData1(vKnikEM == "",:); 
mKnik_2022 = mKnikData(mKnikYear == 2022, :);
mKnik_2023 = mKnikData(mKnikYear == 2023, :);

mLSData = mFullData(vWatershedLocations == "Little Susitna",:);
mLSData1= mFullData(vWatershed == "Little Susitna",:);
mLSDates = vSampleDates(vWatershedLocations == "Little Susitna");
mLSSites = vSampleLocations(vWatershedLocations == "Little Susitna");
mLSYear = vSampleYear(vWatershedLocations == "Little Susitna");
vLSEM = vEndmember(vWatershed == "Little Susitna");
mLSSpring = mLSData1(vLSEM == "Spring",:);
mLSTrib = mLSData1(vLSEM == "Tributary",:);
mLSLake = mLSData1(vLSEM == "Glacial Lake",:);
mLSSub = mLSData1(vLSEM == "Subglacial",:);
mLSSupra = mLSData1(vLSEM == "Supraglacial",:);
mLSPeri = mLSData1(vLSEM == "Periglacial",:);
mLSMain = mLSData1(vLSEM == "",:);
mLS_2022 = mLSData(mLSYear == 2022, :);
mLS_2023 = mLSData(mLSYear == 2023, :);

mGulkanaData = mFullData(vWatershedLocations == "Gulkana",:);
mGulkanaData1= mFullData(vWatershed == "Gulkana",:);
mGulkanaDates = vSampleDates(vWatershedLocations == "Gulkana");
mGulkanaSites = vSampleLocations(vWatershedLocations == "Gulkana");
mGulkanaYear = vSampleYear(vWatershedLocations == "Gulkana");
vGulkanaEM = vEndmember(vWatershed == "Gulkana");
mGulkanaSub = mGulkanaData1(vGulkanaEM == "Subglacial",:);
mGulkanaTrib = mGulkanaData1(vGulkanaEM == "Tributary",:);
mGulkanaMain = mGulkanaData1(vGulkanaEM == "",:);
mGulkanaSupra = mGulkanaData1(vGulkanaEM == "Supraglacial",:);
mGulkana_2022 = mGulkanaData(mGulkanaYear == 2022, :);
mGulkana_2023 = mGulkanaData(mGulkanaYear == 2023, :);

mCanwellData = mFullData(vWatershedLocations == "Canwell",:);
mCanwellDates = vSampleDates(vWatershedLocations == "Canwell");
mCanwellSites = vSampleLocations(vWatershedLocations == "Canwell");
mCanwellYear = vSampleYear(vWatershedLocations == "Canwell");
mCanwell_2022 = mCanwellData(mCanwellYear == 2022, :);
mCanwell_2023 = mCanwellData(mCanwellYear == 2023, :);

mCastnerData = mFullData(vWatershedLocations == "Castner",:);
mCastnerData1= mFullData(vWatershed == "Castner",:);
mCastnerDates1 = vSampleDates(vWatershed == "Castner",:);
mCastnerDates = vSampleDates(vWatershedLocations == "Castner");
mCastnerSites = vSampleLocations(vWatershedLocations == "Castner");
mCastnerYear = vSampleYear(vWatershedLocations == "Castner");
vCastnerEM = vEndmember(vWatershed == "Castner");
mCastnerSub = mCastnerData1(vCastnerEM == "Subglacial",:);
mCastnerMain = mCastnerData1(vCastnerEM == "",:);
mCastner_2022 = mCastnerData(mCastnerYear == 2022, :);
mCastner_2023 = mCastnerData(mCastnerYear == 2023, :);

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

%% Load logarithmic axis v. time plot: USGS sites (RUN ALL ELEMENTS) - done
columnNames = vColumnLabels; %adjusted
folderName = 'U:/GoA plots/NewPlots';
for dPosition = 1:numel(columnNames)
    sInput = columnNames{dPosition};
    sTitle = ['Element ', sInput];
    h = figure('Visible', 'off');
    figure;
    K22_conc = mKR3_2022(:, dPosition);
    K22_dis = mKR3_2022(:, 1);
    y_KR3_22 = K22_conc .* 86.4 .* K22_dis;
    y_KR3_22 = fillmissing(y_KR3_22, 'linear');
    %plot(mKR3_Date22, y_KR3_22, 'o-', 'DisplayName', 'K3','Color','k','MarkerFaceColor',[0.35,0.70,0.90]);
    plot(mKR3_Date22, y_KR3_22, 'o-', 'DisplayName', 'K3','Color',[0.35,0.70,0.90],'MarkerFaceColor',[0.35,0.70,0.90],'MarkerEdgeColor', 'k');
    set(gca, 'YScale', 'log');
    hold on;
    K23_conc = mKR3_2023(:, dPosition);
    K23_dis = mKR3_2023(:, 1);
    y_KR3_23 = K23_conc .* 86.4 .* K23_dis;
    y_KR3_23 = fillmissing(y_KR3_23, 'linear');
    %plot(mKR3_Date23, y_KR3_23, 'o-','Color','k','MarkerFaceColor',[0.35,0.70,0.90],'HandleVisibility', 'off');
    plot(mKR3_Date23, y_KR3_23, 'o-','Color',[0.35,0.70,0.90],'MarkerFaceColor',[0.35,0.70,0.90],'HandleVisibility','off','MarkerEdgeColor','k');
    M22_conc = mMR4_2022(:, dPosition);
    M22_dis = mMR4_2022(:,1);
    y_MR4_2022 = M22_conc.*86.4.*M22_dis;
    y_MR4_2022 = fillmissing(y_MR4_2022, 'linear');
    %plot(mMR4_Date22, y_MR4_2022, 'o-', 'DisplayName', 'M4','Color','k','MarkerFaceColor',[0.95,0.90,0.25]);
    plot(mMR4_Date22, y_MR4_2022, 'o-', 'DisplayName', 'M4','Color',[0.95,0.90,0.25],'MarkerFaceColor',[0.95,0.90,0.25],'MarkerEdgeColor','k');
    M23_conc = mMR4_2023(:, dPosition);
    M23_dis = mMR4_2023(:,1);
    y_MR4_2023 = M23_conc.*86.4.*M23_dis;
    y_MR4_2023 = fillmissing(y_MR4_2023, 'linear');
    %plot(mMR4_Date23, y_MR4_2023, 'o-','Color','k','MarkerFaceColor',[0.95,0.90,0.25],'HandleVisibility', 'off');
    plot(mMR4_Date23, y_MR4_2023, 'o-','Color',[0.95,0.90,0.25],'MarkerFaceColor',[0.95,0.90,0.25],'HandleVisibility', 'off','MarkerEdgeColor','k');
    MC22_conc = mMC_2022(:, dPosition);
    MC22_dis = mMC_2022(:,1);
    y_MC_2022 = MC22_conc.*86.4.* MC22_dis;
    y_MC_2022 = fillmissing(y_MC_2022, 'linear');
    %plot(mMC_Date22, y_MC_2022, 'o-', 'DisplayName', 'MC','Color','k','MarkerFaceColor',[0.90,0.60,0]);
    plot(mMC_Date22, y_MC_2022, 'o-', 'DisplayName', 'MC','Color',[0.90,0.60,0],'MarkerFaceColor',[0.90,0.60,0],'MarkerEdgeColor','k');
    MC23_conc = mMC_2023(:, dPosition);
    MC23_dis = mMC_2023(:,1);
    y_MC_2023 = MC23_conc.*86.4.* MC23_dis;
    y_MC_2023 = fillmissing(y_MC_2023, 'linear');
    %plot(mMC_Date23, y_MC_2023, 'o-','Color','k','MarkerFaceColor',[0.90,0.60,0],'HandleVisibility', 'off');
    plot(mMC_Date23, y_MC_2023, 'o-','Color',[0.90,0.60,0],'MarkerFaceColor',[0.90,0.60,0],'HandleVisibility', 'off','MarkerEdgeColor','k');
    L22_conc = mLS2_2022(:, dPosition);
    L22_dis = mLS2_2022(:,1);
    y_LS2_2022 = L22_conc.*86.4.*L22_dis;
    y_LS2_2022 = fillmissing(y_LS2_2022, 'linear');
    %plot(mLS2_Date22, y_LS2_2022, 'o-', 'DisplayName', 'LS2','Color','k','MarkerFaceColor',[0,0.60,0.50]);
    plot(mLS2_Date22, y_LS2_2022, 'o-', 'DisplayName', 'LS2','Color',[0,0.60,0.50],'MarkerFaceColor',[0,0.60,0.50],'MarkerEdgeColor','k');
    L23_conc = mLS2_2023(:, dPosition);
    L23_dis = mLS2_2023(:,1);
    y_LS2_2023 = L23_conc.*86.4.*L23_dis;
    y_LS2_2023 = fillmissing(y_LS2_2023, 'linear');
    %plot(mLS2_Date23, y_LS2_2023, 'o-','Color','k','MarkerFaceColor',[0,0.60,0.50],'HandleVisibility', 'off');
    plot(mLS2_Date23, y_LS2_2023, 'o-','Color',[0,0.60,0.50],'MarkerFaceColor',[0,0.60,0.50],'HandleVisibility', 'off','MarkerEdgeColor','k');
    %xlabel('Time');  % Label for x-axis
    ylabel('Load (kg/day)');  % Label for y-axis
    title(sInput); 
    legend('Location', 'none');
    legend('Position', [0.5 0.5 0.1 0.1], 'Location', 'none');
    grid on;
    hold off;
    sInputChar = char(sInput);
    formatFileName = ['load_log_t_usgs_no_title_', sInputChar, '.svg'];
    fileName = sprintf(formatFileName);
    fullFilePath = fullfile(folderName, fileName);
    saveas(gcf, fullFilePath, 'svg');
end
%% Load v. time plot: USGS sites (RUN ALL ELEMENTS) - done
columnNames = vColumnLabels; 
folderName = 'U:/GoA plots/NewPlots';
for dPosition = 1:numel(columnNames)
    sInput = columnNames{dPosition};
    sTitle = ['Element ', sInput];
    figure; %start if not iterating
    K22_conc = mKR3_2022(:, dPosition);
    K22_dis = mKR3_2022(:,1);
    y_KR3_2022 = K22_conc.*86.4.*K22_dis;
    y_KR3_2022 = fillmissing(y_KR3_2022, 'linear');
    %plot(mKR3_Date22, y_KR3_2022, 'o-', 'DisplayName', 'K3','Color','k','MarkerFaceColor',[0.35,0.70,0.90]);
    plot(mKR3_Date22, y_KR3_2022, 'o-', 'DisplayName', 'K3','Color',[0.35,0.70,0.90],'MarkerFaceColor',[0.35,0.70,0.90],'MarkerEdgeColor','k');
    hold on;
    K23_conc = mKR3_2023(:, dPosition);
    K23_dis = mKR3_2023(:,1);
    y_KR3_2023 = K23_conc.*86.4.*K23_dis;
    y_KR3_2023 = fillmissing(y_KR3_2023, 'linear');
    %plot(mKR3_Date23, y_KR3_2023, 'o-','Color','k','MarkerFaceColor',[0.35,0.70,0.90],'HandleVisibility', 'off');
    plot(mKR3_Date23, y_KR3_2023, 'o-','Color',[0.35,0.70,0.90],'MarkerFaceColor',[0.35,0.70,0.90],'HandleVisibility', 'off','MarkerEdgeColor','k');
    M22_conc = mMR4_2022(:, dPosition);
    M22_dis = mMR4_2022(:,1);
    y_MR4_2022 = M22_conc.*86.4.*M22_dis;
    y_MR4_2022 = fillmissing(y_MR4_2022, 'linear');
    %plot(mMR4_Date22, y_MR4_2022, 'o-', 'DisplayName', 'M4','Color','k','MarkerFaceColor',[0.95,0.90,0.25]);
    plot(mMR4_Date22, y_MR4_2022, 'o-', 'DisplayName', 'M4','Color',[0.95,0.90,0.25],'MarkerFaceColor',[0.95,0.90,0.25],'MarkerEdgeColor','k');
    M23_conc = mMR4_2023(:, dPosition);
    M23_dis = mMR4_2023(:,1);
    y_MR4_2023 = M23_conc.*86.4.*M23_dis;
    y_MR4_2023 = fillmissing(y_MR4_2023, 'linear');
    %plot(mMR4_Date23, y_MR4_2023, 'o-','Color','k','MarkerFaceColor',[0.95,0.90,0.25],'HandleVisibility', 'off');
    plot(mMR4_Date23, y_MR4_2023, 'o-','Color',[0.95,0.90,0.25],'MarkerFaceColor',[0.95,0.90,0.25],'HandleVisibility', 'off','MarkerEdgeColor','k');
    MC22_conc = mMC_2022(:, dPosition);
    MC22_dis = mMC_2022(:,1);
    y_MC_2022 = MC22_conc.*86.4.*MC22_dis;
    y_MC_2022 = fillmissing(y_MC_2022, 'linear');
    %plot(mMC_Date22, y_MC_2022, 'o-', 'DisplayName', 'MC','Color','k','MarkerFaceColor',[0.90,0.60,0]);
    plot(mMC_Date22, y_MC_2022, 'o-', 'DisplayName', 'MC','Color',[0.90,0.60,0],'MarkerFaceColor',[0.90,0.60,0],'MarkerEdgeColor','k');
    MC23_conc = mMC_2023(:, dPosition);
    MC23_dis = mMC_2023(:,1);
    y_MC_2023 = MC23_conc.*86.4.*MC23_dis;
    y_MC_2023 = fillmissing(y_MC_2023, 'linear');
    %plot(mMC_Date23, y_MC_2023, 'o-','Color','k','MarkerFaceColor',[0.90,0.60,0],'HandleVisibility', 'off');
    plot(mMC_Date23, y_MC_2023, 'o-','Color',[0.90,0.60,0],'MarkerFaceColor',[0.90,0.60,0],'HandleVisibility', 'off','MarkerEdgeColor','k');
    L22_conc = mLS2_2022(:, dPosition);
    L22_dis = mLS2_2022(:,1);
    y_LS2_2022 = L22_conc.*86.4.*L22_dis;
    y_LS2_2022 = fillmissing(y_LS2_2022, 'linear');
    %plot(mLS2_Date22, y_LS2_2022, 'o-', 'DisplayName', 'LS2','Color','k','MarkerFaceColor',[0,0.60,0.50]);
    plot(mLS2_Date22, y_LS2_2022, 'o-', 'DisplayName', 'LS2','Color',[0,0.60,0.50],'MarkerFaceColor',[0,0.60,0.50],'MarkerEdgeColor','k');
    L23_conc = mLS2_2023(:, dPosition);
    L23_dis = mLS2_2023(:,1);
    y_LS2_2023 = L23_conc.*86.4.*L23_dis;
    y_LS2_2023 = fillmissing(y_LS2_2023, 'linear');
    %plot(mLS2_Date23, y_LS2_2023, 'o-','Color','k','MarkerFaceColor',[0,0.60,0.50],'HandleVisibility', 'off');
    plot(mLS2_Date23, y_LS2_2023, 'o-','Color',[0,0.60,0.50],'MarkerFaceColor',[0,0.60,0.50],'HandleVisibility', 'off','MarkerEdgeColor','k');
    xlabel('Time');  % Label for x-axis
    ylabel('Load in kg/day');  % Label for y-axis
    title(sInput);
    legend('Location', 'none');
    legend('Position', [0.5 0.5 0.1 0.1], 'Location', 'none');
    grid on;
    hold off;
    sInputChar = char(sInput);% add if not iterating: folderName = 'U:/GoA plots/';
    formatFileName = "load_t_usgs_%s.jpg";
    fileName = sprintf(formatFileName,sInputChar);
    fullFilePath = fullfile(folderName, fileName);
    saveas(gcf, fullFilePath, 'jpg');%end if not iterating
end
%% Load/ area logarithmic axis v. time plot: USGS sites (RUN ALL ELEMENTS) - done
columnNames = vColumnLabels;
folderName = 'U:/GoA plots/NewPlots';
for dPosition = 1:numel(columnNames)
    sInput = columnNames{dPosition};
    sTitle = ['Element ', sInput];  % Modify the title based on the current sInput
    h = figure('Visible', 'off');
    figure;
    K22_conc = mKR3_2022(:, dPosition);
    K22_dis = mKR3_2022(:, 1);
    y_KR3_2022 = K22_conc .* 86.4 .* K22_dis;
    y_KR3_2022 = y_KR3_2022 ./ 3156.7014; %area 3156.7014 in km^2
    y_KR3_2022 = fillmissing(y_KR3_2022, 'linear');
    %plot(mKR3_Date22, y_KR3_2022, 'o-', 'DisplayName', 'K3','Color','k','MarkerFaceColor',[0.35,0.70,0.90]);
    plot(mKR3_Date22, y_KR3_2022, 'o-', 'DisplayName', 'K3','Color',[0.35,0.70,0.90],'MarkerFaceColor',[0.35,0.70,0.90],'MarkerEdgeColor','k');
    set(gca, 'YScale', 'log');
    hold on;
    K23_conc = mKR3_2023(:, dPosition);
    K23_dis = mKR3_2023(:, 1);
    y_KR3_2023 = K23_conc .* 86.4 .* K23_dis;
    y_KR3_2023 = y_KR3_2023 ./ 3156.7014; %area 3156.7014 in km^2
    y_KR3_2023 = fillmissing(y_KR3_2023, 'linear');
    %plot(mKR3_Date23, y_KR3_2023, 'o-','Color','k','MarkerFaceColor',[0.35,0.70,0.90],'HandleVisibility', 'off');
    plot(mKR3_Date23, y_KR3_2023, 'o-','Color',[0.35,0.70,0.90],'MarkerFaceColor',[0.35,0.70,0.90],'HandleVisibility', 'off','MarkerEdgeColor','k');
    M22_conc = mMR4_2022(:, dPosition);
    M22_dis = mMR4_2022(:,1);
    y_MR4_2022 = M22_conc.*86.4.*M22_dis;
    y_MR4_2022 = y_MR4_2022./ 5329.1646; % area 5329.1646 in km^2
    y_MR4_2022 = fillmissing(y_MR4_2022, 'linear');
    %plot(mMR4_Date22, y_MR4_2022, 'o-', 'DisplayName', 'M4','Color','k','MarkerFaceColor',[0.95,0.90,0.25]);
    plot(mMR4_Date22, y_MR4_2022, 'o-', 'DisplayName', 'M4','Color',[0.95,0.90,0.25],'MarkerFaceColor',[0.95,0.90,0.25],'MarkerEdgeColor','k');
    M23_conc = mMR4_2023(:, dPosition);
    M23_dis = mMR4_2023(:,1);
    y_MR4_2023 = M23_conc.*86.4.* M23_dis;
    y_MR4_2023 = y_MR4_2023./ 5329.1646; % area 5329.1646 in km^2
    y_MR4_2023 = fillmissing(y_MR4_2023, 'linear');
    %plot(mMR4_Date23, y_MR4_2023, 'o-','Color','k','MarkerFaceColor',[0.95,0.90,0.25],'HandleVisibility', 'off');
    plot(mMR4_Date23, y_MR4_2023, 'o-','Color',[0.95,0.90,0.25],'MarkerFaceColor',[0.95,0.90,0.25],'HandleVisibility', 'off','MarkerEdgeColor','k');
    MC22_conc = mMC_2022(:, dPosition);
    MC22_dis = mMC_2022(:,1);
    y_MC_2022 = MC22_conc.*86.4.*MC22_dis;
    y_MC_2022 = y_MC_2022 ./ 125.4186; %area 125.4186 in km^2
    y_MC_2022 = fillmissing(y_MC_2022, 'linear');
    %plot(mMC_Date22, y_MC_2022, 'o-', 'DisplayName', 'MC','Color','k','MarkerFaceColor',[0.90,0.60,0]);
    plot(mMC_Date22, y_MC_2022, 'o-', 'DisplayName', 'MC','Color',[0.90,0.60,0],'MarkerFaceColor',[0.90,0.60,0],'MarkerEdgeColor','k');
    MC23_conc = mMC_2023(:, dPosition);
    MC23_dis = mMC_2023(:,1);
    y_MC_2023 = MC23_conc.*86.4.*MC23_dis;
    y_MC_2023 = y_MC_2023 ./ 125.4186; %area 125.4186 in km^2
    y_MC_2023 = fillmissing(y_MC_2023, 'linear');
    %plot(mMC_Date23, y_MC_2023, 'o-','Color','k','MarkerFaceColor',[0.90,0.60,0],'HandleVisibility', 'off');
    plot(mMC_Date23, y_MC_2023, 'o-','Color',[0.90,0.60,0],'MarkerFaceColor',[0.90,0.60,0],'HandleVisibility', 'off','MarkerEdgeColor','k');
    L22_conc = mLS2_2022(:, dPosition);
    L22_dis = mLS2_2022(:,1);
    y_LS2_2022 = L22_conc.*86.4.*L22_dis;
    y_LS2_2022 = y_LS2_2022 ./ 160.4475; %area 160.4475 in km^2
    y_LS2_2022 = fillmissing(y_LS2_2022, 'linear');
    %plot(mLS2_Date22, y_LS2_2022, 'o-', 'DisplayName', 'LS2','Color','k','MarkerFaceColor',[0,0.60,0.50]);
    plot(mLS2_Date22, y_LS2_2022, 'o-', 'DisplayName', 'LS2','Color',[0,0.60,0.50],'MarkerFaceColor',[0,0.60,0.50],'MarkerEdgeColor','k');
    L23_conc = mLS2_2023(:, dPosition);
    L23_dis = mLS2_2023(:,1);
    y_LS2_2023 = L23_conc.*86.4.*L23_dis;
    y_LS2_2023 = y_LS2_2023 ./ 160.4475; %area 160.4475 in km^2
    y_LS2_2023 = fillmissing(y_LS2_2023, 'linear');
    %plot(mLS2_Date23, y_LS2_2023, 'o-','Color','k','MarkerFaceColor',[0,0.60,0.50],'HandleVisibility', 'off');
    plot(mLS2_Date23, y_LS2_2023, 'o-','Color',[0,0.60,0.50],'MarkerFaceColor',[0,0.60,0.50],'HandleVisibility', 'off','MarkerEdgeColor','k');
    %xlabel('Time');  % Label for x-axis
    ylabel('Load/area (kg/day/km^{2})');  % Label for y-axis
    title(sInput);
    legend('Location', 'none');
    legend('Position', [0.5 0.5 0.1 0.1], 'Location', 'none');
    grid on;
    hold off;
    sInputChar = char(sInput);
    formatFileName = ['load_area_by_time', sInputChar, '.svg'];
    fileName = sprintf(formatFileName);
    fullFilePath = fullfile(folderName, fileName);
    saveas(gcf, fullFilePath, 'svg');
end
%% Load/area v. time plot: USGS sites (RUN ALL ELEMENTS) - done
columnNames = vColumnLabels;
folderName = 'U:/GoA plots/NewPlots';
for dPosition = 1:numel(columnNames)
    sInput = columnNames{dPosition};
    sTitle = ['Element ', sInput];
    figure; %start if not iterating
    K22_conc = mKR3_2022(:, dPosition);
    K22_dis = mKR3_2022(:,1);
    y_KR3_2022 = K22_conc.*86.4.*K22_dis;
    y_KR3_2022 = y_KR3_2022 ./ 3156.7014; %area 3156.7014 in km^2
    y_KR3_2022 = fillmissing(y_KR3_2022, 'linear');
    %plot(mKR3_Date22, y_KR3_2022, 'o-', 'DisplayName', 'K3','Color','k','MarkerFaceColor',[0.35,0.70,0.90]);
    plot(mKR3_Date22, y_KR3_2022, 'o-', 'DisplayName', 'K3','Color',[0.35,0.70,0.90],'MarkerFaceColor',[0.35,0.70,0.90],'MarkerEdgeColor','k');
    hold on;
    K23_conc = mKR3_2023(:, dPosition);
    K23_dis = mKR3_2023(:,1);
    y_KR3_2023 = K23_conc.*86.4.*K23_dis;
    y_KR3_2023 = y_KR3_2023 ./ 3156.7014; %area 3156.7014 in km^2
    y_KR3_2023 = fillmissing(y_KR3_2023, 'linear');
    %plot(mKR3_Date23, y_KR3_2023, 'o-','Color','k','MarkerFaceColor',[0.35,0.70,0.90],'HandleVisibility', 'off');
    plot(mKR3_Date23, y_KR3_2023, 'o-','Color',[0.35,0.70,0.90],'MarkerFaceColor',[0.35,0.70,0.90],'HandleVisibility', 'off','MarkerEdgeColor','k');
    M22_conc = mMR4_2022(:, dPosition);
    M22_dis = mMR4_2022(:,1);
    y_MR4_2022 = M22_conc.*86.4.*M22_dis;
    y_MR4_2022 = y_MR4_2022./ 5329.1646; % area 5329.1646 in km^2
    y_MR4_2022 = fillmissing(y_MR4_2022, 'linear');
    %plot(mMR4_Date22, y_MR4_2022, 'o-', 'DisplayName', 'M4','Color','k','MarkerFaceColor',[0.95,0.90,0.25]);
    plot(mMR4_Date22, y_MR4_2022, 'o-', 'DisplayName', 'M4','Color',[0.95,0.90,0.25],'MarkerFaceColor',[0.95,0.90,0.25],'MarkerEdgeColor','k');
    M23_conc = mMR4_2023(:, dPosition);
    M23_dis = mMR4_2023(:,1);
    y_MR4_2023 = M23_conc.*86.4.*M23_dis;
    y_MR4_2023 = y_MR4_2023./ 5329.1646; % area 5329.1646 in km^2
    y_MR4_2023 = fillmissing(y_MR4_2023, 'linear');
    %plot(mMR4_Date23, y_MR4_2023, 'o-','Color','k','MarkerFaceColor',[0.95,0.90,0.25],'HandleVisibility', 'off');
    plot(mMR4_Date23, y_MR4_2023, 'o-','Color',[0.95,0.90,0.25],'MarkerFaceColor',[0.95,0.90,0.25],'HandleVisibility', 'off','MarkerEdgeColor','k');
    MC22_conc = mMC_2022(:, dPosition);
    MC22_dis = mMC_2022(:,1);
    y_MC_2022 = MC22_conc.*86.4.*MC22_dis;
    y_MC_2022 = y_MC_2022 ./ 125.4186; %area 125.4186 in km^2
    y_MC_2022 = fillmissing(y_MC_2022, 'linear');
    %plot(mMC_Date22, y_MC_2022, 'o-', 'DisplayName', 'MC','Color','k','MarkerFaceColor',[0.90,0.60,0]);
    plot(mMC_Date22, y_MC_2022, 'o-', 'DisplayName', 'MC','Color',[0.90,0.60,0],'MarkerFaceColor',[0.90,0.60,0],'MarkerEdgeColor','k');
    MC23_conc = mMC_2023(:, dPosition);
    MC23_dis = mMC_2023(:,1);
    y_MC_2023 = MC23_conc.*86.4.*MC23_dis;
    y_MC_2023 = y_MC_2023 ./ 125.4186; %area 125.4186 in km^2
    y_MC_2023 = fillmissing(y_MC_2023, 'linear');
    %plot(mMC_Date23, y_MC_2023, 'o-','Color','k','MarkerFaceColor',[0.90,0.60,0],'HandleVisibility', 'off');
    plot(mMC_Date23, y_MC_2023, 'o-','Color',[0.90,0.60,0],'MarkerFaceColor',[0.90,0.60,0],'HandleVisibility', 'off','MarkerEdgeColor','k');
    L22_conc = mLS2_2022(:, dPosition);
    L22_dis = mLS2_2022(:,1);
    y_LS2_2022 = L22_conc.*86.4.*L22_dis;
    y_LS2_2022 = y_LS2_2022 ./ 160.4475; %area 160.4475 in km^2
    y_LS2_2022 = fillmissing(y_LS2_2022, 'linear');
    %plot(mLS2_Date22, y_LS2_2022, 'o-', 'DisplayName', 'LS2','Color','k','MarkerFaceColor',[0,0.60,0.50]);
    plot(mLS2_Date22, y_LS2_2022, 'o-', 'DisplayName', 'LS2','Color',[0,0.60,0.50],'MarkerFaceColor',[0,0.60,0.50],'MarkerEdgeColor','k');
    L23_conc = mLS2_2023(:, dPosition);
    L23_dis = mLS2_2023(:,1);
    y_LS2_2023 = L23_conc.*86.4.*L23_dis;
    y_LS2_2023 = y_LS2_2023 ./ 160.4475; %area 160.4475 in km^2
    y_LS2_2023 = fillmissing(y_LS2_2023, 'linear');
    %plot(mLS2_Date23, y_LS2_2023, 'o-','Color','k','MarkerFaceColor',[0,0.60,0.50],'HandleVisibility', 'off');
    plot(mLS2_Date23, y_LS2_2023, 'o-','Color',[0,0.60,0.50],'MarkerFaceColor',[0,0.60,0.50],'HandleVisibility', 'off','MarkerEdgeColor','k');
    xlabel('Time');  % Label for x-axis
    ylabel('Load in kg/day');  % Label for y-axis
    title(sInput);
    legend('Location', 'none');
    legend('Position', [0.5 0.5 0.1 0.1], 'Location', 'none');
    grid on;
    hold off;
    sInputChar = char(sInput);% add if not iterating: folderName = 'U:/GoA plots/';
    formatFileName = "loadByArea_t_usgs_%s.jpg";
    fileName = sprintf(formatFileName,sInputChar);
    fullFilePath = fullfile(folderName, fileName);
    saveas(gcf, fullFilePath, 'jpg');%end if not iterating
end

%% Conc v discharge - done
columnNames = vColumnLabels;
folderName = 'U:/GoA plots/NewPlots';
for dPosition = 1:numel(columnNames)
    sInput = columnNames{dPosition};
    sTitle = ['Element ', sInput];
    figure;
    y_KR3_2022 = mKR3_2022(:, dPosition);
    x_KR3_2022 = mKR3_2022(:, 1);
    y_KR3_2022 = fillmissing(y_KR3_2022, 'linear');
    %plot(x_KR3_2022, y_KR3_2022, 'o', 'DisplayName', 'K3','Color','k','MarkerFaceColor',[0.35,0.70,0.90]);
    plot(x_KR3_2022, y_KR3_2022, 'o', 'DisplayName', 'K3','Color',[0.35,0.70,0.90],'MarkerFaceColor',[0.35,0.70,0.90],'MarkerEdgeColor','k');
    hold on;
    y_KR3_2023 = mKR3_2023(:, dPosition);
    x_KR3_2023 = mKR3_2023(:, 1);
    y_KR3_2023 = fillmissing(y_KR3_2023, 'linear');
    %plot(x_KR3_2023, y_KR3_2023, 'o','Color','k','MarkerFaceColor',[0.35,0.70,0.90],'HandleVisibility', 'off');
    plot(x_KR3_2023, y_KR3_2023, 'o','Color',[0.35,0.70,0.90],'MarkerFaceColor',[0.35,0.70,0.90],'HandleVisibility', 'off','MarkerEdgeColor','k');
    y_MR4_2022 = mMR4_2022(:, dPosition);
    x_MR4_2022 = mMR4_2022(:, 1);
    y_MR4_2022 = fillmissing(y_MR4_2022, 'linear');
    %plot(x_MR4_2022, y_MR4_2022, 'o', 'DisplayName', 'M4','Color','k','MarkerFaceColor',[0.95,0.90,0.25]);
    plot(x_MR4_2022, y_MR4_2022, 'o', 'DisplayName', 'M4','Color',[0.95,0.90,0.25],'MarkerFaceColor',[0.95,0.90,0.25],'MarkerEdgeColor','k');
    y_MR4_2023 = mMR4_2023(:, dPosition);
    x_MR4_2023 = mMR4_2023(:, 1);
    y_MR4_2023 = fillmissing(y_MR4_2023, 'linear');
    %plot(x_MR4_2023, y_MR4_2023, 'o','Color','k','MarkerFaceColor',[0.95,0.90,0.25],'HandleVisibility', 'off');
    plot(x_MR4_2023, y_MR4_2023, 'o','Color',[0.95,0.90,0.25],'MarkerFaceColor',[0.95,0.90,0.25],'HandleVisibility', 'off','MarkerEdgeColor','k');
    y_MC_2022 = mMC_2022(:, dPosition);
    x_MC_2022 = mMC_2022(:, 1);
    y_MC_2022 = fillmissing(y_MC_2022, 'linear');
    %plot(x_MC_2022, y_MC_2022, 'o', 'DisplayName', 'MC','Color','k','MarkerFaceColor',[0.90,0.60,0]);
    plot(x_MC_2022, y_MC_2022, 'o', 'DisplayName', 'MC','Color',[0.90,0.60,0],'MarkerFaceColor',[0.90,0.60,0],'MarkerEdgeColor','k');
    y_MC_2023 = mMC_2023(:, dPosition);
    x_MC_2023 = mMC_2023(:, 1);
    y_MC_2023 = fillmissing(y_MC_2023, 'linear');
    %plot(x_MC_2023, y_MC_2023, 'o','Color','k','MarkerFaceColor',[0.90,0.60,0],'HandleVisibility', 'off');
    plot(x_MC_2023, y_MC_2023, 'o','Color',[0.90,0.60,0],'MarkerFaceColor',[0.90,0.60,0],'HandleVisibility', 'off','MarkerEdgeColor','k');
    y_LS2_2022 = mLS2_2022(:, dPosition);
    x_LS2_2022 = mLS2_2022(:, 1);
    y_LS2_2022 = fillmissing(y_LS2_2022, 'linear');
    %plot(x_LS2_2022, y_LS2_2022, 'o', 'DisplayName', 'LS2','Color','k','MarkerFaceColor',[0,0.60,0.50]);
    plot(x_LS2_2022, y_LS2_2022, 'o', 'DisplayName', 'LS2','Color',[0,0.60,0.50],'MarkerFaceColor',[0,0.60,0.50],'MarkerEdgeColor','k');
    y_LS2_2023 = mLS2_2023(:, dPosition);
    x_LS2_2023 = mLS2_2023(:, 1);
    y_LS2_2023 = fillmissing(y_LS2_2023, 'linear');
    %plot(x_LS2_2023, y_LS2_2023, 'o','Color','k','MarkerFaceColor',[0,0.60,0.50],'HandleVisibility', 'off');
    plot(x_LS2_2023, y_LS2_2023, 'o','Color',[0,0.60,0.50],'MarkerFaceColor',[0,0.60,0.50],'HandleVisibility', 'off','MarkerEdgeColor','k');
    set(gca, 'XScale', 'log');
    set(gca, 'YScale', 'log');
    xlabel('Discharge (m^{3}/s)');  % Label for x-axis
    ylabel('Concentration (mg/L)');  % Label for y-axis
    title(sInput);
    legend('Location', 'eastoutside');
    grid on;
    hold off;
    sInputChar = char(sInput);
    formatFileName = "conc_dis__log_usgs_%s.svg";
    fileName = sprintf(formatFileName,sInputChar);
    fullFilePath = fullfile(folderName, fileName);
    saveas(gcf, fullFilePath, 'svg');
end
%% Load v time with concentration - done
columnNames = vColumnLabels;
folderName = 'U:/GoA plots/NewPlots';
for dPosition = 1:numel(columnNames)
    sInput = columnNames{dPosition};
    sTitle = ['Element ', sInput];
    figure;
    yyaxis left;
    K22_conc = mKR3_2022(:, dPosition);
    K22_dis = mKR3_2022(:,1);
    y_KR3_2022 = K22_conc.*86.4.*K22_dis;
    y_KR3_2022 = fillmissing(y_KR3_2022, 'linear');
    %plot(mKR3_Date22, y_KR3_2022, 'o-', 'DisplayName', 'K3 Load','Color','k','MarkerFaceColor',[0.35,0.70,0.90]);
    plot(mKR3_Date22, y_KR3_2022, 'o-', 'DisplayName', 'K3 Load','Color',[0.35,0.70,0.90],'MarkerFaceColor',[0.35,0.70,0.90],'MarkerEdgeColor','k');
    hold on;
    K23_conc = mKR3_2023(:, dPosition);
    K23_dis = mKR3_2023(:,1);
    y_KR3_2023 = K23_conc.*86.4.*K23_dis;
    y_KR3_2023 = fillmissing(y_KR3_2023, 'linear');
    %plot(mKR3_Date23, y_KR3_2023, 'o-','Color','k','MarkerFaceColor',[0.35,0.70,0.90],'HandleVisibility', 'off');
    plot(mKR3_Date23, y_KR3_2023, 'o-','Color',[0.35,0.70,0.90],'MarkerFaceColor',[0.35,0.70,0.90],'HandleVisibility', 'off','MarkerEdgeColor','k');
    M22_conc = mMR4_2022(:, dPosition);
    M22_dis = mMR4_2022(:,1);
    y_MR4_2022 = M22_conc.*86.4.*M22_dis;
    y_MR4_2022 = fillmissing(y_MR4_2022, 'linear');
    %plot(mMR4_Date22, y_MR4_2022, 'o-', 'DisplayName', 'M4 Load','Color','k','MarkerFaceColor',[0.95,0.90,0.25]);
    plot(mMR4_Date22, y_MR4_2022, 'o-', 'DisplayName', 'M4 Load','Color',[0.95,0.90,0.25],'MarkerFaceColor',[0.95,0.90,0.25],'MarkerEdgeColor','k');
    M23_conc = mMR4_2023(:, dPosition);
    M23_dis = mMR4_2023(:,1);
    y_MR4_2023 = M23_conc.*86.4.*M23_dis;
    y_MR4_2023 = fillmissing(y_MR4_2023, 'linear');
    %plot(mMR4_Date23, y_MR4_2023, 'o-','Color','k','MarkerFaceColor',[0.95,0.90,0.25],'HandleVisibility', 'off');
    plot(mMR4_Date23, y_MR4_2023, 'o-','Color',[0.95,0.90,0.25],'MarkerFaceColor',[0.95,0.90,0.25],'HandleVisibility', 'off','MarkerEdgeColor','k');
    MC22_conc = mMC_2022(:, dPosition);
    MC22_dis = mMC_2022(:,1);
    y_MC_2022 = MC22_conc.*86.4.*MC22_dis;
    y_MC_2022 = fillmissing(y_MC_2022, 'linear');
    %plot(mMC_Date22, y_MC_2022, 'o-', 'DisplayName', 'MC Load','Color','k','MarkerFaceColor',[0.90,0.60,0]);
    plot(mMC_Date22, y_MC_2022, 'o-', 'DisplayName', 'MC Load','Color',[0.90,0.60,0],'MarkerFaceColor',[0.90,0.60,0],'MarkerEdgeColor','k');
    MC23_conc = mMC_2023(:, dPosition);
    MC23_dis = mMC_2023(:,1);
    y_MC_2023 = MC23_conc.*86.4.*MC23_dis;
    y_MC_2023 = fillmissing(y_MC_2023, 'linear');
    %plot(mMC_Date23, y_MC_2023, 'o-','Color','k','MarkerFaceColor',[0.90,0.60,0],'HandleVisibility', 'off');
    plot(mMC_Date23, y_MC_2023, 'o-','Color',[0.90,0.60,0],'MarkerFaceColor',[0.90,0.60,0],'HandleVisibility', 'off','MarkerEdgeColor','k');
    L22_conc = mLS2_2022(:, dPosition);
    L22_dis = mLS2_2022(:,1);
    y_LS2_2022 = L22_conc.*86.4.*L22_dis;
    y_LS2_2022 = fillmissing(y_LS2_2022, 'linear');
    %plot(mLS2_Date22, y_LS2_2022, 'o-', 'DisplayName', 'LS2 Load','Color','k','MarkerFaceColor',[0,0.60,0.50]);
    plot(mLS2_Date22, y_LS2_2022, 'o-', 'DisplayName', 'LS2 Load','Color',[0,0.60,0.50],'MarkerFaceColor',[0,0.60,0.50],'MarkerEdgeColor','k');
    L23_conc = mLS2_2023(:, dPosition);
    L23_dis = mLS2_2023(:,1);
    y_LS2_2023 = L23_conc.*86.4.*L23_dis;
    y_LS2_2023 = fillmissing(y_LS2_2023, 'linear');
    %plot(mLS2_Date23, y_LS2_2023, 'o-','Color','k','MarkerFaceColor',[0,0.60,0.50],'HandleVisibility', 'off');
    plot(mLS2_Date23, y_LS2_2023, 'o-','Color',[0,0.60,0.50],'MarkerFaceColor',[0,0.60,0.50],'HandleVisibility', 'off','MarkerEdgeColor','k');
    xlabel('Time');  % Label for x-axis
    ylabel('Load in kg/day','Color','k');  % Label for y-axis
    title(sInput);
    legend('Location', 'Best');
    yyaxis left;
    yyaxis('left');
    ax = gca;
    ax.YColor = 'k';
    grid on;
    yyaxis right;
    y_KR3_2022 = mKR3_2022(:, 6);
    y_KR3_2022 = fillmissing(y_KR3_2022, 'linear');
    %plot(mKR3_Date22, y_KR3_2022, 'o--', 'DisplayName', 'K3 Conductivity','Color','k','MarkerFaceColor',[0.35,0.70,0.90]);
    plot(mKR3_Date22, y_KR3_2022, 'o--', 'DisplayName', 'K3 Conductivity','Color',[0.35,0.70,0.90],'MarkerFaceColor',[0.35,0.70,0.90],'MarkerEdgeColor','k');
    hold on;
    y_KR3_2023 = mKR3_2023(:, 6);
    y_KR3_2023 = fillmissing(y_KR3_2023, 'linear');
    %plot(mKR3_Date23, y_KR3_2023, 'o--','Color','k','MarkerFaceColor',[0.35,0.70,0.90],'HandleVisibility', 'off');
    plot(mKR3_Date23, y_KR3_2023, 'o--','Color',[0.35,0.70,0.90],'MarkerFaceColor',[0.35,0.70,0.90],'HandleVisibility', 'off','MarkerEdgeColor','k');
    y_MR4_2022 = mMR4_2022(:, 6);
    y_MR4_2022 = fillmissing(y_MR4_2022, 'linear');
    %plot(mMR4_Date22, y_MR4_2022, 'o--', 'DisplayName', 'M4 Conductivity','Color','k','MarkerFaceColor',[0.95,0.90,0.25]);
    plot(mMR4_Date22, y_MR4_2022, 'o--', 'DisplayName', 'M4 Conductivity','Color',[0.95,0.90,0.25],'MarkerFaceColor',[0.95,0.90,0.25],'MarkerEdgeColor','k');
    y_MR4_2023 = mMR4_2023(:, 6);
    y_MR4_2023 = fillmissing(y_MR4_2023, 'linear');
    %plot(mMR4_Date23, y_MR4_2023, 'o--','Color','k','MarkerFaceColor',[0.95,0.90,0.25],'HandleVisibility', 'off');
    plot(mMR4_Date23, y_MR4_2023, 'o--','Color',[0.95,0.90,0.25],'MarkerFaceColor',[0.95,0.90,0.25],'HandleVisibility', 'off','MarkerEdgeColor','k');
    y_MC_2022 = mMC_2022(:, 6);
    y_MC_2022 = fillmissing(y_MC_2022, 'linear');
    %plot(mMC_Date22, y_MC_2022, 'o--', 'DisplayName', 'MC Conductivity','Color','k','MarkerFaceColor',[0.90,0.60,0]);
    plot(mMC_Date22, y_MC_2022, 'o--', 'DisplayName', 'MC Conductivity','Color',[0.90,0.60,0],'MarkerFaceColor',[0.90,0.60,0],'MarkerEdgeColor','k');
    y_MC_2023 = mMC_2023(:, 6);
    y_MC_2023 = fillmissing(y_MC_2023, 'linear');
    %plot(mMC_Date23, y_MC_2023, 'o--','Color','k','MarkerFaceColor',[0.90,0.60,0],'HandleVisibility', 'off');
    plot(mMC_Date23, y_MC_2023, 'o--','Color',[0.90,0.60,0],'MarkerFaceColor',[0.90,0.60,0],'HandleVisibility', 'off','MarkerEdgeColor','k');
    y_LS2_2022 = mLS2_2022(:, 6);
    y_LS2_2022 = fillmissing(y_LS2_2022, 'linear');
    %plot(mLS2_Date22, y_LS2_2022, 'o--', 'DisplayName', 'LS2 Conductivity','Color','k','MarkerFaceColor',[0,0.60,0.50]);
    plot(mLS2_Date22, y_LS2_2022, 'o--', 'DisplayName', 'LS2 Conductivity','Color',[0,0.60,0.50],'MarkerFaceColor',[0,0.60,0.50],'MarkerEdgeColor','k');
    y_LS2_2023 = mLS2_2023(:, 6);
    y_LS2_2023 = fillmissing(y_LS2_2023, 'linear');
    %plot(mLS2_Date23, y_LS2_2023, 'o--','Color','k','MarkerFaceColor',[0,0.60,0.50],'HandleVisibility', 'off');
    plot(mLS2_Date23, y_LS2_2023, 'o--','Color',[0,0.60,0.50],'MarkerFaceColor',[0,0.60,0.50],'HandleVisibility', 'off','MarkerEdgeColor','k');
    xlabel('Time');  % Label for x-axis
    ylabel('Conductivity in uS/cm','Color','k');  
    yyaxis right;
    yyaxis('right');
    ax = gca;
    ax.YColor = 'k';
    fig = gcf; % Get the current figure
    figPosition = get(fig, 'Position'); % Get current figure position
    figPosition(3) = figPosition(3) * 1.3; % Increase the width by 30%
    set(fig, 'Position', figPosition);
    legend('Location', 'none');
    legend('Position', [0.5 0.5 0.1 0.1], 'Location', 'none');
    grid on;
    hold off;
    sInputChar = char(sInput);
    formatFileName = "load_t_conc_usgs_%s.jpg";
    fileName = sprintf(formatFileName,sInputChar);
    fullFilePath = fullfile(folderName, fileName);
    saveas(gcf, fullFilePath, 'jpg');
end
%% Conc v time with discharge: USGS - done
columnNames = vColumnLabels;
folderName = 'U:/GoA plots/NewPlots';
for dPosition = 1:numel(columnNames)
    sInput = columnNames{dPosition};
    sTitle = ['Element ', sInput];    
    figure;
    set(gcf, 'Position', [600 100 1500 1000]);
    yyaxis left;
    y_KR3_2022 = mKR3_2022(:, dPosition);
    y_KR3_2022 = fillmissing(y_KR3_2022, 'linear');
    %plot(mKR3_Date22, y_KR3_2022, 'o-', 'DisplayName', 'K3 Concentration','Color','k','MarkerFaceColor',[0.35,0.70,0.90]);
    plot(mKR3_Date22, y_KR3_2022, 'o-', 'DisplayName', 'K3 Concentration','Color',[0.35,0.70,0.90],'MarkerFaceColor',[0.35,0.70,0.90],'MarkerEdgeColor','k');
    hold on;
    y_KR3_2023 = mKR3_2023(:, dPosition);
    y_KR3_2023 = fillmissing(y_KR3_2023, 'linear');
    %plot(mKR3_Date23, y_KR3_2023, 'o-','Color','k','MarkerFaceColor',[0.35,0.70,0.90],'HandleVisibility', 'off');
    plot(mKR3_Date23, y_KR3_2023, 'o-','Color',[0.35,0.70,0.90],'MarkerFaceColor',[0.35,0.70,0.90],'HandleVisibility', 'off','MarkerEdgeColor','k');
    y_MR4_2022 = mMR4_2022(:, dPosition);
    y_MR4_2022 = fillmissing(y_MR4_2022, 'linear');
    %plot(mMR4_Date22, y_MR4_2022, 'o-', 'DisplayName', 'M4 Concentration','Color','k','MarkerFaceColor',[0.95,0.90,0.25]);
    plot(mMR4_Date22, y_MR4_2022, 'o-', 'DisplayName', 'M4 Concentration','Color',[0.95,0.90,0.25],'MarkerFaceColor',[0.95,0.90,0.25],'MarkerEdgeColor','k');
    y_MR4_2023 = mMR4_2023(:, dPosition);
    y_MR4_2023 = fillmissing(y_MR4_2023, 'linear');
    %plot(mMR4_Date23, y_MR4_2023, 'o-','Color','k','MarkerFaceColor',[0.95,0.90,0.25],'HandleVisibility', 'off');
    plot(mMR4_Date23, y_MR4_2023, 'o-','Color',[0.95,0.90,0.25],'MarkerFaceColor',[0.95,0.90,0.25],'HandleVisibility', 'off','MarkerEdgeColor','k');
    y_MC_2022 = mMC_2022(:, dPosition);
    y_MC_2022 = fillmissing(y_MC_2022, 'linear');
    %plot(mMC_Date22, y_MC_2022, 'o-', 'DisplayName', 'MC Concentration','Color','k','MarkerFaceColor',[0.90,0.60,0]);
    plot(mMC_Date22, y_MC_2022, 'o-', 'DisplayName', 'MC Concentration','Color',[0.90,0.60,0],'MarkerFaceColor',[0.90,0.60,0],'MarkerEdgeColor','k');
    y_MC_2023 = mMC_2023(:, dPosition);
    y_MC_2023 = fillmissing(y_MC_2023, 'linear');
    %plot(mMC_Date23, y_MC_2023, 'o-','Color','k','MarkerFaceColor',[0.90,0.60,0],'HandleVisibility', 'off');
    plot(mMC_Date23, y_MC_2023, 'o-','Color',[0.90,0.60,0],'MarkerFaceColor',[0.90,0.60,0],'HandleVisibility', 'off','MarkerEdgeColor','k');
    y_LS2_2022 = mLS2_2022(:, dPosition);
    y_LS2_2022 = fillmissing(y_LS2_2022, 'linear');
    %plot(mLS2_Date22, y_LS2_2022, 'o-', 'DisplayName', 'LS2 Concentration','Color','k','MarkerFaceColor',[0,0.60,0.50]);
    plot(mLS2_Date22, y_LS2_2022, 'o-', 'DisplayName', 'LS2 Concentration','Color',[0,0.60,0.50],'MarkerFaceColor',[0,0.60,0.50],'MarkerEdgeColor','k');
    y_LS2_2023 = mLS2_2023(:, dPosition);
    y_LS2_2023 = fillmissing(y_LS2_2023, 'linear');
    %plot(mLS2_Date23, y_LS2_2023, 'o-','Color','k','MarkerFaceColor',[0,0.60,0.50],'HandleVisibility', 'off');
    plot(mLS2_Date23, y_LS2_2023, 'o-','Color',[0,0.60,0.50],'MarkerFaceColor',[0,0.60,0.50],'HandleVisibility', 'off','MarkerEdgeColor','k');
    xlabel('Time');  % Label for x-axis
    ylabel('Concentration in mg/L');  % Label for y-axis
    title(sInput);
    legend('Location', 'Best');
    yyaxis left;
    yyaxis('left');
    ax = gca;
    ax.YColor = 'k';
    grid on;
    yyaxis right;
    y_KR3_2022 = mKR3_2022(:, 1);
    y_KR3_2022 = fillmissing(y_KR3_2022, 'linear');
    %plot(mKR3_Date22, y_KR3_2022, 'o--', 'DisplayName', 'K3 Discharge','Color','k','MarkerFaceColor',[0.35,0.70,0.90]);
    plot(mKR3_Date22, y_KR3_2022, 'o--', 'DisplayName', 'K3 Discharge','Color',[0.35,0.70,0.90],'MarkerFaceColor',[0.35,0.70,0.90],'MarkerEdgeColor','k');
    hold on;
    y_KR3_2023 = mKR3_2023(:, 1);
    y_KR3_2023 = fillmissing(y_KR3_2023, 'linear');
    %plot(mKR3_Date23, y_KR3_2023, 'o--','Color','k','MarkerFaceColor',[0.35,0.70,0.90],'HandleVisibility', 'off');
    plot(mKR3_Date23, y_KR3_2023, 'o--','Color',[0.35,0.70,0.90],'MarkerFaceColor',[0.35,0.70,0.90],'HandleVisibility', 'off','MarkerEdgeColor','k');
    y_MR4_2022 = mMR4_2022(:, 1);
    y_MR4_2022 = fillmissing(y_MR4_2022, 'linear');
    %plot(mMR4_Date22, y_MR4_2022, 'o--', 'DisplayName', 'M4 Discharge','Color','k','MarkerFaceColor',[0.95,0.90,0.25]);
    plot(mMR4_Date22, y_MR4_2022, 'o--', 'DisplayName', 'M4 Discharge','Color',[0.95,0.90,0.25],'MarkerFaceColor',[0.95,0.90,0.25],'MarkerEdgeColor','k');
    y_MR4_2023 = mMR4_2023(:, 1);
    y_MR4_2023 = fillmissing(y_MR4_2023, 'linear');
    %plot(mMR4_Date23, y_MR4_2023, 'o--','Color','k','MarkerFaceColor',[0.95,0.90,0.25],'HandleVisibility', 'off');
    plot(mMR4_Date23, y_MR4_2023, 'o--','Color',[0.95,0.90,0.25],'MarkerFaceColor',[0.95,0.90,0.25],'HandleVisibility', 'off','MarkerEdgeColor','k');
    y_MC_2022 = mMC_2022(:, 1);
    y_MC_2022 = fillmissing(y_MC_2022, 'linear');
    %plot(mMC_Date22, y_MC_2022, 'o--', 'DisplayName', 'MC Discharge','Color','k','MarkerFaceColor',[0.90,0.60,0]);
    plot(mMC_Date22, y_MC_2022, 'o--', 'DisplayName', 'MC Discharge','Color',[0.90,0.60,0],'MarkerFaceColor',[0.90,0.60,0],'MarkerEdgeColor','k');
    y_MC_2023 = mMC_2023(:, 1);
    y_MC_2023 = fillmissing(y_MC_2023, 'linear');
    %plot(mMC_Date23, y_MC_2023, 'o--','Color','k','MarkerFaceColor',[0.90,0.60,0],'HandleVisibility', 'off');
    plot(mMC_Date23, y_MC_2023, 'o--','Color',[0.90,0.60,0],'MarkerFaceColor',[0.90,0.60,0],'HandleVisibility', 'off','MarkerEdgeColor','k');
    y_LS2_2022 = mLS2_2022(:, 1);
    y_LS2_2022 = fillmissing(y_LS2_2022, 'linear');
    %plot(mLS2_Date22, y_LS2_2022, 'o--', 'DisplayName', 'LS2 Discharge','Color','k','MarkerFaceColor',[0,0.60,0.50]);
    plot(mLS2_Date22, y_LS2_2022, 'o--', 'DisplayName', 'LS2 Discharge','Color',[0,0.60,0.50],'MarkerFaceColor',[0,0.60,0.50],'MarkerEdgeColor','k');
    y_LS2_2023 = mLS2_2023(:, 1);
    y_LS2_2023 = fillmissing(y_LS2_2023, 'linear');
    %plot(mLS2_Date23, y_LS2_2023, 'o--','Color','k','MarkerFaceColor',[0,0.60,0.50],'HandleVisibility', 'off');
    plot(mLS2_Date23, y_LS2_2023, 'o--','Color',[0,0.60,0.50],'MarkerFaceColor',[0,0.60,0.50],'HandleVisibility', 'off','MarkerEdgeColor','k');
    xlabel('Time');  % Label for x-axis
    ylabel('Discharge m^3/s','Color','k');  
    yyaxis right;
    yyaxis('right');
    ax = gca;
    ax.YColor = 'k';
    legend('Location', 'none');
    legend('Position', [0.5 0.5 0.1 0.1], 'Location', 'none');
    grid on;
    hold off;
    sInputChar = char(sInput);
    formatFileName = "conc_t_dis_usgs_%s.jpg";
    fileName = sprintf(formatFileName,sInputChar);
    fullFilePath = fullfile(folderName, fileName);
    saveas(gcf, fullFilePath, 'jpg');
end
%% Time-series figure: USGS - done
columnNames = vColumnLabels;
folderName = 'U:/GoA plots/NewPlots';
for dPosition = 1:numel(columnNames)
    sInput = columnNames{dPosition};
    sTitle = ['Element ', sInput];
    figure;
    y_KR3_2022 = mKR3_2022(:, dPosition);
    y_KR3_2022 = fillmissing(y_KR3_2022, 'linear');
    %plot(mKR3_Date22, y_KR3_2022, 'o-', 'DisplayName', 'K3','Color','k','MarkerFaceColor',[0.35,0.70,0.90]);
    plot(mKR3_Date22, y_KR3_2022, 'o-', 'DisplayName', 'K3','Color',[0.35,0.70,0.90],'MarkerFaceColor',[0.35,0.70,0.90],'MarkerEdgeColor','k');
    hold on;
    y_KR3_2023 = mKR3_2023(:, dPosition);
    y_KR3_2023 = fillmissing(y_KR3_2023, 'linear');
    %plot(mKR3_Date23, y_KR3_2023, 'o-','Color','k','MarkerFaceColor',[0.35,0.70,0.90],'HandleVisibility', 'off');
    plot(mKR3_Date23, y_KR3_2023, 'o-','Color',[0.35,0.70,0.90],'MarkerFaceColor',[0.35,0.70,0.90],'HandleVisibility', 'off','MarkerEdgeColor','k');
    y_MR4_2022 = mMR4_2022(:, dPosition);
    y_MR4_2022 = fillmissing(y_MR4_2022, 'linear');
    %plot(mMR4_Date22, y_MR4_2022, 'o-', 'DisplayName', 'M4','Color','k','MarkerFaceColor',[0.95,0.90,0.25]);
    plot(mMR4_Date22, y_MR4_2022, 'o-', 'DisplayName', 'M4','Color',[0.95,0.90,0.25],'MarkerFaceColor',[0.95,0.90,0.25],'MarkerEdgeColor','k');
    y_MR4_2023 = mMR4_2023(:, dPosition);
    y_MR4_2023 = fillmissing(y_MR4_2023, 'linear');
    %plot(mMR4_Date23, y_MR4_2023, 'o-','Color','k','MarkerFaceColor',[0.95,0.90,0.25],'HandleVisibility', 'off');
    plot(mMR4_Date23, y_MR4_2023, 'o-','Color',[0.95,0.90,0.25],'MarkerFaceColor',[0.95,0.90,0.25],'HandleVisibility', 'off','MarkerEdgeColor','k');
    y_MC_2022 = mMC_2022(:, dPosition);
    y_MC_2022 = fillmissing(y_MC_2022, 'linear');
    %plot(mMC_Date22, y_MC_2022, 'o-', 'DisplayName', 'MC','Color','k','MarkerFaceColor',[0.90,0.60,0]);
    plot(mMC_Date22, y_MC_2022, 'o-', 'DisplayName', 'MC','Color',[0.90,0.60,0],'MarkerFaceColor',[0.90,0.60,0],'MarkerEdgeColor','k');
    y_MC_2023 = mMC_2023(:, dPosition);
    y_MC_2023 = fillmissing(y_MC_2023, 'linear');
    %plot(mMC_Date23, y_MC_2023, 'o-','Color','k','MarkerFaceColor',[0.90,0.60,0],'HandleVisibility', 'off');
    plot(mMC_Date23, y_MC_2023, 'o-','Color',[0.90,0.60,0],'MarkerFaceColor',[0.90,0.60,0],'HandleVisibility', 'off','MarkerEdgeColor','k');
    y_LS2_2022 = mLS2_2022(:, dPosition);
    y_LS2_2022 = fillmissing(y_LS2_2022, 'linear');
    %plot(mLS2_Date22, y_LS2_2022, 'o-', 'DisplayName', 'LS2','Color','k','MarkerFaceColor',[0,0.60,0.50]);
    plot(mLS2_Date22, y_LS2_2022, 'o-', 'DisplayName', 'LS2','Color',[0,0.60,0.50],'MarkerFaceColor',[0,0.60,0.50],'MarkerEdgeColor','k');
    y_LS2_2023 = mLS2_2023(:, dPosition);
    y_LS2_2023 = fillmissing(y_LS2_2023, 'linear');
    %plot(mLS2_Date23, y_LS2_2023, 'o-','Color','k','MarkerFaceColor',[0,0.60,0.50],'HandleVisibility', 'off');
    plot(mLS2_Date23, y_LS2_2023, 'o-','Color',[0,0.60,0.50],'MarkerFaceColor',[0,0.60,0.50],'HandleVisibility', 'off','MarkerEdgeColor','k');
    xlabel('Time');  % Label for x-axis
    ylabel('Concentration in mg/L');  % Label for y-axis
    title(sInput);
    legend('Location', 'none');
    legend('Position', [0.5 0.5 0.1 0.1], 'Location', 'none');
    grid on;
    hold off;
    sInputChar = char(sInput);
    formatFileName = "time_series_usgs_no_title_%s.svg";
    fileName = sprintf(formatFileName,sInputChar);
    fullFilePath = fullfile(folderName, fileName);
    saveas(gcf, fullFilePath, 'svg');
end
%% Time-series figure: Matanuska Watershed - done
columnNames = vColumnLabels;
folderName = 'U:/GoA plots/NewPlots';
for dPosition = 1:numel(columnNames)
    sInput = columnNames{dPosition};
    sTitle = ['Element ', sInput];
    figure;
    y_MR1_2022 = mMR1_2022(:, dPosition);
    y_MR1_2022 = fillmissing(y_MR1_2022, 'linear');
    %plot(mMR1_Date22, y_MR1_2022, 'o-', 'DisplayName', 'R1','Color','k','MarkerFaceColor',[0.8,0.4,0]);
    plot(mMR1_Date22, y_MR1_2022, 'o-', 'DisplayName', 'M1','Color',[0.8,0.4,0],'MarkerFaceColor',[0.8,0.4,0],'MarkerEdgeColor','k');
    hold on;
    y_MR1_2023 = mMR1_2023(:, dPosition);
    y_MR1_2023 = fillmissing(y_MR1_2023, 'linear');
    %plot(mMR1_Date23, y_MR1_2023, 'o-','Color','k','MarkerFaceColor',[0.8,0.4,0],'HandleVisibility', 'off');
    plot(mMR1_Date23, y_MR1_2023, 'o-','Color',[0.8,0.4,0],'MarkerFaceColor',[0.8,0.4,0],'HandleVisibility', 'off','MarkerEdgeColor','k');
    y_MR2_2022 = mMR2_2022(:, dPosition);
    y_MR2_2022 = fillmissing(y_MR2_2022, 'linear');
    %plot(mMR2_Date22, y_MR2_2022, 'o-', 'DisplayName', 'M2','Color','k','MarkerFaceColor',[0.8,0.6,0.7]);
    plot(mMR2_Date22, y_MR2_2022, 'o-', 'DisplayName', 'M2','Color',[0.8,0.6,0.7],'MarkerFaceColor',[0.8,0.6,0.7],'MarkerEdgeColor','k');
    y_MR2_2023 = mMR2_2023(:, dPosition);
    y_MR2_2023 = fillmissing(y_MR2_2023, 'linear');
    %plot(mMR2_Date23, y_MR2_2023, 'o-','Color','k','MarkerFaceColor',[0.8,0.6,0.7],'HandleVisibility', 'off');
    plot(mMR2_Date23, y_MR2_2023, 'o-','Color',[0.8,0.6,0.7],'MarkerFaceColor',[0.8,0.6,0.7],'HandleVisibility', 'off','MarkerEdgeColor','k');
    y_MR3_2022 = mMR3_2022(:, dPosition);
    y_MR3_2022 = fillmissing(y_MR3_2022, 'linear');
    %plot(mMR3_Date22, y_MR3_2022, 'o-', 'DisplayName', 'M3','Color','k','MarkerFaceColor',[0.95,0.90,0.25]);
    plot(mMR3_Date22, y_MR3_2022, 'o-', 'DisplayName', 'M3','Color',[0.95,0.90,0.25],'MarkerFaceColor',[0.95,0.90,0.25],'MarkerEdgeColor','k');
    y_MR3_2023 = mMR3_2023(:, dPosition);
    y_MR3_2023 = fillmissing(y_MR3_2023, 'linear');
    %plot(mMR3_Date23, y_MR3_2023, 'o-','Color','k','MarkerFaceColor',[0.95,0.90,0.25],'HandleVisibility', 'off');
    plot(mMR3_Date23, y_MR3_2023, 'o-','Color',[0.95,0.90,0.25],'MarkerFaceColor',[0.95,0.90,0.25],'HandleVisibility', 'off','MarkerEdgeColor','k');
    y_MR4_2022 = mMR4_2022(:, dPosition);
    y_MR4_2022 = fillmissing(y_MR4_2022, 'linear');
    %plot(mMR4_Date22, y_MR4_2022, 'o-', 'DisplayName', 'M4','Color','k','MarkerFaceColor',[0,0.6,0.5]);
    plot(mMR4_Date22, y_MR4_2022, 'o-', 'DisplayName', 'M4','Color',[0,0.6,0.5],'MarkerFaceColor',[0,0.6,0.5],'MarkerEdgeColor','k');
    y_MR4_2023 = mMR4_2023(:, dPosition);
    y_MR4_2023 = fillmissing(y_MR4_2023, 'linear');
    %plot(mMR4_Date23, y_MR4_2023, 'o-','Color','k','MarkerFaceColor',[0,0.6,0.5],'HandleVisibility', 'off');
    plot(mMR4_Date23, y_MR4_2023, 'o-','Color',[0,0.6,0.5],'MarkerFaceColor',[0,0.6,0.5],'HandleVisibility', 'off','MarkerEdgeColor','k');
    y_MR5_2022 = mMR5_2022(:, dPosition);
    y_MR5_2022 = fillmissing(y_MR5_2022, 'linear');
    %plot(mMR5_Date22, y_MR5_2022, 'o-', 'DisplayName', 'M5','Color','k','MarkerFaceColor',[0.35,0.70,0.9]);
    plot(mMR5_Date22, y_MR5_2022, 'o-', 'DisplayName', 'M5','Color',[0.35,0.70,0.9],'MarkerFaceColor',[0.35,0.70,0.9],'MarkerEdgeColor','k');
    y_MR5_2023 = mMR5_2023(:, dPosition);
    y_MR5_2023 = fillmissing(y_MR5_2023, 'linear');
    %plot(mMR5_Date23, y_MR5_2023, 'o-','Color','k','MarkerFaceColor',[0.35,0.70,0.9],'HandleVisibility', 'off');
    plot(mMR5_Date23, y_MR5_2023, 'o-','Color',[0.35,0.70,0.9],'MarkerFaceColor',[0.35,0.70,0.9],'HandleVisibility', 'off','MarkerEdgeColor','k');
    xlabel('Time');  % Label for x-axis
    ylabel('Concentration in mg/L');  % Label for y-axis
    title(sInput);
    legend('Location', 'none');
    legend('Position', [0.5 0.5 0.1 0.1], 'Location', 'none');
    grid on;
    hold off;
    sInputChar = char(sInput);
    formatFileName = "time_series_mat_%s.jpg";
    fileName = sprintf(formatFileName,sInputChar);
    fullFilePath = fullfile(folderName, fileName);
    saveas(gcf, fullFilePath, 'jpg');
end
%% Time-series figure: Knik Watershed - done
columnNames = vColumnLabels; 
folderName = 'U:/GoA plots/NewPlots';
for dPosition = 1:numel(columnNames)
    sInput = columnNames{dPosition};
    sTitle = ['Element ', sInput];
    figure;
    y_KR1_2022 = mKR1_2022(:, dPosition);
    y_KR1_2022 = fillmissing(y_KR1_2022, 'linear');
    %plot(mKR1_Date22, y_KR1_2022, 'o-','DisplayName', 'K1','Color','k','MarkerFaceColor',[0.95,0.90,0.25]);
    plot(mKR1_Date22, y_KR1_2022, 'o-','DisplayName', 'K1','Color',[0.95,0.90,0.25],'MarkerFaceColor',[0.95,0.90,0.25],'MarkerEdgeColor','k');
    hold on;
    y_KR1_2023 = mKR1_2023(:, dPosition);
    y_KR1_2023 = fillmissing(y_KR1_2023, 'linear');
    %plot(mKR1_Date23, y_KR1_2023, 'o-','Color','k','MarkerFaceColor',[0.95,0.90,0.25],'HandleVisibility', 'off');
    plot(mKR1_Date23, y_KR1_2023, 'o-','Color',[0.95,0.90,0.25],'MarkerFaceColor',[0.95,0.90,0.25],'HandleVisibility', 'off','MarkerEdgeColor','k');
    y_KR2_2022 = mKR2_2022(:, dPosition);
    y_KR2_2022 = fillmissing(y_KR2_2022, 'linear');
    %plot(mKR2_Date22, y_KR2_2022, 'o-', 'DisplayName', 'K2','Color','k','MarkerFaceColor',[0.8,0.6,0.7]);
    plot(mKR2_Date22, y_KR2_2022, 'o-', 'DisplayName', 'K2','Color',[0.8,0.6,0.7],'MarkerFaceColor',[0.8,0.6,0.7],'MarkerEdgeColor','k');
    y_KR2_2023 = mKR2_2023(:, dPosition);
    y_KR2_2023 = fillmissing(y_KR2_2023, 'linear');
    %plot(mKR2_Date23, y_KR2_2023, 'o-','Color','k','MarkerFaceColor',[0.8,0.6,0.7],'HandleVisibility', 'off');
    plot(mKR2_Date23, y_KR2_2023, 'o-','Color',[0.8,0.6,0.7],'MarkerFaceColor',[0.8,0.6,0.7],'HandleVisibility', 'off','MarkerEdgeColor','k');
    y_KR3_2022 = mKR3_2022(:, dPosition);
    y_KR3_2022 = fillmissing(y_KR3_2022, 'linear');
    %plot(mKR3_Date22, y_KR3_2022, 'o-', 'DisplayName', 'K3','Color','k','MarkerFaceColor',[0.8,0.4,0]);
    plot(mKR3_Date22, y_KR3_2022, 'o-', 'DisplayName', 'K3','Color',[0.8,0.4,0],'MarkerFaceColor',[0.8,0.4,0],'MarkerEdgeColor','k');
    y_KR3_2023 = mKR3_2023(:, dPosition);
    y_KR3_2023 = fillmissing(y_KR3_2023, 'linear');
    %plot(mKR3_Date23, y_KR3_2023, 'o-','Color','k','MarkerFaceColor',[0.8,0.4,0],'HandleVisibility', 'off');
    plot(mKR3_Date23, y_KR3_2023, 'o-','Color',[0.8,0.4,0],'MarkerFaceColor',[0.8,0.4,0],'HandleVisibility', 'off','MarkerEdgeColor','k');
    y_KR4_2022 = mKR4_2022(:, dPosition);
    y_KR4_2022 = fillmissing(y_KR4_2022, 'linear');
    %plot(mKR4_Date22, y_KR4_2022, 'o-', 'DisplayName', 'K4','Color','k','MarkerFaceColor',[0,0.6,0.5]);
    plot(mKR4_Date22, y_KR4_2022, 'o-', 'DisplayName', 'K4','Color',[0,0.6,0.5],'MarkerFaceColor',[0,0.6,0.5],'MarkerEdgeColor','k');
    y_KR4_2023 = mKR4_2023(:, dPosition);
    y_KR4_2023 = fillmissing(y_KR4_2023, 'linear');
    %plot(mKR4_Date23, y_KR4_2023, 'o-','Color','k','MarkerFaceColor',[0,0.6,0.5],'HandleVisibility', 'off');
    plot(mKR4_Date23, y_KR4_2023, 'o-','Color',[0,0.6,0.5],'MarkerFaceColor',[0,0.6,0.5],'HandleVisibility', 'off','MarkerEdgeColor','k');
    xlabel('Time');  % Label for x-axis
    ylabel('Concentration in mg/L');  % Label for y-axis
    title(sInput);
    legend('Location', 'none');
    legend('Position', [0.5 0.5 0.1 0.1], 'Location', 'none');
    grid on;
    hold off;
    sInputChar = char(sInput);
    formatFileName = "time_series_knik_%s.jpg";
    fileName = sprintf(formatFileName,sInputChar);
    fullFilePath = fullfile(folderName, fileName);
    saveas(gcf, fullFilePath, 'jpg');
end
%% Time-series figure: Little Susitna Watershed - done
columnNames = vColumnLabels;
folderName = 'U:/GoA plots/NewPlots';
for dPosition = 1:numel(columnNames)
    sInput = columnNames{dPosition};
    sTitle = ['Element ', sInput];
    figure;
    y_LS1_2022 = mLS1_2022(:, dPosition);
    y_LS1_2022 = fillmissing(y_LS1_2022, 'linear');
    %plot(mLS1_Date22, y_LS1_2022, 'o-','DisplayName', 'LS1','Color','k','MarkerFaceColor',[0.8,0.4,0]);
    plot(mLS1_Date22, y_LS1_2022, 'o-','DisplayName', 'LS1','Color',[0.8,0.4,0],'MarkerFaceColor',[0.8,0.4,0],'MarkerEdgeColor','k');
    hold on;
    y_LS1_2023 = mLS1_2023(:, dPosition);
    y_LS1_2023 = fillmissing(y_LS1_2023, 'linear');
    %plot(mLS1_Date23, y_LS1_2023, 'o-','Color','k','MarkerFaceColor',[0.8,0.4,0],'HandleVisibility', 'off');
    plot(mLS1_Date23, y_LS1_2023, 'o-','Color',[0.8,0.4,0],'MarkerFaceColor',[0.8,0.4,0],'HandleVisibility', 'off','MarkerEdgeColor','k');
    y_LS15_2022 = mLS15_2022(:, dPosition);
    y_LS15_2022 = fillmissing(y_LS15_2022, 'linear');
    %plot(mLS15_Date22, y_LS15_2022, 'o-', 'DisplayName', 'LS1.5','Color','k','MarkerFaceColor',[0.8,0.6,0.7]);
    plot(mLS15_Date22, y_LS15_2022, 'o-', 'DisplayName', 'LS1.5','Color',[0.8,0.6,0.7],'MarkerFaceColor',[0.8,0.6,0.7],'MarkerEdgeColor','k');
    y_LS15_2023 = mLS15_2023(:, dPosition);
    y_LS15_2023 = fillmissing(y_LS15_2023, 'linear');
    %plot(mLS15_Date23, y_LS15_2023, 'o-','Color','k','MarkerFaceColor',[0.8,0.6,0.7],'HandleVisibility', 'off');
    plot(mLS15_Date23, y_LS15_2023, 'o-','Color',[0.8,0.6,0.7],'MarkerFaceColor',[0.8,0.6,0.7],'HandleVisibility', 'off','MarkerEdgeColor','k');
    y_LS2_2022 = mLS2_2022(:, dPosition);
    y_LS2_2022 = fillmissing(y_LS2_2022, 'linear');
    %plot(mLS2_Date22, y_LS2_2022, 'o-', 'DisplayName', 'LR2','Color','k','MarkerFaceColor',[0.95,0.90,0.25]);
    plot(mLS2_Date22, y_LS2_2022, 'o-', 'DisplayName', 'LR2','Color',[0.95,0.90,0.25],'MarkerFaceColor',[0.95,0.90,0.25],'MarkerEdgeColor','k');
    y_LS2_2023 = mLS2_2023(:, dPosition);
    y_LS2_2023 = fillmissing(y_LS2_2023, 'linear');
    %plot(mLS2_Date23, y_LS2_2023, 'o-','Color','k','MarkerFaceColor',[0.95,0.90,0.25],'HandleVisibility', 'off');
    plot(mLS2_Date23, y_LS2_2023, 'o-','Color',[0.95,0.90,0.25],'MarkerFaceColor',[0.95,0.90,0.25],'HandleVisibility', 'off','MarkerEdgeColor','k');
    y_LS3_2022 = mLS3_2022(:, dPosition);
    y_LS3_2022 = fillmissing(y_LS3_2022, 'linear');
    %plot(mLS3_Date22, y_LS3_2022, 'o-', 'DisplayName', 'LS3','Color','k','MarkerFaceColor',[0,0.6,0.5]);
    plot(mLS3_Date22, y_LS3_2022, 'o-', 'DisplayName', 'LS3','Color',[0,0.6,0.5],'MarkerFaceColor',[0,0.6,0.5],'MarkerEdgeColor','k');
    y_LS3_2023 = mLS3_2023(:, dPosition);
    y_LS3_2023 = fillmissing(y_LS3_2023, 'linear');
    %plot(mLS3_Date23, y_LS3_2023, 'o-','Color','k','MarkerFaceColor',[0,0.6,0.5],'HandleVisibility', 'off');
    plot(mLS3_Date23, y_LS3_2023, 'o-','Color',[0,0.6,0.5],'MarkerFaceColor',[0,0.6,0.5],'HandleVisibility', 'off','MarkerEdgeColor','k');
    y_LS4_2022 = mLS4_2022(:, dPosition);
    y_LS4_2022 = fillmissing(y_LS4_2022, 'linear');
    %plot(mLS4_Date22, y_LS4_2022, 'o-', 'DisplayName', 'LS4','Color','k','MarkerFaceColor',[0.35,0.70,0.9]);
    plot(mLS4_Date22, y_LS4_2022, 'o-', 'DisplayName', 'LS4','Color',[0.35,0.70,0.9],'MarkerFaceColor',[0.35,0.70,0.9],'MarkerEdgeColor','k');
    y_LS4_2023 = mLS4_2023(:, dPosition);
    y_LS4_2023 = fillmissing(y_LS4_2023, 'linear');
    %plot(mLS4_Date23, y_LS4_2023, 'o-','Color','k','MarkerFaceColor',[0.35,0.70,0.9],'HandleVisibility', 'off');
    plot(mLS4_Date23, y_LS4_2023, 'o-','Color',[0.35,0.70,0.9],'MarkerFaceColor',[0.35,0.70,0.9],'HandleVisibility', 'off','MarkerEdgeColor','k');
    xlabel('Time');  % Label for x-axis
    ylabel('Concentration in mg/L');  % Label for y-axis
    title(sInput);
    legend('Location', 'none');
    legend('Position', [0.5 0.5 0.1 0.1], 'Location', 'none');
    grid on;
    hold off;
    sInputChar = char(sInput);
    formatFileName = "time_series_LS_%s.jpg";
    fileName = sprintf(formatFileName,sInputChar);
    fullFilePath = fullfile(folderName, fileName);
    saveas(gcf, fullFilePath, 'jpg');
end
%% Time-series figure: Canwell Watershed -done
columnNames = vColumnLabels;
folderName = 'U:/GoA plots/NewPlots';
for dPosition = 1:numel(columnNames)
    sInput = columnNames{dPosition};
    sTitle = ['Element ', sInput];
    figure;
    y_CW1_2022 = mCW1_2022(:, dPosition);
    y_CW1_2022 = fillmissing(y_CW1_2022, 'linear');
    %plot(mCW1_Date22, y_CW1_2022, 'o-','DisplayName', 'Canwell 1','Color','k','MarkerFaceColor',[0.8,0.4,0]);
    plot(mCW1_Date22, y_CW1_2022, 'o-','DisplayName', 'Canwell 1','Color',[0.8,0.4,0],'MarkerFaceColor',[0.8,0.4,0],'MarkerEdgeColor','k');
    hold on;
    y_CW1_2023 = mCW1_2023(:, dPosition);
    y_CW1_2023 = fillmissing(y_CW1_2023, 'linear');
    %plot(mCW1Dates, y_CW1_2023, 'o-','Color','k','MarkerFaceColor',[0.8,0.4,0],'HandleVisibility', 'off');
    plot(mCW1Dates, y_CW1_2023, 'o-','Color',[0.8,0.4,0],'MarkerFaceColor',[0.8,0.4,0],'HandleVisibility', 'off','MarkerEdgeColor','k');
    y_CW2_2022 = mCW2_2022(:, dPosition);
    y_CW2_2022 = fillmissing(y_CW2_2022, 'linear');
    %plot(mCW2_Date22, y_CW2_2022, 'o-', 'DisplayName', 'Canwell 2','Color','k','MarkerFaceColor',[0,0.6,0.5]);
    plot(mCW2_Date22, y_CW2_2022, 'o-', 'DisplayName', 'Canwell 2','Color',[0,0.6,0.5],'MarkerFaceColor',[0,0.6,0.5],'MarkerEdgeColor','k');
    y_CW2_2023 = mCW2_2023(:, dPosition);
    y_CW2_2023 = fillmissing(y_CW2_2023, 'linear');
    %plot(mCW2_Date23, y_CW2_2023, 'o-','Color','k','MarkerFaceColor',[0,0.6,0.5],'HandleVisibility', 'off');
    plot(mCW2_Date23, y_CW2_2023, 'o-','Color',[0,0.6,0.5],'MarkerFaceColor',[0,0.6,0.5],'HandleVisibility', 'off','MarkerEdgeColor','k');
    xlabel('Time');  % Label for x-axis
    ylabel('Concentration in mg/L');  % Label for y-axis
    title(sInput);
    legend('Location', 'none');
    legend('Position', [0.5 0.5 0.1 0.1], 'Location', 'none');
    grid on;
    hold off;
    sInputChar = char(sInput);
    formatFileName = "time_series_canwell_%s.jpg";
    fileName = sprintf(formatFileName,sInputChar);
    fullFilePath = fullfile(folderName, fileName);
    saveas(gcf, fullFilePath, 'jpg');
end
%% Time-series figure: Castner Watershed - done
columnNames = vColumnLabels;
folderName = 'U:/GoA plots/NewPlots';
for dPosition = 1:numel(columnNames)
    sInput = columnNames{dPosition};
    sTitle = ['Element ', sInput];
    figure;
    y_CT1_2022 = mCT1_2022(:, dPosition);
    y_CT1_2022 = fillmissing(y_CT1_2022, 'linear');
    %plot(mCT1_Date22, y_CT1_2022, 'o-','DisplayName', 'Castner 1','Color','k','MarkerFaceColor',[0.8,0.4,0]);
    plot(mCT1_Date22, y_CT1_2022, 'o-','DisplayName', 'Castner 1','Color',[0.8,0.4,0],'MarkerFaceColor',[0.8,0.4,0],'MarkerEdgeColor','k');
    hold on;
    y_CT1_2023 = mCT1_2023(:, dPosition);
    y_CT1_2023 = fillmissing(y_CT1_2023, 'linear');
    %plot(mCT1_Date23, y_CT1_2023, 'o-','Color','k','MarkerFaceColor',[0.8,0.4,0],'HandleVisibility', 'off');
    plot(mCT1_Date23, y_CT1_2023, 'o-','Color',[0.8,0.4,0],'MarkerFaceColor',[0.8,0.4,0],'HandleVisibility', 'off','MarkerEdgeColor','k');
    y_CT2_2022 = mCT2_2022(:, dPosition);
    y_CT2_2022 = fillmissing(y_CT2_2022, 'linear');
    %plot(mCT2_Date22, y_CT2_2022, 'o-', 'DisplayName', 'Castner 2','Color','k','MarkerFaceColor',[0,0.6,0.5]);
    plot(mCT2_Date22, y_CT2_2022, 'o-', 'DisplayName', 'Castner 2','Color',[0,0.6,0.5],'MarkerFaceColor',[0,0.6,0.5],'MarkerEdgeColor','k');
    y_CT2_2023 = mCT2_2023(:, dPosition);
    y_CT2_2023 = fillmissing(y_CT2_2023, 'linear');
    %plot(mCT2_Date23, y_CT2_2023, 'o-','Color','k','MarkerFaceColor',[0,0.6,0.5],'HandleVisibility', 'off');
    plot(mCT2_Date23, y_CT2_2023, 'o-','Color',[0,0.6,0.5],'MarkerFaceColor',[0,0.6,0.5],'HandleVisibility', 'off','MarkerEdgeColor','k');
    xlabel('Time');  
    ylabel('Concentration in mg/L');  
    title(sInput);
    legend('Location', 'none');
    legend('Position', [0.5 0.5 0.1 0.1], 'Location', 'none');
    grid on;
    hold off;
    sInputChar = char(sInput);
    formatFileName = "time_series_castner_%s.jpg";
    fileName = sprintf(formatFileName,sInputChar);
    fullFilePath = fullfile(folderName, fileName);
    saveas(gcf, fullFilePath, 'jpg');
end
%% Time-series figure: Gulkana Watershed - done
columnNames = vColumnLabels;
folderName = 'U:/GoA plots/NewPlots';
for dPosition = 1:numel(columnNames)
    sInput = columnNames{dPosition};
    sTitle = ['Element ', sInput];
    figure;
    y_G1_2022 = mG1_2022(:, dPosition);
    y_G1_2022 = fillmissing(y_G1_2022, 'linear');
    %plot(mG1_Date22, y_G1_2022, 'o-','DisplayName', 'Gulkana 1','Color','k','MarkerFaceColor',[0.8,0.4,0]);
    plot(mG1_Date22, y_G1_2022, 'o-','DisplayName', 'Gulkana 1','Color',[0.8,0.4,0],'MarkerFaceColor',[0.8,0.4,0],'MarkerEdgeColor','k');
    hold on;
    y_G1_2023 = mG1_2023(:, dPosition);
    y_G1_2023 = fillmissing(y_G1_2023, 'linear');
    %plot(mG1_Date23, y_G1_2023, 'o-','Color','k','MarkerFaceColor',[0.8,0.4,0],'HandleVisibility', 'off');
    plot(mG1_Date23, y_G1_2023, 'o-','Color',[0.8,0.4,0],'MarkerFaceColor',[0.8,0.4,0],'HandleVisibility', 'off','MarkerEdgeColor','k');
    y_G2_2022 = mG2_2022(:, dPosition);
    y_G2_2022 = fillmissing(y_G2_2022, 'linear');
    %plot(mG2_Date22, y_G2_2022, 'o-', 'DisplayName', 'Gulkana 2','Color','k','MarkerFaceColor',[0,0.6,0.5]);
    plot(mG2_Date22, y_G2_2022, 'o-', 'DisplayName', 'Gulkana 2','Color',[0,0.6,0.5],'MarkerFaceColor',[0,0.6,0.5],'MarkerEdgeColor','k');
    y_G2_2023 = mG2_2023(:, dPosition);
    y_G2_2023 = fillmissing(y_G2_2023, 'linear');
    %plot(mG2_Date23, y_G2_2023, 'o-','Color','k','MarkerFaceColor',[0,0.6,0.5],'HandleVisibility', 'off');
    plot(mG2_Date23, y_G2_2023, 'o-','Color',[0,0.6,0.5],'MarkerFaceColor',[0,0.6,0.5],'HandleVisibility', 'off','MarkerEdgeColor','k');
    y_G3_2022 = mG3_2022(:, dPosition);
    y_G3_2022 = fillmissing(y_G3_2022, 'linear');
    %plot(mG3_Date22, y_G3_2022, 'o-', 'DisplayName', 'Gulkana 3','Color','k','MarkerFaceColor',[0.95,0.90,0.25]);
    plot(mG3_Date22, y_G3_2022, 'o-', 'DisplayName', 'Gulkana 3','Color',[0.95,0.90,0.25],'MarkerFaceColor',[0.95,0.90,0.25],'MarkerEdgeColor','k');
    y_G3_2023 = mG3_2023(:, dPosition);
    y_G3_2023 = fillmissing(y_G3_2023, 'linear');
    %plot(mG3_Date23, y_G3_2023, 'o-','Color','k','MarkerFaceColor',[0.95,0.90,0.25],'HandleVisibility', 'off');
    plot(mG3_Date23, y_G3_2023, 'o-','Color',[0.95,0.90,0.25],'MarkerFaceColor',[0.95,0.90,0.25],'HandleVisibility', 'off','MarkerEdgeColor','k');
    xlabel('Time');  % Label for x-axis
    ylabel('Concentration in mg/L');  % Label for y-axis
    title(sInput);
    legend('Location', 'none');
    legend('Position', [0.5 0.5 0.1 0.1], 'Location', 'none');
    grid on;
    hold off;
    sInputChar = char(sInput);
    formatFileName = "time_series_gulk_%s.jpg";
    fileName = sprintf(formatFileName,sInputChar);
    fullFilePath = fullfile(folderName, fileName);
    saveas(gcf, fullFilePath, 'jpg');
end
%% Spatial figures: Matanuska Watershed - done
% Extract month from dates
columnNames = vColumnLabels;
folderName = 'U:/GoA plots/NewPlots';
years = year(mMatDates); % Extract years from dates

% Separate data for 2022 and 2023
data2022 = mMatData(years == 2022, :);
data2023 = mMatData(years == 2023, :);

% Define the months for which to create lines
targetMonths = [5, 6, 7, 8, 9];
colors = {[0.8, 0.6, 0.7],[0.9, 0.6, 0],[0.95, 0.9, 0.25],[0, 0.6, 0.5],[0.35, 0.7, 0.9]};
monthLabels = {'May', 'June', 'July', 'August', 'September'};

% Separate data for each year
uniqueYears = unique(years);
for dPosition = 1:numel(columnNames)
    figure;
    yLimits = []; % Initialize an array to store y-axis limits of all subplots
    for yearIndex = 1:numel(uniqueYears) % Loop for each year
        currentYear = uniqueYears(yearIndex);
        currentData = mMatData(years == currentYear, :); % Use currentData instead of data

        % Extract month and year from dates
        months = month(mMatDates(years == currentYear));
        uniqueSites = unique(mMatSites); % Get unique site names and sort them alphabetically
        sortedSites = sort(uniqueSites);
        subplot(1, numel(uniqueYears), yearIndex); % Create subplot for current year
        for i = 1:numel(targetMonths)
            xCoords = [];
            yCoords = [];
            for j = 1:numel(sortedSites)
                currentSite = sortedSites{j};
                monthData = currentData(strcmp(mMatSites(years == currentYear), currentSite) & months == targetMonths(i), dPosition);
                nonNaNIndices = ~isnan(monthData); % Store x and y coordinates, skipping NaN values
                xCoords = [xCoords, repmat(j, 1, sum(nonNaNIndices))];
                yCoords = [yCoords, monthData(nonNaNIndices)];
            end
            if ~isempty(xCoords) % Plot the data
                plot(xCoords, yCoords, 'o-', 'Color', colors{i}, 'MarkerFaceColor', colors{i}, 'MarkerEdgeColor', 'k', 'DisplayName', monthLabels{i});
                hold on;
            end
        end
        xticks(1:numel(sortedSites));
        xticklabels(sortedSites);
        title([num2str(currentYear)]);

        % Get the y-axis limits of the current subplot and store them
        yLimits = [yLimits, ylim];
        
        % Add y-axis label and tick labels only to the first subplot
        if yearIndex == 1
            ylabel('Concentration (mg/L)');
        end
    end
    
    % Set the same y-axis limits for all subplots
    y_min = min(yLimits);
    y_max = max(yLimits);
    for yearIndex = 1:numel(uniqueYears)
        subplot(1, numel(uniqueYears), yearIndex);
        ylim([y_min, y_max]);
        % Add right y-axis labels for the second subplot
        if yearIndex == 2
            ax = gca;
            ax.YAxisLocation = 'right';
            set(gca, 'YTickLabel', []); % Remove numeric tick labels on the right y-axis
            set(get(gca, 'YLabel'), 'String', ''); % Remove the y-axis label
        end
    end
    
    hold off;
    sInput = columnNames{dPosition};
    sInputChar = char(sInput);
    formatFileName = "spatial_mat_%s.svg";
    fileName = sprintf(formatFileName, sInputChar);
    fullFilePath = fullfile(folderName, fileName);
    saveas(gcf, fullFilePath, 'svg');
end
%% Matanuska subplots
% Extract month from dates
columnNames = vColumnLabels;
folderName = 'U:/GoA plots/NewPlots';
years = year(mMatDates); % Extract years from dates

% Define the months for which to create lines
targetMonths = [5, 6, 7, 8, 9];
colors = {[0.8, 0.6, 0.7], [0.9, 0.6, 0], [0.95, 0.9, 0.25], [0, 0.6, 0.5], [0.35, 0.7, 0.9]};
monthLabels = {'May', 'June', 'July', 'August', 'September'};

% Define the dPosition values to iterate over
dPositions = [42, 29, 11, 60]; %for paper[42, 29, 11, 60], cations [26, 25, 28, 29], anions [11, 14, 18, 13, 16]

% Create a figure with 8 subplots (4 rows x 2 columns)
figure;
% Set the figure size: [left, bottom, width, height]
set(gcf, 'Position', [100, 100, 600, 800]); % Example dimensions: width, height

for dpIndex = 1:numel(dPositions)
    dPosition = dPositions(dpIndex);
    sInput = columnNames{dPosition}; % Get the label for the current dPosition
    yLimits = []; % Initialize an array to store y-axis limits of all subplots

    for yearIndex = 1:2 % 1 for 2022, 2 for 2023
        if yearIndex == 1
            currentYear = 2022;
        else
            currentYear = 2023;
        end
        
        % Filter data for the current year and dPosition
        currentData = mMatData(years == currentYear, :);
        currentSites = mMatSites(years == currentYear);
        months = month(mMatDates(years == currentYear));

        % Get unique sites and sort them
        uniqueSites = unique(currentSites); % Get unique site names and sort them alphabetically
        sortedSites = sort(uniqueSites);

        % Calculate the subplot index for 4 rows x 2 columns layout
        subplotIndex = (dpIndex - 1) * 2 + yearIndex;
        subplot(4, 2, subplotIndex); % Create subplot for current year and dPosition

        % Plot data for the current dPosition and year
        for i = 1:numel(targetMonths)
            xCoords = [];
            yCoords = [];
            for j = 1:numel(sortedSites)
                currentSite = sortedSites{j};
                % Filter data for the current site and month
                monthData = currentData(strcmp(currentSites, currentSite) & months == targetMonths(i), dPosition);
                nonNaNIndices = ~isnan(monthData); % Store x and y coordinates, skipping NaN values
                xCoords = [xCoords, repmat(j, 1, sum(nonNaNIndices))];
                yCoords = [yCoords, monthData(nonNaNIndices)];
            end
            if ~isempty(xCoords) % Plot the data
                plot(xCoords, yCoords, 'o-', 'Color', colors{i}, 'MarkerFaceColor', colors{i}, 'MarkerEdgeColor', 'k', 'DisplayName', monthLabels{i});
                hold on;
            end
        end
        
        xticks(1:numel(sortedSites));
        xticklabels(sortedSites);
        % Update title to use sInput and year
        if contains(sInput, 'HCO3')
            sInput = strrep(sInput, 'HCO3', 'HCO_3^-');
        end

        if contains(sInput, 'Ca') && ~contains(sInput, 'Ca^{+2}')
            sInput = strrep(sInput, 'Ca', 'Ca^{+2}');
        end

        title([sInput ' ' num2str(currentYear)], 'Interpreter', 'tex');


        % Get the y-axis limits of the current subplot and store them
        yLimits = [yLimits, ylim];
        
        % Add y-axis label and tick labels only to the first column subplots
        if yearIndex == 1
            ylabel('Concentration (mg/L)');
        end
        
        % Add the legend to the second subplot (2023 plot for dPosition = 42)
        if dpIndex == 1 && yearIndex == 1
            legendHandle = legend('Location', 'northeast');
            set(legendHandle, 'FontSize',5.5);
            set(legendHandle, 'ItemTokenSize', [5,5]);
        end
    end

    % Set the same y-axis limits for all subplots in this pair
    y_min = min(yLimits);
    y_max = max(yLimits);
    for yearIndex = 1:2
        subplotIndex = (dpIndex - 1) * 2 + yearIndex;
        subplot(4, 2, subplotIndex);
        ylim([y_min, y_max]);
        % Add right y-axis labels for the second column of each row
        if yearIndex == 2
            ax = gca;
            ax.YAxisLocation = 'right';
            set(gca, 'YTickLabel', []); % Remove numeric tick labels on the right y-axis
            set(get(gca, 'YLabel'), 'String', ''); % Remove the y-axis label
        end
    end
end

% Save the figure
saveas(gcf, fullfile(folderName, 'spatial_mat_combined.svg'), 'svg');

%% Matanuska subplots - heavy metals
% Extract month from dates
columnNames = vColumnLabels;
folderName = 'U:/GoA plots/NewPlots';
years = year(mMatDates); % Extract years from dates

% Define the months for which to create lines
targetMonths = [5, 6, 7, 8, 9];
colors = {[0.8, 0.6, 0.7], [0.9, 0.6, 0], [0.95, 0.9, 0.25], [0, 0.6, 0.5], [0.35, 0.7, 0.9]};
monthLabels = {'May', 'June', 'July', 'August', 'September'};

% Define the dPosition values to iterate over
dPositions = [32, 37, 33, 36, 38, 63, 45, 42, 35]; %for paper[42, 29, 11, 60], cations [26, 25, 28, 29], anions [11, 14, 18, 13, 16]

% Create a figure with 8 subplots (5 rows x 4 columns)
figure;
% Set the figure size: [left, bottom, width, height]
set(gcf, 'Position', [100, 100, 1200, 5000]); % Example dimensions: width, height

for dpIndex = 1:numel(dPositions)
    dPosition = dPositions(dpIndex);
    sInput = columnNames{dPosition}; % Get the label for the current dPosition
    yLimits = []; % Initialize an array to store y-axis limits of all subplots

    for yearIndex = 1:2 % 1 for 2022, 2 for 2023
        if yearIndex == 1
            currentYear = 2022;
        else
            currentYear = 2023;
        end
        
        % Filter data for the current year and dPosition
        currentData = mMatData(years == currentYear, :);
        currentSites = mMatSites(years == currentYear);
        months = month(mMatDates(years == currentYear));

        % Get unique sites and sort them
        uniqueSites = unique(currentSites); % Get unique site names and sort them alphabetically
        sortedSites = sort(uniqueSites);

        % Calculate the subplot index for 4 rows x 2 columns layout
        subplotIndex = (dpIndex - 1) * 2 + yearIndex;
        if subplotIndex > 18
            break;
        end
        subplot(5, 4, subplotIndex); % Create subplot for current year and dPosition

        % Plot data for the current dPosition and year
        for i = 1:numel(targetMonths)
            xCoords = [];
            yCoords = [];
            for j = 1:numel(sortedSites)
                currentSite = sortedSites{j};
                % Filter data for the current site and month
                monthData = currentData(strcmp(currentSites, currentSite) & months == targetMonths(i), dPosition);
                nonNaNIndices = ~isnan(monthData); % Store x and y coordinates, skipping NaN values
                xCoords = [xCoords, repmat(j, 1, sum(nonNaNIndices))];
                yCoords = [yCoords, monthData(nonNaNIndices)];
            end
            if ~isempty(xCoords) % Plot the data
                plot(xCoords, yCoords, 'o-', 'Color', colors{i}, 'MarkerFaceColor', colors{i}, 'MarkerEdgeColor', 'k', 'DisplayName', monthLabels{i});
                hold on;
            end
        end
        
        xticks(1:numel(sortedSites));
        xticklabels(sortedSites);
        % Update title to use sInput and year
        if contains(sInput, 'HCO3')
            sInput = strrep(sInput, 'HCO3', 'HCO_3');
        end
        title([sInput ' ' num2str(currentYear)], 'Interpreter', 'tex');

        % Get the y-axis limits of the current subplot and store them
        yLimits = [yLimits, ylim];
        
        % Add y-axis label and tick labels only to the first column subplots
        if yearIndex == 1
            ylabel('Concentration (mg/L)');
        end
        
        % Add the legend to the second subplot (2023 plot for dPosition = 42)
        if dpIndex == 1 && yearIndex == 1
            legendHandle = legend('Location', 'northwest');
            set(legendHandle, 'FontSize',5.5);
            set(legendHandle, 'ItemTokenSize', [5,5]);
        end
    end

    % Set the same y-axis limits for all subplots in this pair
    y_min = min(yLimits);
    y_max = max(yLimits);
    for yearIndex = 1:2
        subplotIndex = (dpIndex - 1) * 2 + yearIndex;
        subplot(5, 4, subplotIndex);
        ylim([y_min, y_max]);
        % Add right y-axis labels for the second column of each row
        if yearIndex == 2
            ax = gca;
            ax.YAxisLocation = 'right';
            set(gca, 'YTickLabel', []); % Remove numeric tick labels on the right y-axis
            set(get(gca, 'YLabel'), 'String', ''); % Remove the y-axis label
        end
    end
end

% Save the figure
saveas(gcf, fullfile(folderName, 'spatial_mat_heavymetals.jpg'), 'jpg');

%% Matanuska subplots - REE
% Extract month from dates
columnNames = vColumnLabels;
folderName = 'U:/GoA plots/NewPlots';
years = year(mMatDates); % Extract years from dates

% Define the months for which to create lines
targetMonths = [5, 6, 7, 8, 9];
colors = {[0.8, 0.6, 0.7], [0.9, 0.6, 0], [0.95, 0.9, 0.25], [0, 0.6, 0.5], [0.35, 0.7, 0.9]};
monthLabels = {'May', 'June', 'July', 'August', 'September'};

% Define the dPosition values to iterate over
dPositions = [30,50,57,60,61,49,56,51,55,59,41,52,53,54,58]; %15 elements, 30 plots

% Create a figure with 8 subplots (5 rows x 6 columns)
figure;
% Set the figure size: [left, bottom, width, height]
set(gcf, 'Position', [100, 100, 1200, 5000]); % Example dimensions: width, height

for dpIndex = 1:numel(dPositions)
    dPosition = dPositions(dpIndex);
    sInput = columnNames{dPosition}; % Get the label for the current dPosition
    yLimits = []; % Initialize an array to store y-axis limits of all subplots

    for yearIndex = 1:2 % 1 for 2022, 2 for 2023
        if yearIndex == 1
            currentYear = 2022;
        else
            currentYear = 2023;
        end
        
        % Filter data for the current year and dPosition
        currentData = mMatData(years == currentYear, :);
        currentSites = mMatSites(years == currentYear);
        months = month(mMatDates(years == currentYear));

        % Get unique sites and sort them
        uniqueSites = unique(currentSites); % Get unique site names and sort them alphabetically
        sortedSites = sort(uniqueSites);

        % Calculate the subplot index for 4 rows x 2 columns layout
        subplotIndex = (dpIndex - 1) * 2 + yearIndex;
        subplot(5, 6, subplotIndex); % Create subplot for current year and dPosition

        % Plot data for the current dPosition and year
        for i = 1:numel(targetMonths)
            xCoords = [];
            yCoords = [];
            for j = 1:numel(sortedSites)
                currentSite = sortedSites{j};
                % Filter data for the current site and month
                monthData = currentData(strcmp(currentSites, currentSite) & months == targetMonths(i), dPosition);
                nonNaNIndices = ~isnan(monthData); % Store x and y coordinates, skipping NaN values
                xCoords = [xCoords, repmat(j, 1, sum(nonNaNIndices))];
                yCoords = [yCoords, monthData(nonNaNIndices)];
            end
            if ~isempty(xCoords) % Plot the data
                plot(xCoords, yCoords, 'o-', 'Color', colors{i}, 'MarkerFaceColor', colors{i}, 'MarkerEdgeColor', 'k', 'DisplayName', monthLabels{i});
                hold on;
            end
        end
        
        xticks(1:numel(sortedSites));
        xticklabels(sortedSites);
        % Update title to use sInput and year
        if contains(sInput, 'HCO3')
            sInput = strrep(sInput, 'HCO3', 'HCO_3');
        end
        title([sInput ' ' num2str(currentYear)], 'Interpreter', 'tex');

        % Get the y-axis limits of the current subplot and store them
        yLimits = [yLimits, ylim];
        
        % Add y-axis label and tick labels only to the first column subplots
        if yearIndex == 1
            ylabel('Concentration (mg/L)');
        end
        
        % Add the legend to the second subplot (2023 plot for dPosition = 42)
        if dpIndex == 2 && yearIndex == 1
            legendHandle = legend('Location', 'northeast');
            set(legendHandle, 'FontSize',5.5);
            set(legendHandle, 'ItemTokenSize', [5,5]);
        end
    end

    % Set the same y-axis limits for all subplots in this pair
    y_min = min(yLimits);
    y_max = max(yLimits);
    for yearIndex = 1:2
        subplotIndex = (dpIndex - 1) * 2 + yearIndex;
        subplot(5, 6, subplotIndex);
        ylim([y_min, y_max]);
        % Add right y-axis labels for the second column of each row
        if yearIndex == 2
            ax = gca;
            ax.YAxisLocation = 'right';
            set(gca, 'YTickLabel', []); % Remove numeric tick labels on the right y-axis
            set(get(gca, 'YLabel'), 'String', ''); % Remove the y-axis label
        end
    end
end

% Save the figure
saveas(gcf, fullfile(folderName, 'spatial_mat_REE.svg'), 'svg');

%% Knik
columnNames = vColumnLabels;
folderName = 'U:/GoA plots/NewPlots';
years = year(mKnikDates); % Extract years from dates
data2022 = mKnikData(years == 2022, :);
data2023 = mKnikData(years == 2023, :);
targetMonths = [5, 6, 7, 8, 9];
colors = {[0.8, 0.6, 0.7],[0.9, 0.6, 0],[0.95, 0.9, 0.25],[0, 0.6, 0.5],[0.35, 0.7, 0.9]};
monthLabels = {'May', 'June', 'July', 'August', 'September'};
uniqueYears = unique(years);
for dPosition = 1:numel(columnNames)
    figure;
    yLimits = []; % Initialize an array to store y-axis limits of all subplots
    for yearIndex = 1:numel(uniqueYears) % Loop for each year
        currentYear = uniqueYears(yearIndex);
        currentData = mKnikData(years == currentYear, :);
        months = month(mKnikDates(years == currentYear));
        uniqueSites = unique(mKnikSites); % Get unique site names and sort them alphabetically
        sortedSites = sort(uniqueSites);
        subplot(1, numel(uniqueYears), yearIndex); % Create subplot for current year
        for i = 1:numel(targetMonths)
            xCoords = [];
            yCoords = [];
            for j = 1:numel(sortedSites)
                currentSite = sortedSites{j};
                monthData = currentData(strcmp(mKnikSites(years == currentYear), currentSite) & months == targetMonths(i), dPosition);
                nonNaNIndices = ~isnan(monthData); % Store x and y coordinates, skipping NaN values
                xCoords = [xCoords, repmat(j, 1, sum(nonNaNIndices))];
                yCoords = [yCoords, monthData(nonNaNIndices)];
            end
            if ~isempty(xCoords) % Plot the data
                plot(xCoords, yCoords, 'o-', 'Color', colors{i}, 'MarkerFaceColor', colors{i}, 'MarkerEdgeColor', 'k', 'DisplayName', monthLabels{i});
                hold on;
            end
        end
        xticks(1:numel(sortedSites));
        xticklabels(sortedSites);
        if contains(sInput, 'HCO3')
            sInput = strrep(sInput, 'HCO3', 'HCO_3');
        end
        title([sInput ' - ' num2str(currentYear)], 'Interpreter', 'tex');;

        % Get the y-axis limits of the current subplot and store them
        yLimits = [yLimits, ylim];
        
        % Add y-axis label and tick labels only to the first subplot
        if yearIndex == 1
            ylabel('Concentration (mg/L)');
        end
    end
    y_min = min(yLimits);
    y_max = max(yLimits);
    for yearIndex = 1:numel(uniqueYears)
        subplot(1, numel(uniqueYears), yearIndex);
        ylim([y_min, y_max]);
        % Add right y-axis labels for the second subplot
        if yearIndex == 2
            ax = gca;
            ax.YAxisLocation = 'right';
            set(gca, 'YTickLabel', []); % Remove numeric tick labels on the right y-axis
            set(get(gca, 'YLabel'), 'String', ''); % Remove the y-axis label
        end
    end
    hold off;
    sInput = columnNames{dPosition};
    sInputChar = char(sInput);
    formatFileName = "spatial_knik_%s.svg";
    fileName = sprintf(formatFileName, sInputChar);
    fullFilePath = fullfile(folderName, fileName);
    saveas(gcf, fullFilePath, 'svg');
end
%% Knik subplots
% Extract month from dates
columnNames = vColumnLabels;
folderName = 'U:/GoA plots/NewPlots';
years = year(mKnikDates); % Extract years from dates

% Define the months for which to create lines
targetMonths = [5, 6, 7, 8, 9];
colors = {[0.8, 0.6, 0.7], [0.9, 0.6, 0], [0.95, 0.9, 0.25], [0, 0.6, 0.5], [0.35, 0.7, 0.9]};
monthLabels = {'May', 'June', 'July', 'August', 'September'};

% Define the dPosition values to iterate over
dPositions = [42, 29, 11, 60];

% Create a figure with 8 subplots (4 rows x 2 columns)
figure;
% Set the figure size: [left, bottom, width, height]
set(gcf, 'Position', [100, 100, 600, 800]); % Width = 600, Height = 800

for dpIndex = 1:numel(dPositions)
    dPosition = dPositions(dpIndex);
    sInput = columnNames{dPosition}; % Get the label for the current dPosition
    yLimits = []; % Initialize an array to store y-axis limits of all subplots

    for yearIndex = 1:2 % 1 for 2022, 2 for 2023
        if yearIndex == 1
            currentYear = 2022;
        else
            currentYear = 2023;
        end
        
        % Filter data for the current year and dPosition
        currentData = mKnikData(years == currentYear, :);
        currentSites = mKnikSites(years == currentYear);
        months = month(mKnikDates(years == currentYear));

        % Get unique sites and sort them
        uniqueSites = unique(currentSites); % Get unique site names and sort them alphabetically
        sortedSites = sort(uniqueSites);

        % Calculate the subplot index for 4 rows x 2 columns layout
        subplotIndex = (dpIndex - 1) * 2 + yearIndex;
        subplot(4, 2, subplotIndex); % Create subplot for current year and dPosition

        % Plot data for the current dPosition and year
        hold on;
        for i = 1:numel(targetMonths)
            xCoords = [];
            yCoords = [];
            for j = 1:numel(sortedSites)
                currentSite = sortedSites{j};
                % Filter data for the current site and month
                monthData = currentData(strcmp(currentSites, currentSite) & months == targetMonths(i), dPosition);
                % Check if monthData is empty or contains only NaN values
                if ~isempty(monthData) && any(~isnan(monthData))
                    nonNaNIndices = ~isnan(monthData); % Store x and y coordinates, skipping NaN values
                    xCoords = [xCoords, repmat(j, 1, sum(nonNaNIndices))];
                    yCoords = [yCoords, monthData(nonNaNIndices)];
                end
            end
            if ~isempty(xCoords) % Plot the data
                plot(xCoords, yCoords, 'o-', 'Color', colors{i}, 'MarkerFaceColor', colors{i}, 'MarkerEdgeColor', 'k', 'DisplayName', monthLabels{i});
            end
        end
        
        xticks(1:numel(sortedSites));
        xticklabels(sortedSites);        
        
        if contains(sInput, 'HCO3')
            sInput = strrep(sInput, 'HCO3', 'HCO_3^-');
        end

        if contains(sInput, 'Ca') && ~contains(sInput, 'Ca^{+2}')
            sInput = strrep(sInput, 'Ca', 'Ca^{+2}');
        end

        title([sInput ' ' num2str(currentYear)], 'Interpreter', 'tex');

        % Get the y-axis limits of the current subplot and store them
        yLimits = [yLimits, ylim];
        
        % Add y-axis label and tick labels only to the first column subplots
        if yearIndex == 1
            ylabel('Concentration (mg/L)');
        end
        if dpIndex == 1 && yearIndex == 1
            legendHandle = legend('Location', 'west');
            set(legendHandle, 'FontSize',5.5);
            set(legendHandle, 'ItemTokenSize', [3,3]);
        end
    end

    % Set the same y-axis limits for all subplots in this pair
    y_min = min(yLimits);
    y_max = max(yLimits);
    for yearIndex = 1:2
        subplotIndex = (dpIndex - 1) * 2 + yearIndex;
        subplot(4, 2, subplotIndex);
        ylim([y_min, y_max]);
        % Add right y-axis labels for the second column of each row
        if yearIndex == 2
            ax = gca;
            ax.YAxisLocation = 'right';
            set(gca, 'YTickLabel', []); % Remove numeric tick labels on the right y-axis
            set(get(gca, 'YLabel'), 'String', ''); % Remove the y-axis label
        end
        % Add black outline around each subplot with thinner line width
        box on;
        set(gca, 'LineWidth', 0.75, 'Box', 'on'); % Decrease line width to 0.75
    end

    hold off;
end

% Save the figure
formatFileName = "spatial_knik_subplot.svg";
fileName = sprintf(formatFileName);
fullFilePath = fullfile(folderName, fileName);
saveas(gcf, fullFilePath, 'svg');

%% Knik subplots, heavy metals
% Extract month from dates
columnNames = vColumnLabels;
folderName = 'U:/GoA plots/NewPlots';
years = year(mKnikDates); % Extract years from dates

% Define the months for which to create lines
targetMonths = [5, 6, 7, 8, 9];
colors = {[0.8, 0.6, 0.7], [0.9, 0.6, 0], [0.95, 0.9, 0.25], [0, 0.6, 0.5], [0.35, 0.7, 0.9]};
monthLabels = {'May', 'June', 'July', 'August', 'September'};

% Define the dPosition values to iterate over
dPositions = [32, 37, 33, 36, 38, 63, 45, 42, 35];

% Create a figure with 8 subplots (4 rows x 2 columns)
figure;
% Set the figure size: [left, bottom, width, height]
set(gcf, 'Position', [100, 100, 1200, 6000]); % Width = 600, Height = 800

for dpIndex = 1:numel(dPositions)
    dPosition = dPositions(dpIndex);
    sInput = columnNames{dPosition}; % Get the label for the current dPosition
    yLimits = []; % Initialize an array to store y-axis limits of all subplots

    for yearIndex = 1:2 % 1 for 2022, 2 for 2023
        if yearIndex == 1
            currentYear = 2022;
        else
            currentYear = 2023;
        end
        
        % Filter data for the current year and dPosition
        currentData = mKnikData(years == currentYear, :);
        currentSites = mKnikSites(years == currentYear);
        months = month(mKnikDates(years == currentYear));

        % Get unique sites and sort them
        uniqueSites = unique(currentSites); % Get unique site names and sort them alphabetically
        sortedSites = sort(uniqueSites);

        % Calculate the subplot index for 4 rows x 2 columns layout
        subplotIndex = (dpIndex - 1) * 2 + yearIndex;
        if subplotIndex > 18
            break;
        end
        subplot(5, 4, subplotIndex); % Create subplot for current year and dPosition

        % Plot data for the current dPosition and year
        hold on;
        for i = 1:numel(targetMonths)
            xCoords = [];
            yCoords = [];
            for j = 1:numel(sortedSites)
                currentSite = sortedSites{j};
                % Filter data for the current site and month
                monthData = currentData(strcmp(currentSites, currentSite) & months == targetMonths(i), dPosition);
                % Check if monthData is empty or contains only NaN values
                if ~isempty(monthData) && any(~isnan(monthData))
                    nonNaNIndices = ~isnan(monthData); % Store x and y coordinates, skipping NaN values
                    xCoords = [xCoords, repmat(j, 1, sum(nonNaNIndices))];
                    yCoords = [yCoords, monthData(nonNaNIndices)];
                end
            end
            if ~isempty(xCoords) % Plot the data
                plot(xCoords, yCoords, 'o-', 'Color', colors{i}, 'MarkerFaceColor', colors{i}, 'MarkerEdgeColor', 'k', 'DisplayName', monthLabels{i});
            end
        end
        
        xticks(1:numel(sortedSites));
        xticklabels(sortedSites);        
        
        if contains(sInput, 'HCO3')
            sInput = strrep(sInput, 'HCO3', 'HCO_3');
        end
        title([sInput ' ' num2str(currentYear)], 'Interpreter', 'tex');
        
        % Get the y-axis limits of the current subplot and store them
        yLimits = [yLimits, ylim];
        
        % Add y-axis label and tick labels only to the first column subplots
        if yearIndex == 1
            ylabel('Concentration (mg/L)');
        end
        if dpIndex == 1 && yearIndex == 1
            legendHandle = legend('Location', 'west');
            set(legendHandle, 'FontSize',5.5);
            set(legendHandle, 'ItemTokenSize', [3,3]);
        end
    end

    % Set the same y-axis limits for all subplots in this pair
    y_min = min(yLimits);
    y_max = max(yLimits);
    for yearIndex = 1:2
        subplotIndex = (dpIndex - 1) * 2 + yearIndex;
        subplot(5, 4, subplotIndex);
        ylim([y_min, y_max]);
        % Add right y-axis labels for the second column of each row
        if yearIndex == 2
            ax = gca;
            ax.YAxisLocation = 'right';
            set(gca, 'YTickLabel', []); % Remove numeric tick labels on the right y-axis
            set(get(gca, 'YLabel'), 'String', ''); % Remove the y-axis label
        end
        % Add black outline around each subplot with thinner line width
        box on;
        set(gca, 'LineWidth', 0.75, 'Box', 'on'); % Decrease line width to 0.75
    end

    hold off;
end

% Save the figure
formatFileName = "spatial_knik_heavymetals.jpg";
fileName = sprintf(formatFileName);
fullFilePath = fullfile(folderName, fileName);
saveas(gcf, fullFilePath, 'jpg');

%% Knik subplots, REE
% Extract month from dates
columnNames = vColumnLabels;
folderName = 'U:/GoA plots/NewPlots';
years = year(mKnikDates); % Extract years from dates

% Define the months for which to create lines
targetMonths = [5, 6, 7, 8, 9];
colors = {[0.8, 0.6, 0.7], [0.9, 0.6, 0], [0.95, 0.9, 0.25], [0, 0.6, 0.5], [0.35, 0.7, 0.9]};
monthLabels = {'May', 'June', 'July', 'August', 'September'};

% Define the dPosition values to iterate over
dPositions = [30,50,57,60,61,49,56,51,55,59,41,52,53,54,58];

% Create a figure with 8 subplots (4 rows x 2 columns)
figure;
% Set the figure size: [left, bottom, width, height]
set(gcf, 'Position', [100, 100, 1200, 6000]); % Width = 600, Height = 800

for dpIndex = 1:numel(dPositions)
    dPosition = dPositions(dpIndex);
    sInput = columnNames{dPosition}; % Get the label for the current dPosition
    yLimits = []; % Initialize an array to store y-axis limits of all subplots

    for yearIndex = 1:2 % 1 for 2022, 2 for 2023
        if yearIndex == 1
            currentYear = 2022;
        else
            currentYear = 2023;
        end
        
        % Filter data for the current year and dPosition
        currentData = mKnikData(years == currentYear, :);
        currentSites = mKnikSites(years == currentYear);
        months = month(mKnikDates(years == currentYear));

        % Get unique sites and sort them
        uniqueSites = unique(currentSites); % Get unique site names and sort them alphabetically
        sortedSites = sort(uniqueSites);

        % Calculate the subplot index for 4 rows x 2 columns layout
        subplotIndex = (dpIndex - 1) * 2 + yearIndex;
        subplot(5, 6, subplotIndex); % Create subplot for current year and dPosition

        % Plot data for the current dPosition and year
        hold on;
        for i = 1:numel(targetMonths)
            xCoords = [];
            yCoords = [];
            for j = 1:numel(sortedSites)
                currentSite = sortedSites{j};
                % Filter data for the current site and month
                monthData = currentData(strcmp(currentSites, currentSite) & months == targetMonths(i), dPosition);
                % Check if monthData is empty or contains only NaN values
                if ~isempty(monthData) && any(~isnan(monthData))
                    nonNaNIndices = ~isnan(monthData); % Store x and y coordinates, skipping NaN values
                    xCoords = [xCoords, repmat(j, 1, sum(nonNaNIndices))];
                    yCoords = [yCoords, monthData(nonNaNIndices)];
                end
            end
            if ~isempty(xCoords) % Plot the data
                plot(xCoords, yCoords, 'o-', 'Color', colors{i}, 'MarkerFaceColor', colors{i}, 'MarkerEdgeColor', 'k', 'DisplayName', monthLabels{i});
            end
        end
        
        xticks(1:numel(sortedSites));
        xticklabels(sortedSites);        
        
        if contains(sInput, 'HCO3')
            sInput = strrep(sInput, 'HCO3', 'HCO_3');
        end
        title([sInput ' ' num2str(currentYear)], 'Interpreter', 'tex');
        
        % Get the y-axis limits of the current subplot and store them
        yLimits = [yLimits, ylim];
        
        % Add y-axis label and tick labels only to the first column subplots
        if yearIndex == 1
            ylabel('Concentration (mg/L)');
        end
        if dpIndex == 5 && yearIndex == 2
            legendHandle = legend('Location', 'northwest');
            set(legendHandle, 'FontSize',5.5);
            set(legendHandle, 'ItemTokenSize', [3,3]);
        end
    end

    % Set the same y-axis limits for all subplots in this pair
    y_min = min(yLimits);
    y_max = max(yLimits);
    for yearIndex = 1:2
        subplotIndex = (dpIndex - 1) * 2 + yearIndex;
        subplot(5, 6, subplotIndex);
        ylim([y_min, y_max]);
        % Add right y-axis labels for the second column of each row
        if yearIndex == 2
            ax = gca;
            ax.YAxisLocation = 'right';
            set(gca, 'YTickLabel', []); % Remove numeric tick labels on the right y-axis
            set(get(gca, 'YLabel'), 'String', ''); % Remove the y-axis label
        end
        % Add black outline around each subplot with thinner line width
        box on;
        set(gca, 'LineWidth', 0.75, 'Box', 'on'); % Decrease line width to 0.75
    end

    hold off;
end

% Save the figure
formatFileName = "spatial_knik_REE.svg";
fileName = sprintf(formatFileName);
fullFilePath = fullfile(folderName, fileName);
saveas(gcf, fullFilePath, 'svg');

%% Spatial figures: Little Susitna Watershed
columnNames = vColumnLabels;
folderName = 'U:/GoA plots/NewPlots';
years = year(mLSDates); % Extract years from dates
data2022 = mLSData(years == 2022, :);
data2023 = mLSData(years == 2023, :);
targetMonths = [5, 6, 7, 8, 9];
colors = {[0.8, 0.6, 0.7],[0.9, 0.6, 0],[0.95, 0.9, 0.25],[0, 0.6, 0.5],[0.35, 0.7, 0.9]};
monthLabels = {'May', 'June', 'July', 'August', 'September'};
uniqueYears = unique(years);
for dPosition = 1:numel(columnNames)
    figure;
    yLimits = []; % Initialize an array to store y-axis limits of all subplots
    for yearIndex = 1:numel(uniqueYears) % Loop for each year
        currentYear = uniqueYears(yearIndex);
        currentData = mLSData(years == currentYear, :);
        months = month(mLSDates(years == currentYear));
        uniqueSites = unique(mLSSites); % Get unique site names and sort them alphabetically
        sortedSites = sort(uniqueSites);
        subplot(1, numel(uniqueYears), yearIndex); % Create subplot for current year
        for i = 1:numel(targetMonths)
            xCoords = [];
            yCoords = [];
            for j = 1:numel(sortedSites)
                currentSite = sortedSites{j};
                monthData = currentData(strcmp(mLSSites(years == currentYear), currentSite) & months == targetMonths(i), dPosition);
                nonNaNIndices = ~isnan(monthData); % Store x and y coordinates, skipping NaN values
                xCoords = [xCoords, repmat(j, 1, sum(nonNaNIndices))];
                yCoords = [yCoords, monthData(nonNaNIndices)];
            end
            if ~isempty(xCoords) % Plot the data
                plot(xCoords, yCoords, 'o-', 'Color', colors{i}, 'MarkerFaceColor', colors{i}, 'MarkerEdgeColor', 'k', 'DisplayName', monthLabels{i});
                hold on;
            end
        end
        xticks(1:numel(sortedSites));
        xticklabels(sortedSites);
        title([num2str(currentYear)]);
        yLimits = [yLimits, ylim];
        if yearIndex == 1
            ylabel('Concentration (mg/L)');
        end
    end
    y_min = min(yLimits);
    y_max = max(yLimits);
    for yearIndex = 1:numel(uniqueYears)
        subplot(1, numel(uniqueYears), yearIndex);
        ylim([y_min, y_max]);
        if yearIndex == 2
            ax = gca;
            ax.YAxisLocation = 'right';
            set(gca, 'YTickLabel', []); % Remove numeric tick labels on the right y-axis
            set(get(gca, 'YLabel'), 'String', ''); % Remove the y-axis label
        end
    end
    hold off;
    sInput = columnNames{dPosition};
    sInputChar = char(sInput);
    formatFileName = "spatial_ls_%s.svg";
    fileName = sprintf(formatFileName, sInputChar);
    fullFilePath = fullfile(folderName, fileName);
    saveas(gcf, fullFilePath, 'svg');
end
%% Little Susitna subplots
% Extract month from dates
columnNames = vColumnLabels;
folderName = 'U:/GoA plots/NewPlots';
years = year(mLSDates); % Extract years from dates

% Define the months for which to create lines
targetMonths = [5, 6, 7, 8, 9];
colors = {[0.8, 0.6, 0.7], [0.9, 0.6, 0], [0.95, 0.9, 0.25], [0, 0.6, 0.5], [0.35, 0.7, 0.9]};
monthLabels = {'May', 'June', 'July', 'August', 'September'};

% Define the dPosition values to iterate over
dPositions = [42, 29, 11, 60];

% Create a figure with 8 subplots (4 rows x 2 columns)
figure;
% Set the figure size: [left, bottom, width, height]
set(gcf, 'Position', [100, 100, 600, 800]); % Width = 600, Height = 800

for dpIndex = 1:numel(dPositions)
    dPosition = dPositions(dpIndex);
    sInput = columnNames{dPosition}; % Get the label for the current dPosition
    yLimits = []; % Initialize an array to store y-axis limits of all subplots

    for yearIndex = 1:2 % 1 for 2022, 2 for 2023
        if yearIndex == 1
            currentYear = 2022;
        else
            currentYear = 2023;
        end
        
        % Filter data for the current year and dPosition
        currentData = mLSData(years == currentYear, :);
        currentSites = mLSSites(years == currentYear);
        months = month(mLSDates(years == currentYear));

        % Get unique sites and sort them
        uniqueSites = unique(currentSites); % Get unique site names and sort them alphabetically
        sortedSites = sort(uniqueSites);

        % Calculate the subplot index for 4 rows x 2 columns layout
        subplotIndex = (dpIndex - 1) * 2 + yearIndex;
        subplot(4, 2, subplotIndex); % Create subplot for current year and dPosition

        % Plot data for the current dPosition and year
        hold on;
        for i = 1:numel(targetMonths)
            xCoords = [];
            yCoords = [];
            for j = 1:numel(sortedSites)
                currentSite = sortedSites{j};
                % Filter data for the current site and month
                monthData = currentData(strcmp(currentSites, currentSite) & months == targetMonths(i), dPosition);
                % Check if monthData is empty or contains only NaN values
                if ~isempty(monthData) && any(~isnan(monthData))
                    nonNaNIndices = ~isnan(monthData); % Store x and y coordinates, skipping NaN values
                    xCoords = [xCoords, repmat(j, 1, sum(nonNaNIndices))];
                    yCoords = [yCoords, monthData(nonNaNIndices)];
                end
            end
            if ~isempty(xCoords) % Plot the data
                plot(xCoords, yCoords, 'o-', 'Color', colors{i}, 'MarkerFaceColor', colors{i}, 'MarkerEdgeColor', 'k', 'DisplayName', monthLabels{i});
            end
        end
        
        xticks(1:numel(sortedSites));
        xticklabels(sortedSites);
        
        if contains(sInput, 'HCO3')
            sInput = strrep(sInput, 'HCO3', 'HCO_3^-');
        end

        if contains(sInput, 'Ca') && ~contains(sInput, 'Ca^{+2}')
            sInput = strrep(sInput, 'Ca', 'Ca^{+2}');
        end

        title([sInput ' ' num2str(currentYear)], 'Interpreter', 'tex');


        % Get the y-axis limits of the current subplot and store them
        yLimits = [yLimits, ylim];
        
        % Add y-axis label and tick labels only to the first column subplots
        if yearIndex == 1
            ylabel('Concentration (mg/L)');
        end
        if dpIndex == 1 && yearIndex == 1
            legendHandle = legend('Location', 'northeast');
            set(legendHandle, 'FontSize',7);
            set(legendHandle, 'ItemTokenSize', [5,5]);
        end
    end

    % Set the same y-axis limits for all subplots in this pair
    y_min = min(yLimits);
    y_max = max(yLimits);
    for yearIndex = 1:2
        subplotIndex = (dpIndex - 1) * 2 + yearIndex;
        subplot(4, 2, subplotIndex);
        ylim([y_min, y_max]);
        % Add right y-axis labels for the second column of each row
        if yearIndex == 2
            ax = gca;
            ax.YAxisLocation = 'right';
            set(gca, 'YTickLabel', []); % Remove numeric tick labels on the right y-axis
            set(get(gca, 'YLabel'), 'String', ''); % Remove the y-axis label
        end
        % Add black outline around each subplot with thinner line width
        box on;
        set(gca, 'LineWidth', 0.75, 'Box', 'on'); % Decrease line width to 0.75
    end

    hold off;
end
% Save the figure
formatFileName = "spatial_ls_subplot.svg";
fileName = sprintf(formatFileName);
fullFilePath = fullfile(folderName, fileName);
saveas(gcf, fullFilePath, 'svg');
%% Little Susitna subplots, heavy metals
% Extract month from dates
columnNames = vColumnLabels;
folderName = 'U:/GoA plots/NewPlots';
years = year(mLSDates); % Extract years from dates

% Define the months for which to create lines
targetMonths = [5, 6, 7, 8, 9];
colors = {[0.8, 0.6, 0.7], [0.9, 0.6, 0], [0.95, 0.9, 0.25], [0, 0.6, 0.5], [0.35, 0.7, 0.9]};
monthLabels = {'May', 'June', 'July', 'August', 'September'};

% Define the dPosition values to iterate over
dPositions = [32, 37, 33, 36, 38, 63, 45, 42, 35];

% Create a figure with 8 subplots (4 rows x 2 columns)
figure;
% Set the figure size: [left, bottom, width, height]
set(gcf, 'Position', [100, 100, 1200, 6000]); % Width = 600, Height = 800

for dpIndex = 1:numel(dPositions)
    dPosition = dPositions(dpIndex);
    sInput = columnNames{dPosition}; % Get the label for the current dPosition
    yLimits = []; % Initialize an array to store y-axis limits of all subplots

    for yearIndex = 1:2 % 1 for 2022, 2 for 2023
        if yearIndex == 1
            currentYear = 2022;
        else
            currentYear = 2023;
        end
        
        % Filter data for the current year and dPosition
        currentData = mLSData(years == currentYear, :);
        currentSites = mLSSites(years == currentYear);
        months = month(mLSDates(years == currentYear));

        % Get unique sites and sort them
        uniqueSites = unique(currentSites); % Get unique site names and sort them alphabetically
        sortedSites = sort(uniqueSites);

        % Calculate the subplot index for 4 rows x 2 columns layout
        subplotIndex = (dpIndex - 1) * 2 + yearIndex;
        if subplotIndex > 18
            break;
        end
        subplot(5, 4, subplotIndex); % Create subplot for current year and dPosition

        % Plot data for the current dPosition and year
        hold on;
        for i = 1:numel(targetMonths)
            xCoords = [];
            yCoords = [];
            for j = 1:numel(sortedSites)
                currentSite = sortedSites{j};
                % Filter data for the current site and month
                monthData = currentData(strcmp(currentSites, currentSite) & months == targetMonths(i), dPosition);
                % Check if monthData is empty or contains only NaN values
                if ~isempty(monthData) && any(~isnan(monthData))
                    nonNaNIndices = ~isnan(monthData); % Store x and y coordinates, skipping NaN values
                    xCoords = [xCoords, repmat(j, 1, sum(nonNaNIndices))];
                    yCoords = [yCoords, monthData(nonNaNIndices)];
                end
            end
            if ~isempty(xCoords) % Plot the data
                plot(xCoords, yCoords, 'o-', 'Color', colors{i}, 'MarkerFaceColor', colors{i}, 'MarkerEdgeColor', 'k', 'DisplayName', monthLabels{i});
            end
        end
        
        xticks(1:numel(sortedSites));
        xticklabels(sortedSites);
        if contains(sInput, 'HCO3')
            sInput = strrep(sInput, 'HCO3', 'HCO_3');
        end
        title([sInput ' ' num2str(currentYear)], 'Interpreter', 'tex');

        % Get the y-axis limits of the current subplot and store them
        yLimits = [yLimits, ylim];
        
        % Add y-axis label and tick labels only to the first column subplots
        if yearIndex == 1
            ylabel('Concentration (mg/L)');
        end
        if dpIndex == 3 && yearIndex == 1
            legendHandle = legend('Location', 'northwest');
            set(legendHandle, 'FontSize',7);
            set(legendHandle, 'ItemTokenSize', [5,5]);
        end
    end

    % Set the same y-axis limits for all subplots in this pair
    y_min = min(yLimits);
    y_max = max(yLimits);
    for yearIndex = 1:2
        subplotIndex = (dpIndex - 1) * 2 + yearIndex;
        subplot(5, 4, subplotIndex);
        ylim([y_min, y_max]);
        % Add right y-axis labels for the second column of each row
        if yearIndex == 2
            ax = gca;
            ax.YAxisLocation = 'right';
            set(gca, 'YTickLabel', []); % Remove numeric tick labels on the right y-axis
            set(get(gca, 'YLabel'), 'String', ''); % Remove the y-axis label
        end
        % Add black outline around each subplot with thinner line width
        box on;
        set(gca, 'LineWidth', 0.75, 'Box', 'on'); % Decrease line width to 0.75
    end

    hold off;
end
% Save the figure
formatFileName = "spatial_ls_heavymetals.jpg";
fileName = sprintf(formatFileName);
fullFilePath = fullfile(folderName, fileName);
saveas(gcf, fullFilePath, 'jpg');
%% Little Susitna subplots, REE
% Extract month from dates
columnNames = vColumnLabels;
folderName = 'U:/GoA plots/NewPlots';
years = year(mLSDates); % Extract years from dates

% Define the months for which to create lines
targetMonths = [5, 6, 7, 8, 9];
colors = {[0.8, 0.6, 0.7], [0.9, 0.6, 0], [0.95, 0.9, 0.25], [0, 0.6, 0.5], [0.35, 0.7, 0.9]};
monthLabels = {'May', 'June', 'July', 'August', 'September'};

% Define the dPosition values to iterate over
dPositions = [30,50,57,60,61,49,56,51,55,59,41,52,53,54,58];

% Create a figure with 8 subplots (4 rows x 2 columns)
figure;
% Set the figure size: [left, bottom, width, height]
set(gcf, 'Position', [100, 100, 1200, 6000]); % Width = 600, Height = 800

for dpIndex = 1:numel(dPositions)
    dPosition = dPositions(dpIndex);
    sInput = columnNames{dPosition}; % Get the label for the current dPosition
    yLimits = []; % Initialize an array to store y-axis limits of all subplots

    for yearIndex = 1:2 % 1 for 2022, 2 for 2023
        if yearIndex == 1
            currentYear = 2022;
        else
            currentYear = 2023;
        end
        
        % Filter data for the current year and dPosition
        currentData = mLSData(years == currentYear, :);
        currentSites = mLSSites(years == currentYear);
        months = month(mLSDates(years == currentYear));

        % Get unique sites and sort them
        uniqueSites = unique(currentSites); % Get unique site names and sort them alphabetically
        sortedSites = sort(uniqueSites);

        % Calculate the subplot index for 4 rows x 2 columns layout
        subplotIndex = (dpIndex - 1) * 2 + yearIndex;
        subplot(5, 6, subplotIndex); % Create subplot for current year and dPosition

        % Plot data for the current dPosition and year
        hold on;
        for i = 1:numel(targetMonths)
            xCoords = [];
            yCoords = [];
            for j = 1:numel(sortedSites)
                currentSite = sortedSites{j};
                % Filter data for the current site and month
                monthData = currentData(strcmp(currentSites, currentSite) & months == targetMonths(i), dPosition);
                % Check if monthData is empty or contains only NaN values
                if ~isempty(monthData) && any(~isnan(monthData))
                    nonNaNIndices = ~isnan(monthData); % Store x and y coordinates, skipping NaN values
                    xCoords = [xCoords, repmat(j, 1, sum(nonNaNIndices))];
                    yCoords = [yCoords, monthData(nonNaNIndices)];
                end
            end
            if ~isempty(xCoords) % Plot the data
                plot(xCoords, yCoords, 'o-', 'Color', colors{i}, 'MarkerFaceColor', colors{i}, 'MarkerEdgeColor', 'k', 'DisplayName', monthLabels{i});
            end
        end
        
        xticks(1:numel(sortedSites));
        xticklabels(sortedSites);
        if contains(sInput, 'HCO3')
            sInput = strrep(sInput, 'HCO3', 'HCO_3');
        end
        title([sInput ' ' num2str(currentYear)], 'Interpreter', 'tex');

        % Get the y-axis limits of the current subplot and store them
        yLimits = [yLimits, ylim];
        
        % Add y-axis label and tick labels only to the first column subplots
        if yearIndex == 1
            ylabel('Concentration (mg/L)');
        end
        if dpIndex == 2 && yearIndex == 2
            legendHandle = legend('Location', 'northwest');
            set(legendHandle, 'FontSize',7);
            set(legendHandle, 'ItemTokenSize', [5,5]);
        end
    end

    % Set the same y-axis limits for all subplots in this pair
    y_min = min(yLimits);
    y_max = max(yLimits);
    for yearIndex = 1:2
        subplotIndex = (dpIndex - 1) * 2 + yearIndex;
        subplot(5, 6, subplotIndex);
        ylim([y_min, y_max]);
        % Add right y-axis labels for the second column of each row
        if yearIndex == 2
            ax = gca;
            ax.YAxisLocation = 'right';
            set(gca, 'YTickLabel', []); % Remove numeric tick labels on the right y-axis
            set(get(gca, 'YLabel'), 'String', ''); % Remove the y-axis label
        end
        % Add black outline around each subplot with thinner line width
        box on;
        set(gca, 'LineWidth', 0.75, 'Box', 'on'); % Decrease line width to 0.75
    end

    hold off;
end
% Save the figure
formatFileName = "spatial_ls_REE.svg";
fileName = sprintf(formatFileName);
fullFilePath = fullfile(folderName, fileName);
saveas(gcf, fullFilePath, 'svg');
%% Spatial figures: Gulkana Watershed
columnNames = vColumnLabels;
folderName = 'U:/GoA plots/NewPlots';
years = year(mGulkanaDates); % Extract years from dates
data2022 = mGulkanaData(years == 2022, :);
data2023 = mGulkanaData(years == 2023, :);
targetMonths = [5, 6, 7, 8, 9];
colors = {[0.8, 0.6, 0.7],[0.9, 0.6, 0],[0.95, 0.9, 0.25],[0, 0.6, 0.5],[0.35, 0.7, 0.9]};
monthLabels = {'May', 'June', 'July', 'August', 'September'};
uniqueYears = unique(years);
for dPosition = 1:numel(columnNames)
    figure;
    yLimits = []; % Initialize an array to store y-axis limits of all subplots
    for yearIndex = 1:numel(uniqueYears) % Loop for each year
        currentYear = uniqueYears(yearIndex);
        currentData = mGulkanaData(years == currentYear, :);
        months = month(mGulkanaDates(years == currentYear));
        uniqueSites = unique(mGulkanaSites); % Get unique site names and sort them alphabetically
        sortedSites = sort(uniqueSites);
        subplot(1, numel(uniqueYears), yearIndex); % Create subplot for current year
        for i = 1:numel(targetMonths)
            xCoords = [];
            yCoords = [];
            for j = 1:numel(sortedSites)
                currentSite = sortedSites{j};
                monthData = currentData(strcmp(mGulkanaSites(years == currentYear), currentSite) & months == targetMonths(i), dPosition);
                nonNaNIndices = ~isnan(monthData); % Store x and y coordinates, skipping NaN values
                xCoords = [xCoords, repmat(j, 1, sum(nonNaNIndices))];
                yCoords = [yCoords, monthData(nonNaNIndices)];
            end
            if ~isempty(xCoords) % Plot the data
                plot(xCoords, yCoords, 'o-', 'Color', colors{i}, 'MarkerFaceColor', colors{i}, 'MarkerEdgeColor', 'k', 'DisplayName', monthLabels{i});
                hold on;
            end
        end
        xticks(1:numel(sortedSites));
        xticklabels(sortedSites);
        title([num2str(currentYear)]);
        yLimits = [yLimits, ylim];
        if yearIndex == 1
            ylabel('Concentration (mg/L)');
        end
    end
    y_min = min(yLimits);
    y_max = max(yLimits);
    for yearIndex = 1:numel(uniqueYears)
        subplot(1, numel(uniqueYears), yearIndex);
        ylim([y_min, y_max]);
        if yearIndex == 2
            ax = gca;
            ax.YAxisLocation = 'right';
            set(gca, 'YTickLabel', []); % Remove numeric tick labels on the right y-axis
            set(get(gca, 'YLabel'), 'String', ''); % Remove the y-axis label
        end
    end
    hold off;
    sInput = columnNames{dPosition};
    sInputChar = char(sInput);
    formatFileName = "spatial_gulk_%s.jpg";
    fileName = sprintf(formatFileName, sInputChar);
    fullFilePath = fullfile(folderName, fileName);
    saveas(gcf, fullFilePath, 'jpg');
end
%% Spatial figures: Canwell Watershed
columnNames = vColumnLabels;
folderName = 'U:/GoA plots/NewPlots';
years = year(mCanwellDates); % Extract years from dates
data2022 = mCanwellData(years == 2022, :);
data2023 = mCanwellData(years == 2023, :);
targetMonths = [5, 6, 7, 8, 9];
colors = {[0.8, 0.6, 0.7],[0.9, 0.6, 0],[0.95, 0.9, 0.25],[0, 0.6, 0.5],[0.35, 0.7, 0.9]};
monthLabels = {'May', 'June', 'July', 'August', 'September'}; %pink, orange, yellow, green, blue
uniqueYears = unique(years);
for dPosition = 1:numel(columnNames)
    figure;
    yLimits = []; % Initialize an array to store y-axis limits of all subplots
    for yearIndex = 1:numel(uniqueYears) % Loop for each year
        currentYear = uniqueYears(yearIndex);
        currentData = mCanwellData(years == currentYear, :);
        months = month(mCanwellDates(years == currentYear));
        uniqueSites = unique(mCanwellSites); % Get unique site names and sort them alphabetically
        sortedSites = sort(uniqueSites);
        subplot(1, numel(uniqueYears), yearIndex); % Create subplot for current year
        for i = 1:numel(targetMonths)
            xCoords = [];
            yCoords = [];
            for j = 1:numel(sortedSites)
                currentSite = sortedSites{j};
                monthData = currentData(strcmp(mCanwellSites(years == currentYear), currentSite) & months == targetMonths(i), dPosition);
                nonNaNIndices = ~isnan(monthData); % Store x and y coordinates, skipping NaN values
                xCoords = [xCoords, repmat(j, 1, sum(nonNaNIndices))];
                yCoords = [yCoords, monthData(nonNaNIndices)];
            end
            if ~isempty(xCoords) % Plot the data
                plot(xCoords, yCoords, 'o-', 'Color', colors{i}, 'MarkerFaceColor', colors{i}, 'MarkerEdgeColor', 'k', 'DisplayName', monthLabels{i});
                hold on;
            end
        end
        xticks(1:numel(sortedSites));
        xticklabels(sortedSites);
        title([num2str(currentYear)]);
        yLimits = [yLimits, ylim];
        if yearIndex == 1
            ylabel('Concentration (mg/L)');
        end
    end
    y_min = min(yLimits);
    y_max = max(yLimits);
    for yearIndex = 1:numel(uniqueYears)
        subplot(1, numel(uniqueYears), yearIndex);
        ylim([y_min, y_max]);
        if yearIndex == 2
            ax = gca;
            ax.YAxisLocation = 'right';
            set(gca, 'YTickLabel', []); % Remove numeric tick labels on the right y-axis
            set(get(gca, 'YLabel'), 'String', ''); % Remove the y-axis label
        end
    end
    hold off;
    sInput = columnNames{dPosition};
    sInputChar = char(sInput);
    formatFileName = "spatial_canwell_%s.jpg";
    fileName = sprintf(formatFileName, sInputChar);
    fullFilePath = fullfile(folderName, fileName);
    saveas(gcf, fullFilePath, 'jpg');
end
%% Spatial figures: Castner Watershed
columnNames = vColumnLabels;
folderName = 'U:/GoA plots/NewPlots';
years = year(mCastnerDates); % Extract years from dates
data2022 = mCastnerData(years == 2022, :);
data2023 = mCastnerData(years == 2023, :);
targetMonths = [5, 6, 7, 8, 9];
colors = {[0.8, 0.6, 0.7],[0.9, 0.6, 0],[0.95, 0.9, 0.25],[0, 0.6, 0.5],[0.35, 0.7, 0.9]};
monthLabels = {'May', 'June', 'July', 'August', 'September'};
uniqueYears = unique(years);
for dPosition = 1:numel(columnNames)
    figure;
    yLimits = []; % Initialize an array to store y-axis limits of all subplots
    for yearIndex = 1:numel(uniqueYears) % Loop for each year
        currentYear = uniqueYears(yearIndex);
        currentData = mCastnerData(years == currentYear, :);
        months = month(mCastnerDates(years == currentYear));
        uniqueSites = unique(mCastnerSites); % Get unique site names and sort them alphabetically
        sortedSites = sort(uniqueSites);
        subplot(1, numel(uniqueYears), yearIndex); % Create subplot for current year
        for i = 1:numel(targetMonths)
            xCoords = [];
            yCoords = [];
            for j = 1:numel(sortedSites)
                currentSite = sortedSites{j};
                monthData = currentData(strcmp(mCastnerSites(years == currentYear), currentSite) & months == targetMonths(i), dPosition);
                nonNaNIndices = ~isnan(monthData); % Store x and y coordinates, skipping NaN values
                xCoords = [xCoords, repmat(j, 1, sum(nonNaNIndices))];
                yCoords = [yCoords, monthData(nonNaNIndices)];
            end
            if ~isempty(xCoords) % Plot the data
                plot(xCoords, yCoords, 'o-', 'Color', colors{i}, 'MarkerFaceColor', colors{i}, 'MarkerEdgeColor', 'k', 'DisplayName', monthLabels{i});
                hold on;
            end
        end
        xticks(1:numel(sortedSites));
        xticklabels(sortedSites);
        title([num2str(currentYear)]);
        yLimits = [yLimits, ylim];
        if yearIndex == 1
            ylabel('Concentration (mg/L)');
        end
    end
    y_min = min(yLimits);
    y_max = max(yLimits);
    for yearIndex = 1:numel(uniqueYears)
        subplot(1, numel(uniqueYears), yearIndex);
        ylim([y_min, y_max]);
        if yearIndex == 2
            ax = gca;
            ax.YAxisLocation = 'right';
            set(gca, 'YTickLabel', []); % Remove numeric tick labels on the right y-axis
            set(get(gca, 'YLabel'), 'String', ''); % Remove the y-axis label
        end
    end
    hold off;
    sInput = columnNames{dPosition};
    sInputChar = char(sInput);
    formatFileName = "spatial_castner_%s.jpg";
    fileName = sprintf(formatFileName, sInputChar);
    fullFilePath = fullfile(folderName, fileName);
    saveas(gcf, fullFilePath, 'jpg');
end
%% %% Stable isotope comparison: USGS - done
figure;
plot(mKR3Data(:,9), mKR3Data(:,10), 'o', 'DisplayName', 'KR3','Color','k','MarkerFaceColor',[0.35,0.70,0.90]);
hold on;
plot(mMR4Data(:,9), mMR4Data(:,10), 'o', 'DisplayName', 'MR4','Color','k','MarkerFaceColor',[0.95,0.90,0.25]);
plot(mMCData(:,9), mMCData(:,10), 'o', 'DisplayName', 'Moose Creek','Color','k','MarkerFaceColor',[0.90,0.60,0]);
plot(mLS2Data(:,9), mLS2Data(:,10),'o', 'DisplayName', 'LS2','Color','k','MarkerFaceColor',[0,0.60,0.50]);
xlabel('DO'); 
ylabel('DD');  
title('USGS stable isotope comparison');
legend('Location', 'Best');
grid on;
hold off;
folderName = 'U:/GoA plots/NewPlots';
fileName = 'isotopes_usgs_no_title.svg';
fullFilePath = fullfile(folderName, fileName);
saveas(gcf, fullFilePath, 'svg');
%% Stable isotope comparison: all - done
figure;
plot(mKR1Data(:,9), mKR1Data(:,10), 'd', 'DisplayName', 'KR1','Color','k','MarkerFaceColor',[0.35,0.70,0.90]);
hold on;
plot(mKR2Data(:,9), mKR2Data(:,10), 's', 'DisplayName', 'KR2','Color','k','MarkerFaceColor',[0.35,0.70,0.90]);
plot(mKR3Data(:,9), mKR3Data(:,10), 'o', 'DisplayName', 'KR3','Color','k','MarkerFaceColor',[0.35,0.70,0.90]);
plot(mKR4Data(:,9), mKR4Data(:,10),'^', 'DisplayName', 'KR4','Color','k','MarkerFaceColor',[0.35,0.70,0.90]);
plot(mMR1Data(:,9), mMR1Data(:,10), '^', 'DisplayName', 'MR1','Color','k','MarkerFaceColor',[0.95,0.90,0.25]);
plot(mMR2Data(:,9), mMR2Data(:,10), 's', 'DisplayName', 'MR2','Color','k','MarkerFaceColor',[0.95,0.90,0.25]);
plot(mMR3Data(:,9), mMR3Data(:,10), 'd', 'DisplayName', 'MR3','Color','k','MarkerFaceColor',[0.95,0.90,0.25]);
plot(mMR4Data(:,9), mMR4Data(:,10), 'o', 'DisplayName', 'MR4','Color','k','MarkerFaceColor',[0.95,0.90,0.25]);
plot(mMR5Data(:,9), mMR5Data(:,10),'pentagram', 'DisplayName', 'MR5','Color','k','MarkerFaceColor',[0.95,0.90,0.25]);
plot(mMCData(:,9), mMCData(:,10), 'o', 'DisplayName', 'Moose Creek','Color','k','MarkerFaceColor',[0.90,0.60,0]);
plot(mLS1Data(:,9), mLS1Data(:,10), '^', 'DisplayName', 'LS1','Color','k','MarkerFaceColor',[0,0.60,0.50])
plot(mLS15Data(:,9), mLS15Data(:,10),'s', 'DisplayName', 'LS1.5','Color','k','MarkerFaceColor',[0,0.60,0.50])
plot(mLS2Data(:,9), mLS2Data(:,10), 'o', 'DisplayName', 'LS2','Color','k','MarkerFaceColor',[0,0.60,0.50])
plot(mLS3Data(:,9), mLS3Data(:,10),'d', 'DisplayName', 'LS3','Color','k','MarkerFaceColor',[0,0.60,0.50])
plot(mLS4Data(:,9), mLS4Data(:,10), 'pentagram', 'DisplayName', 'LS4','Color','k','MarkerFaceColor',[0,0.60,0.50])
plot(mG1Data(:,9), mG1Data(:,10),'o', 'DisplayName', 'G1','Color','k','MarkerFaceColor',[0.80,0.40,0]);
plot(mG2Data(:,9), mG2Data(:,10),'s', 'DisplayName', 'G2','Color','k','MarkerFaceColor',[0.80,0.40,0]);
plot(mG3Data(:,9), mG3Data(:,10),'d', 'DisplayName', 'G3','Color','k','MarkerFaceColor',[0.80,0.40,0]);
plot(mCT1Data(:,9), mCT1Data(:,10),'o', 'DisplayName', 'CT1','Color','k','MarkerFaceColor',[0,0.45,0.70]);
plot(mCT2Data(:,9), mCT2Data(:,10),'s', 'DisplayName', 'CT2','Color','k','MarkerFaceColor',[0,0.45,0.70]);
plot(mCW1Data(:,9), mCW1Data(:,10),'o', 'DisplayName', 'CW1','Color','k','MarkerFaceColor',[0.80,0.60,0.70]);
plot(mCW2Data(:,9), mCW2Data(:,10),'s', 'DisplayName', 'CW2','Color','k','MarkerFaceColor',[0.80,0.60,0.70]);
xlabel('Average DO'); 
ylabel('Average DD');  
title('Stable isotope comparison');
legend('Location', 'Best');
grid on;
legend('Location', 'eastoutside', 'Orientation', 'vertical');
hold off;
folderName = 'U:/GoA plots/NewPlots';
fileName = 'isotope_comparison_all.jpg';
fullFilePath = fullfile(folderName, fileName);
saveas(gcf, fullFilePath, 'jpg');
%% Stable Isotopes over time: USGS - done
figure;
plot(mKR3_Date22, mKR3_2022(:,9), 'o-', 'DisplayName', 'KR3 average DO','Color',[0.35,0.70,0.90],'MarkerFaceColor',[0.35,0.70,0.90],'MarkerEdgeColor','k');
hold on;
plot(mKR3_Date23, mKR3_2023(:,9), 'o-', 'DisplayName', 'KR3 average DO','Color',[0.35,0.70,0.90],'MarkerFaceColor',[0.35,0.70,0.90],'MarkerEdgeColor','k');
plot(mMR4_Date22, mMR4_2022(:,9), 'o-', 'DisplayName', 'MR4 average DO','Color',[0.95,0.90,0.25],'MarkerFaceColor',[0.95,0.90,0.25],'MarkerEdgeColor','k');
plot(mMR4_Date23, mMR4_2023(:,9), 'o-', 'DisplayName', 'MR4 average DO','Color',[0.95,0.90,0.25],'MarkerFaceColor',[0.95,0.90,0.25],'MarkerEdgeColor','k');
plot(mMC_Date22, mMC_2022(:,9), 'o-', 'DisplayName', 'MC average DO','Color',[0.90,0.60,0],'MarkerFaceColor',[0.90,0.60,0],'MarkerEdgeColor','k');
plot(mMC_Date23, mMC_2023(:,9), 'o-', 'DisplayName', 'MC average DO','Color',[0.90,0.60,0],'MarkerFaceColor',[0.90,0.60,0],'MarkerEdgeColor','k');
plot(mLS2_Date22, mLS2_2022(:,9),'o-', 'DisplayName', 'LS2 average DO','Color',[0,0.60,0.50],'MarkerFaceColor',[0,0.60,0.50],'MarkerEdgeColor','k');
plot(mLS2_Date23, mLS2_2023(:,9),'o-', 'DisplayName', 'LS2 average DO','Color',[0,0.60,0.50],'MarkerFaceColor',[0,0.60,0.50],'MarkerEdgeColor','k');
plot(mKR3_Date22, mKR3_2022(:,10), '--o', 'DisplayName', 'KR3 average DD','Color',[0.35,0.70,0.90],'MarkerFaceColor',[0.35,0.70,0.90],'MarkerEdgeColor','k');
plot(mKR3_Date23, mKR3_2023(:,10), '--o', 'DisplayName', 'KR3 average DD','Color',[0.35,0.70,0.90],'MarkerFaceColor',[0.35,0.70,0.90],'MarkerEdgeColor','k');
plot(mMR4_Date22, mMR4_2022(:,10), '--o', 'DisplayName', 'MR4 average DD','Color',[0.95,0.90,0.25],'MarkerFaceColor',[0.95,0.90,0.25],'MarkerEdgeColor','k');
plot(mMR4_Date23, mMR4_2023(:,10), '--o', 'DisplayName', 'MR4 average DD','Color',[0.95,0.90,0.25],'MarkerFaceColor',[0.95,0.90,0.25],'MarkerEdgeColor','k');
plot(mMC_Date22, mMC_2022(:,10), '--o', 'DisplayName', 'MC average DD','Color',[0.90,0.60,0],'MarkerFaceColor',[0.90,0.60,0],'MarkerEdgeColor','k');
plot(mMC_Date23, mMC_2023(:,10), '--o', 'DisplayName', 'MC average DD','Color',[0.90,0.60,0],'MarkerFaceColor',[0.90,0.60,0],'MarkerEdgeColor','k');
plot(mLS2_Date22, mLS2_2022(:,10),'--o', 'DisplayName', 'LS2 average DD','Color',[0,0.60,0.50],'MarkerFaceColor',[0,0.60,0.50],'MarkerEdgeColor','k');
plot(mLS2_Date23, mLS2_2023(:,10),'--o', 'DisplayName', 'LS2 average DD','Color',[0,0.60,0.50],'MarkerFaceColor',[0,0.60,0.50],'MarkerEdgeColor','k');
xlabel('Time');  
ylabel('Stable isotope value');  
title('Stable isotope variations over time');
legend('Location', 'eastoutside');
grid on;
hold off;
folderName = 'U:/GoA plots/NewPlots';
fileName = 'isotopes_time_usgs.jpg';
fullFilePath = fullfile(folderName, fileName);
saveas(gcf, fullFilePath, 'jpg');

%% Stable isotope comparison: watershed - looks good
figure;
ms=7;
msS=10;
% Define the range for DO (δ¹⁸O) based on your x-axis range
DO_range = linspace(-26,-17, 100); 
% Compute DD (δ²H) using the meteoric water line equation
DD_meteoric = 8 * DO_range + 10;
% Plot the meteoric water line
plot(DO_range, DD_meteoric, 'k--', 'DisplayName', 'GMWL', 'LineWidth', 0.5);
hold on;
plot(mKnikMain(:,9), mKnikMain(:,10), 'o', 'DisplayName', 'Knik','MarkerEdgeColor', [0.35,0.70,0.90],'MarkerSize', ms, 'LineWidth', 2.5);
plot(mKnikSpring(:,9), mKnikSpring(:,10), 'kd', 'DisplayName', 'Knik Springs','MarkerFaceColor', [0.35,0.70,0.90],'MarkerSize', ms, 'LineWidth', 0.5);
plot(mKnikTrib(:,9), mKnikTrib(:,10), 'ks', 'DisplayName', 'Knik Tributaries','MarkerFaceColor', [0.35,0.70,0.90],'MarkerSize', msS, 'LineWidth', 0.5);
plot(mMatMain(:,9), mMatMain(:,10), 'o', 'DisplayName', 'Matanuska','MarkerEdgeColor', [0.95,0.90,0.25],'MarkerSize', ms, 'LineWidth', 2.5);
plot(mMatSpring(:,9), mMatSpring(:,10), 'kd', 'DisplayName', 'Matanuska Springs','MarkerFaceColor', [0.95,0.90,0.25],'MarkerSize', ms, 'LineWidth', 0.5);
plot(mMatTrib(:,9), mMatTrib(:,10), 'ks', 'DisplayName', 'Matanuska Tributaries','MarkerFaceColor', [0.95,0.90,0.25],'MarkerSize', msS, 'LineWidth', 0.5);
plot(mMCData(:,9), mMCData(:,10), 'o', 'DisplayName', 'Moose Creek','MarkerEdgeColor', [0.90,0.60,0],'MarkerSize', ms, 'LineWidth', 0.5);
plot(mLSMain(:,9), mLSMain(:,10),'o', 'DisplayName', 'Little Susitna','MarkerEdgeColor', [0,0.60,0.50],'MarkerSize', ms, 'LineWidth', 2.5);
plot(mLSSpring(:,9), mLSSpring(:,10),'kd', 'DisplayName', 'Little Susitna Springs','MarkerFaceColor', [0,0.60,0.50],'MarkerSize', ms, 'LineWidth', 0.5);
plot(mLSTrib(:,9), mLSTrib(:,10),'ks', 'DisplayName', 'Little Susitna Tributaries','MarkerFaceColor', [0,0.60,0.50],'MarkerSize', msS, 'LineWidth', 0.5);
plot(mLSSub(:,9), mLSSub(:,10),'kp', 'DisplayName', 'Little Susitna Subglacial','MarkerFaceColor', [0,0.60,0.50],'MarkerSize', msS, 'LineWidth', 0.5);
plot(mLSPeri(:,9), mLSPeri(:,10),'k^', 'DisplayName', 'Little Susitna Periglacial','MarkerFaceColor', [0,0.60,0.50],'MarkerSize', ms, 'LineWidth', 0.5);
plot(mLSSupra(:,9), mLSSupra(:,10),'kv', 'DisplayName', 'Little Susitna Supraglacial','MarkerFaceColor', [0,0.60,0.50],'MarkerSize', ms, 'LineWidth', 0.5);
plot(mLSLake(:,9), mLSLake(:,10),'k>', 'DisplayName', 'Little Susitna Glacial Lake','MarkerFaceColor', [0,0.60,0.50],'MarkerSize', ms, 'LineWidth', 0.5);
plot(mGulkanaMain(:,9), mGulkanaMain(:,10), 'o', 'DisplayName', 'Gulkana','MarkerEdgeColor', [0.80,0.40,0],'MarkerSize', ms, 'LineWidth', 2.5);
plot(mGulkanaSupra(:,9), mGulkanaSupra(:,10), 'kv', 'DisplayName', 'Gulkana Supraglacial','MarkerFaceColor', [0.80,0.40,0],'MarkerSize', ms, 'LineWidth', 0.5);
plot(mGulkanaTrib(:,9), mGulkanaTrib(:,10), 'ks', 'DisplayName', 'Gulkana Tributaries','MarkerFaceColor', [0.80,0.40,0],'MarkerSize', msS, 'LineWidth', 0.5);
plot(mGulkanaSub(:,9), mGulkanaSub(:,10), 'kp', 'DisplayName', 'Gulkana Subglacial','MarkerFaceColor', [0.80,0.40,0],'MarkerSize', msS, 'LineWidth', 0.5);
plot(mCastnerMain(:,9), mCastnerMain(:,10),'o', 'DisplayName', 'Castner','MarkerEdgeColor', [0,0.45,0.70],'MarkerSize', ms, 'LineWidth', 2.5);
plot(mCastnerSub(:,9), mCastnerSub(:,10),'kp', 'DisplayName', 'Castner Subglacial','MarkerFaceColor', [0,0.45,0.70],'MarkerSize', msS, 'LineWidth', 0.5);
plot(mCanwellData(:,9), mCanwellData(:,10), 'o', 'DisplayName', 'Canwell','MarkerEdgeColor', [0.80,0.60,0.70],'MarkerSize', ms, 'LineWidth', 2.5);

xlabel('δ¹⁸O (‰)', 'FontSize', 20); 
ylabel('δ²H (‰)', 'FontSize', 20);  
%title('Stable isotope comparison');
ax = gca;
ax.FontSize = 15; 
legendHandle = legend('Location', 'northwest');
set(legendHandle, 'FontSize',11, 'NumColumns', 2);
set(legendHandle, 'ItemTokenSize', [20,20]);
grid off;
hold off;
set(gcf, 'Position', [100, 100, 1400, 1000]);
folderName = 'U:/GoA plots/NewPlots';
fileName = 'isotope_comparison_watershed.svg';
fullFilePath = fullfile(folderName, fileName);
saveas(gcf, fullFilePath, 'svg');

%% Stable isotope comparison: watershed - looks good
figure;
ms=7;
m=8;
% Define the range for DO (δ¹⁸O) based on your x-axis range
DO_range = linspace(-26,-17, 100); 
% Compute DD (δ²H) using the meteoric water line equation
DD_meteoric = 8 * DO_range + 10;
% Plot the meteoric water line
plot(DO_range, DD_meteoric, 'k--', 'DisplayName', 'GMWL', 'LineWidth', 0.5,'HandleVisibility', 'off');
hold on;
plot(mKnikMain(:,9), mKnikMain(:,10), '^', 'DisplayName', 'Knik', 'Color', 'k', 'MarkerFaceColor', [0.35,0.70,0.90],'MarkerSize',m);
plot(mMatMain(:,9), mMatMain(:,10), 's', 'DisplayName', 'Matanuska', 'Color', 'k', 'MarkerFaceColor', [0.95,0.90,0.25],'MarkerSize',m);
plot(mLSMain(:,9), mLSMain(:,10), 'd', 'DisplayName', 'Little Susitna', 'Color', 'k', 'MarkerFaceColor', [0,0.60,0.50],'MarkerSize',m);
plot(mMCData(:,9), mMCData(:,10), '>', 'DisplayName', 'Moose Creek', 'Color', 'k', 'MarkerFaceColor', [0.90,0.60,0],'MarkerSize',m);
plot(mCastnerMain(:,9), mCastnerMain(:,10), 'v', 'DisplayName', 'Castner', 'Color', 'k', 'MarkerFaceColor', [0,0.45,0.70],'MarkerSize',m);
plot(mCanwellData(:,9), mCanwellData(:,10), 'pentagram', 'DisplayName', 'Canwell', 'Color', 'k', 'MarkerFaceColor', [0.80,0.60,0.70],'MarkerSize',m);
plot(mGulkanaMain(:,9), mGulkanaMain(:,10), 'o', 'DisplayName', 'Gulkana', 'Color', 'k', 'MarkerFaceColor', [0.80,0.40,0],'MarkerSize',m);

xlabel('δ¹⁸O (‰)', 'FontSize', 20); 
ylabel('δ²H (‰)', 'FontSize', 20);  
%title('Stable isotope comparison');
ax = gca;
ax.FontSize = 15; 
legendHandle = legend('Location', 'northwest');
set(legendHandle, 'FontSize',11);
set(legendHandle, 'ItemTokenSize', [30,20]);
grid off;
hold off;
set(gcf, 'Position', [100, 100, 1000, 800]);
folderName = 'U:/GoA plots/NewPlots';
fileName = 'isotope_comparison_watershed_main_working.jpg';
fullFilePath = fullfile(folderName, fileName);
saveas(gcf, fullFilePath, 'jpg');
%% Subplots by month, isotopic data dH v dO
% Define months of interest
months = {'May', 'June', 'July', 'August', 'September/October'};

% Define the date range for each month (adjust the exact date ranges accordingly)
monthRanges = [5, 6, 7, 8, 9]; % Representing May, June, etc.

% Define the year of interest
years = [2022, 2023];

% Create figure for 5 rows (months) and 2 columns (years)
figure;
ms=7;
for i = 1:length(monthRanges)
    for j = 1:length(years)
        % Filter data based on year and month using logical indexing
        currentYear = years(j);
        currentMonth = monthRanges(i);
        
        % Filter by year and month for each dataset
        % For instance, filtering mMatMain data
        idxMat = year(mMatDates) == currentYear & month(mMatDates) == currentMonth;
        idxCanwell = year(mCanwellDates) == currentYear & month(mCanwellDates) == currentMonth;
        idxCastner = year(mCastnerDates) == currentYear & month(mCanwellDates) == currentMonth;
        idxGulkana = year(mGulkanaDates) == currentYear & month(mGulkanaDates) == currentMonth;
        idxMC = year(mMCDates) == currentYear & month(mMCDates) == currentMonth;
        idxKnik = year(mKnikDates) == currentYear & month(mKnikDates) == currentMonth;
        idxLS = year(mLSDates) == currentYear & month(mLSDates) == currentMonth;

        % Subplot for the corresponding month and year
        subplot(5, 2, (i-1)*2 + j);
        
        % Plot the meteoric water line
        DO_range = linspace(-26, -17, 100);
        DD_meteoric = 8 * DO_range + 10;
        plot(DO_range, DD_meteoric, 'k--', 'DisplayName', 'GMWL', 'LineWidth', 0.5);
        hold on;
        
        % Plot filtered data
        plot(mKnikData(idxKnik, 9), mKnikData(idxKnik, 10), 'ko', 'DisplayName', 'Knik', 'MarkerFaceColor', [0.35,0.70,0.90], 'MarkerSize', ms, 'LineWidth', 0.5);
        plot(mMatData(idxMat, 9), mMatData(idxMat, 10), 'ko', 'DisplayName', 'Matanuska', 'MarkerFaceColor', [0.95,0.90,0.25], 'MarkerSize', ms, 'LineWidth', 0.5);
        plot(mMCData(idxMC, 9), mMCData(idxMC, 10), 'ko', 'DisplayName', 'Moose Creek', 'MarkerFaceColor', [0.90,0.60,0], 'MarkerSize', ms, 'LineWidth', 0.5);
        plot(mLSData(idxLS, 9), mLSData(idxLS, 10), 'ko', 'DisplayName', 'Little Susitna', 'MarkerFaceColor', [0,0.60,0.50], 'MarkerSize', ms, 'LineWidth', 0.5);
        plot(mGulkanaData(idxGulkana, 9), mGulkanaData(idxGulkana, 10), 'ko', 'DisplayName', 'Gulkana', 'MarkerFaceColor', [0.80,0.40,0], 'MarkerSize', ms, 'LineWidth', 0.5);
        plot(mCastnerData(idxCastner, 9), mCastnerData(idxCastner, 10), 'ko', 'DisplayName', 'Castner', )
        plot(mCanwellData(idxCanwell, 9), mCanwellData(idxCanwell, 10), 'ko', 'DisplayName', 'Canwell', 'MarkerFaceColor', [0.80,0.60,0.70], 'MarkerSize', ms, 'LineWidth', 0.5);
        
        % Set labels and other plot details
        xlabel('δ¹⁸O (‰)', 'FontSize', 12); 
        ylabel('δ²H (‰)', 'FontSize', 12);  
        title([months{i}, ' ', num2str(currentYear)]);
        
        % Adjust legend and grid settings
        legendHandle = legend('Location', 'northwest');
        set(legendHandle, 'FontSize', 9);
        grid off;
        hold off;
    end
end

% Adjust figure size
set(gcf, 'Position', [100, 100, 1200, 1000]);

% Save the figure
folderName = 'U:/GoA plots/NewPlots';
fileName = 'isotope_comparison_month.svg';
fullFilePath = fullfile(folderName, fileName);
saveas(gcf, fullFilePath, 'svg');

%% Median and Averages for isotope endmembers
avgKMO = mean(mKnikMain(:,9), 'omitnan');
avgKMH = mean(mKnikMain(:,10), 'omitnan');
avgKSO = mean(mKnikSpring(:,9), 'omitnan');
avgKSH = mean(mKnikSpring(:,10), 'omitnan');
avgKTO = mean(mKnikTrib(:,9), 'omitnan');
avgKTH = mean(mKnikTrib(:,10), 'omitnan');
avgMMO = mean(mMatMain(:,9), 'omitnan');
avgMMH = mean(mMatMain(:,10), 'omitnan');
avgMSO = mean(mMatSpring(:,9), 'omitnan');
avgMSH = mean(mMatSpring(:,10), 'omitnan');
avgMTO = mean(mMatTrib(:,9), 'omitnan');
avgMTH = mean(mMatTrib(:,10), 'omitnan');
avgMCO = mean(mMCData(:,9), 'omitnan');
avgMCH = mean(mMCData(:,10), 'omitnan');
avgLSMO = mean(mLSMain(:,9), 'omitnan');
avgLSMH = mean(mLSMain(:,10), 'omitnan');
avgLSSO = mean(mLSSpring(:,9), 'omitnan');
avgLSSH = mean(mLSSpring(:,10), 'omitnan');
avgLSTO = mean(mLSTrib(:,9), 'omitnan');
avgLSTH = mean(mLSTrib(:,10), 'omitnan');
avgLSSubO = mean(mLSSub(:,9), 'omitnan');
avgLSSubH = mean(mLSSub(:,10), 'omitnan');
avgLSPeriO = mean(mLSPeri(:,9), 'omitnan');
avgLSPeriH = mean(mLSPeri(:,10), 'omitnan');
avgLSSupraO = mean(mLSSupra(:,9), 'omitnan');
avgLSSupraH = mean(mLSSupra(:,10), 'omitnan');
avgLSLO = mean(mLSLake(:,9), 'omitnan');
avgLSLH = mean(mLSLake(:,10), 'omitnan');
avgGMO = mean(mGulkanaMain(:,9), 'omitnan');
avgGMH = mean(mGulkanaMain(:,10), 'omitnan');
avgGTO = mean(mGulkanaTrib(:,9), 'omitnan');
avgGTH = mean(mGulkanaTrib(:,10), 'omitnan');
avgGSubO = mean(mGulkanaSub(:,9), 'omitnan');
avgGSubH = mean(mGulkanaSub(:,10), 'omitnan');
avgCTMO = mean(mCastnerMain(:,9), 'omitnan');
avgCTMH = mean(mCastnerMain(:,10), 'omitnan');
avgCTSubO = mean(mCastnerSub(:,9), 'omitnan');
avgCTSubH = mean(mCastnerSub(:,10), 'omitnan');
avgCWMO = mean(mCanwellData(:,9), 'omitnan');
avgCWMH = mean(mCanwellData(:,10), 'omitnan');
%% Stable isotope comparison: watershed - looks good
figure;
m=14;
% Define the range for DO (δ¹⁸O) based on your x-axis range
DO_range = linspace(-26,-17, 100); 
% Compute DD (δ²H) using the meteoric water line equation
DD_meteoric = 8 * DO_range + 10;
% Plot the meteoric water line
plot(DO_range, DD_meteoric, 'k--', 'DisplayName', 'GMWL', 'LineWidth', 0.5,'HandleVisibility', 'off');
hold on;
plot(avgKMO,avgKMH, '^', 'DisplayName', 'Knik', 'Color', 'k', 'MarkerFaceColor', [0.35,0.70,0.90],'MarkerSize',m);
plot(avgMMO, avgMMH, 's', 'DisplayName', 'Matanuska', 'Color', 'k', 'MarkerFaceColor', [0.95,0.90,0.25],'MarkerSize',m);
plot(avgLSMO, avgLSMH, 'd', 'DisplayName', 'Little Susitna', 'Color', 'k', 'MarkerFaceColor', [0,0.60,0.50],'MarkerSize',m);
plot(avgMCO, avgMCH, '>', 'DisplayName', 'Moose Creek', 'Color', 'k', 'MarkerFaceColor', [0.90,0.60,0],'MarkerSize',m);
plot(avgCTMO, avgCTMH, 'v', 'DisplayName', 'Castner', 'Color', 'k', 'MarkerFaceColor', [0,0.45,0.70],'MarkerSize',m);
plot(avgCWMO, avgCWMH, 'pentagram', 'DisplayName', 'Canwell', 'Color', 'k', 'MarkerFaceColor', [0.80,0.60,0.70],'MarkerSize',m);
plot(avgGMO, avgGMH, 'o', 'DisplayName', 'Gulkana', 'Color', 'k', 'MarkerFaceColor', [0.80,0.40,0],'MarkerSize',m);

xlabel('δ¹⁸O (‰)', 'FontSize', 30); 
ylabel('δ²H (‰)', 'FontSize', 30);  
%title('Stable isotope comparison');
ax = gca;
ax.FontSize = 15; 
legendHandle = legend('Location', 'northwest');
set(legendHandle, 'FontSize',18);
set(legendHandle, 'ItemTokenSize', [30,30]);
grid off;
hold off;
set(gcf, 'Position', [100, 100, 1000, 800]);
folderName = 'U:/GoA plots/NewPlots';
fileName = 'isotope_comparison_watershed_avg.svg';
fullFilePath = fullfile(folderName, fileName);
saveas(gcf, fullFilePath, 'svg');
%% Stable isotope comparison: watershed with averages
figure;
ms=12;
msS=16;
% Define the range for DO (δ¹⁸O) based on your x-axis range
DO_range = linspace(-24,-18, 100); 
% Compute DD (δ²H) using the meteoric water line equation
DD_meteoric = 8 * DO_range + 10;
% Plot the meteoric water line
plot(DO_range, DD_meteoric, 'k--', 'DisplayName', 'GMWL', 'LineWidth', 0.5);
hold on;
plot(avgKMO, avgKMH, 'o', 'DisplayName', 'Knik','MarkerEdgeColor', [0.35,0.70,0.90],'MarkerSize', ms, 'LineWidth', 2.5);
plot(avgKSO, avgKSH, 'kd', 'DisplayName', 'Knik Springs','MarkerFaceColor', [0.35,0.70,0.90],'MarkerSize', ms, 'LineWidth', 0.5);
plot(avgKTO, avgKTH, 'ks', 'DisplayName', 'Knik Tributaries','MarkerFaceColor', [0.35,0.70,0.90],'MarkerSize', msS, 'LineWidth', 0.5);
plot(avgMMO, avgMMH, 'o', 'DisplayName', 'Matanuska','MarkerEdgeColor', [0.95,0.90,0.25],'MarkerSize', ms, 'LineWidth', 2.5);
plot(avgMSO, avgMSH, 'kd', 'DisplayName', 'Matanuska Springs','MarkerFaceColor', [0.95,0.90,0.25],'MarkerSize', ms, 'LineWidth', 0.5);
plot(avgMTO, avgMTH, 'ks', 'DisplayName', 'Matanuska Tributaries','MarkerFaceColor', [0.95,0.90,0.25],'MarkerSize', msS, 'LineWidth', 0.5);
plot(avgMCO, avgMCH, 'o', 'DisplayName', 'Moose Creek','MarkerEdgeColor', [0.90,0.60,0],'MarkerSize', ms, 'LineWidth', 2.5);
plot(avgLSMO, avgLSMH,'o', 'DisplayName', 'Little Susitna','MarkerEdgeColor', [0,0.60,0.50],'MarkerSize', ms, 'LineWidth', 2.5);
plot(avgLSSO, avgLSSH,'kd', 'DisplayName', 'Little Susitna Springs','MarkerFaceColor', [0,0.60,0.50],'MarkerSize', ms, 'LineWidth', 0.5);
plot(avgLSTO, avgLSTH,'ks', 'DisplayName', 'Little Susitna Tributaries','MarkerFaceColor', [0,0.60,0.50],'MarkerSize', msS, 'LineWidth', 0.5);
plot(avgLSSubO, avgLSSubH,'kp', 'DisplayName', 'Little Susitna Subglacial','MarkerFaceColor', [0,0.60,0.50],'MarkerSize', msS, 'LineWidth', 0.5);
plot(avgLSPeriO, avgLSPeriH,'k^', 'DisplayName', 'Little Susitna Periglacial','MarkerFaceColor', [0,0.60,0.50],'MarkerSize', ms, 'LineWidth', 0.5);
plot(avgLSSupraO, avgLSSupraH,'kv', 'DisplayName', 'Little Susitna Supraglacial','MarkerFaceColor', [0,0.60,0.50],'MarkerSize', ms, 'LineWidth', 0.5);
plot(avgLSLO, avgLSLH,'k>', 'DisplayName', 'Little Susitna Glacial Lake','MarkerFaceColor', [0,0.60,0.50],'MarkerSize', ms, 'LineWidth', 0.5);
plot(avgGMO,avgGMH, 'o', 'DisplayName', 'Gulkana','MarkerEdgeColor', [0.80,0.40,0],'MarkerSize', ms, 'LineWidth', 2.5);
plot(mGulkanaSupra(:,9), mGulkanaSupra(:,10), 'kv', 'DisplayName', 'Gulkana Supraglacial','MarkerFaceColor', [0.80,0.40,0],'MarkerSize', ms, 'LineWidth', 0.5);
plot(avgGTO, avgGTH, 'ks', 'DisplayName', 'Gulkana Tributaries','MarkerFaceColor', [0.80,0.40,0],'MarkerSize', msS, 'LineWidth', 0.5);
plot(avgGSubO, avgGSubH, 'kp', 'DisplayName', 'Gulkana Subglacial','MarkerFaceColor', [0.80,0.40,0],'MarkerSize', msS, 'LineWidth', 0.5);
plot(avgCTMO, avgCTMH,'o', 'DisplayName', 'Castner','MarkerEdgeColor', [0,0.45,0.70],'MarkerSize', ms, 'LineWidth', 2.5);
plot(avgCTSubO, avgCTSubH,'kp', 'DisplayName', 'Castner Subglacial','MarkerFaceColor', [0,0.45,0.70],'MarkerSize', msS, 'LineWidth', 0.5);
plot(avgCWMO, avgCWMH, 'o', 'DisplayName', 'Canwell','MarkerEdgeColor', [0.80,0.60,0.70],'MarkerSize', ms, 'LineWidth', 2.5);
xlabel('δ¹⁸O (‰)', 'FontSize', 20); 
ylabel('δ²H (‰)', 'FontSize', 20);  
%title('Stable isotope comparison');
ax = gca;
ax.FontSize = 12; 
% legendHandle = legend('Location', 'southeast');
% set(legendHandle, 'FontSize',9, 'NumColumns', 2, 'ItemTokenSize', [8,8]);
grid off;
hold off;
set(gcf, 'Position', [100, 100, 600, 600]);
folderName = 'U:/GoA plots/NewPlots';
fileName = 'isotope_comparison_watershed_averages.jpg';
fullFilePath = fullfile(folderName, fileName);
saveas(gcf, fullFilePath, 'jpg');

%% Isotope script with local meteoric water lines
DO_range = linspace(-26, -17, 100); % Define the range for DO (δ¹⁸O) for regression lines
DD_meteoric = 8 * DO_range + 10; % Compute Global Meteoric Water Line (GMWL)

figure;
plot(DO_range, DD_meteoric, 'k--', 'DisplayName', 'GMWL', 'LineWidth', 0.5);
hold on;
ms = 7; % Marker sizes
msS = 10;

% Knik
Knik_combined = [mKnikMain; mKnikTrib];
Knik_combined = Knik_combined(~any(isnan(Knik_combined(:, [9, 10])), 2), :);
coeff_Knik = polyfit(Knik_combined(:, 9), Knik_combined(:, 10), 1);
plot(DO_range, polyval(coeff_Knik, DO_range), '-', 'Color', [0.35, 0.70, 0.90], ...
    'LineWidth', 1.5, 'DisplayName', 'Knik LMWL','HandleVisibility', 'off');

plot(mKnikMain(:,9), mKnikMain(:,10), 'o', 'DisplayName', 'Knik','MarkerEdgeColor', [0.35,0.70,0.90],'MarkerSize', ms, 'LineWidth', 2.5);
plot(mKnikSpring(:,9), mKnikSpring(:,10), 'kd', 'DisplayName', 'Knik Springs','MarkerFaceColor', [0.35,0.70,0.90],'MarkerSize', ms, 'LineWidth', 0.5);
plot(mKnikTrib(:,9), mKnikTrib(:,10), 'ks', 'DisplayName', 'Knik Tributaries','MarkerFaceColor', [0.35,0.70,0.90],'MarkerSize', msS, 'LineWidth', 0.5);

% Matanuska: Combine datasets
Matanuska_combined = [mMatMain; mMatTrib];
Matanuska_combined = Matanuska_combined(~any(isnan(Matanuska_combined(:, [9, 10])), 2), :);
coeff_Matanuska = polyfit(Matanuska_combined(:, 9), Matanuska_combined(:, 10), 1);
plot(DO_range, polyval(coeff_Matanuska, DO_range), '-', 'Color', [0.95, 0.90, 0.25], ...
    'LineWidth', 1.5, 'DisplayName', 'Matanuska LMWL','HandleVisibility', 'off');

plot(mMatMain(:,9), mMatMain(:,10), 'o', 'DisplayName', 'Matanuska','MarkerEdgeColor', [0.95,0.90,0.25],'MarkerSize', ms, 'LineWidth', 2.5);
plot(mMatSpring(:,9), mMatSpring(:,10), 'kd', 'DisplayName', 'Matanuska Springs','MarkerFaceColor', [0.95,0.90,0.25],'MarkerSize', ms, 'LineWidth', 0.5);
plot(mMatTrib(:,9), mMatTrib(:,10), 'ks', 'DisplayName', 'Matanuska Tributaries','MarkerFaceColor', [0.95,0.90,0.25],'MarkerSize', msS, 'LineWidth', 0.5);

% plot(mMC_2022(:,9), mMC_2022(:,10), 'o', 'DisplayName', 'Moose Creek','MarkerEdgeColor', [0.90,0.60,0],'MarkerSize', ms, 'LineWidth', 2.5);
% plot(mMC_2023(:,9), mMC_2023(:,10), 'o', 'DisplayName', 'Moose Creek','MarkerEdgeColor', [0.90,0.60,0],'MarkerSize', ms, 'LineWidth', 2.5,'HandleVisibility', 'off');
% 
% % Moose Creek: Single dataset
% MC_combined = [mMC_2022, mMC_2023];
% MC_combined = MC_combined(~any(isnan(MC_combined(:, [9, 10])), 2), :);
% coeff_MC = polyfit(MC_combined(:, 9), MC_combined(:, 10), 1);
% plot(DO_range, polyval(coeff_MC, DO_range), '-', 'Color', [0.90, 0.60, 0], ...
%     'LineWidth', 1.5, 'DisplayName', 'Moose Creek LMWL','HandleVisibility', 'off');

% Little Susitna: Combine datasets
LS_combined = [mLSMain; mLSTrib; mLSSub; mLSPeri; mLSSupra; mLSLake];
LS_combined = LS_combined(~any(isnan(LS_combined(:, [9, 10])), 2), :);
coeff_LS = polyfit(LS_combined(:, 9), LS_combined(:, 10), 1);
plot(DO_range, polyval(coeff_LS, DO_range), '-', 'Color', [0, 0.60, 0.50], ...
    'LineWidth', 1.5, 'DisplayName', 'Little Susitna LMWL','HandleVisibility', 'off');

plot(mLSMain(:,9), mLSMain(:,10),'o', 'DisplayName', 'Little Susitna','MarkerEdgeColor', [0,0.60,0.50],'MarkerSize', ms, 'LineWidth', 2.5);
plot(mLSSpring(:,9), mLSSpring(:,10),'kd', 'DisplayName', 'Little Susitna Springs','MarkerFaceColor', [0,0.60,0.50],'MarkerSize', ms, 'LineWidth', 0.5);
plot(mLSTrib(:,9), mLSTrib(:,10),'ks', 'DisplayName', 'Little Susitna Tributaries','MarkerFaceColor', [0,0.60,0.50],'MarkerSize', msS, 'LineWidth', 0.5);
plot(mLSSub(:,9), mLSSub(:,10),'kp', 'DisplayName', 'Little Susitna Subglacial','MarkerFaceColor', [0,0.60,0.50],'MarkerSize', msS, 'LineWidth', 0.5);
plot(mLSPeri(:,9), mLSPeri(:,10),'k^', 'DisplayName', 'Little Susitna Periglacial','MarkerFaceColor', [0,0.60,0.50],'MarkerSize', ms, 'LineWidth', 0.5);
plot(mLSSupra(:,9), mLSSupra(:,10),'kv', 'DisplayName', 'Little Susitna Supraglacial','MarkerFaceColor', [0,0.60,0.50],'MarkerSize', ms, 'LineWidth', 0.5);
plot(mLSLake(:,9), mLSLake(:,10),'k>', 'DisplayName', 'Little Susitna Glacial Lake','MarkerFaceColor', [0,0.60,0.50],'MarkerSize', ms, 'LineWidth', 0.5);

% Gulkana: Combine datasets
Gulkana_combined = [mGulkanaMain; mGulkanaTrib; mGulkanaSub; mGulkanaSupra];
Gulkana_combined = Gulkana_combined(~any(isnan(Gulkana_combined(:, [9, 10])), 2), :);
coeff_Gulkana = polyfit(Gulkana_combined(:, 9), Gulkana_combined(:, 10), 1);
plot(DO_range, polyval(coeff_Gulkana, DO_range), '-', 'Color', [0.80, 0.40, 0], ...
    'LineWidth', 1.5, 'DisplayName', 'Gulkana LMWL','HandleVisibility', 'off');

plot(mGulkanaMain(:,9), mGulkanaMain(:,10), 'o', 'DisplayName', 'Gulkana','MarkerEdgeColor', [0.80,0.40,0],'MarkerSize', ms, 'LineWidth', 2.5);
plot(mGulkanaSupra(:,9), mGulkanaSupra(:,10), 'kv', 'DisplayName', 'Gulkana Supraglacial','MarkerFaceColor', [0.80,0.40,0],'MarkerSize', ms, 'LineWidth', 0.5);
plot(mGulkanaTrib(:,9), mGulkanaTrib(:,10), 'ks', 'DisplayName', 'Gulkana Tributaries','MarkerFaceColor', [0.80,0.40,0],'MarkerSize', msS, 'LineWidth', 0.5);
plot(mGulkanaSub(:,9), mGulkanaSub(:,10), 'kp', 'DisplayName', 'Gulkana Subglacial','MarkerFaceColor', [0.80,0.40,0],'MarkerSize', msS, 'LineWidth', 0.5);

% Castner: Combine datasets
Castner_combined = [mCastnerMain; mCastnerSub];
Castner_combined = Castner_combined(~any(isnan(Castner_combined(:, [9, 10])), 2), :);
coeff_Castner = polyfit(Castner_combined(:, 9), Castner_combined(:, 10), 1);
plot(DO_range, polyval(coeff_Castner, DO_range), '-', 'Color', [0, 0.45, 0.70], ...
    'LineWidth', 1.5, 'DisplayName', 'Castner LMWL','HandleVisibility', 'off');

plot(mCastnerMain(:,9), mCastnerMain(:,10),'o', 'DisplayName', 'Castner','MarkerEdgeColor', [0,0.45,0.70],'MarkerSize', ms, 'LineWidth', 2.5);
plot(mCastnerSub(:,9), mCastnerSub(:,10),'kp', 'DisplayName', 'Castner Subglacial','MarkerFaceColor', [0,0.45,0.70],'MarkerSize', msS, 'LineWidth', 0.5);

% Canwell: Single dataset
coeff_Canwell = polyfit(mCanwellData(:, 9), mCanwellData(:, 10), 1);
plot(DO_range, polyval(coeff_Canwell, DO_range), '-', 'Color', [0.80, 0.60, 0.70], ...
    'LineWidth', 1.5, 'DisplayName', 'Canwell LMWL','HandleVisibility', 'off');

plot(mCanwellData(:,9), mCanwellData(:,10), 'o', 'DisplayName', 'Canwell','MarkerEdgeColor', [0.80,0.60,0.70],'MarkerSize', ms, 'LineWidth', 2.5);

% Add labels, legend, and save
xlabel('δ¹⁸O (‰)', 'FontSize', 20); 
ylabel('δ²H (‰)', 'FontSize', 20);  
ax = gca;
ax.FontSize = 15; 
legendHandle = legend('Location', 'northwest');
set(legendHandle, 'FontSize', 11, 'NumColumns', 2);
set(legendHandle, 'ItemTokenSize', [20, 20]);
grid off;
hold off;
set(gcf, 'Position', [100, 100, 1400, 1000]);

% Save the plot
folderName = 'U:/GoA plots/NewPlots';
fileName = 'isotope_LMWLs.svg';
fullFilePath = fullfile(folderName, fileName);
saveas(gcf, fullFilePath, 'svg');

%% LMWL plot
DO_range = linspace(-26, -17, 100); % Define the range for DO (δ¹⁸O) for regression lines
DD_meteoric = 8 * DO_range + 10; % Compute Global Meteoric Water Line (GMWL)

figure;
plot(DO_range, DD_meteoric, 'k--', 'DisplayName', 'GMWL', 'LineWidth', 0.5);
hold on;
ms = 7; % Marker sizes
msS = 10;
% Knik
% Knik_combined = [mKnikMain; mKnikTrib];
% Knik_combined = Knik_combined(~any(isnan(Knik_combined(:, [9, 10])), 2), :);
% coeff_Knik = polyfit(Knik_combined(:, 9), Knik_combined(:, 10), 1);
% plot(DO_range, polyval(coeff_Knik, DO_range), '-', 'Color', [0.35, 0.70, 0.90], ...
%     'LineWidth', 1.5, 'DisplayName', 'Knik LMWL');
% Matanuska: Combine datasets
% Matanuska_combined = [mMatMain; mMatTrib];
% Matanuska_combined = Matanuska_combined(~any(isnan(Matanuska_combined(:, [9, 10])), 2), :);
% coeff_Matanuska = polyfit(Matanuska_combined(:, 9), Matanuska_combined(:, 10), 1);
% plot(DO_range, polyval(coeff_Matanuska, DO_range), '-', 'Color', [0.95, 0.90, 0.25], ...
%     'LineWidth', 1.5, 'DisplayName', 'Matanuska LMWL');
% % Moose Creek: Single dataset
% MC_combined = [mMC_2022, mMC_2023];
% MC_combined = MC_combined(~any(isnan(MC_combined(:, [9, 10])), 2), :);
% coeff_MC = polyfit(MC_combined(:, 9), MC_combined(:, 10), 1);
% plot(DO_range, polyval(coeff_MC, DO_range), '-', 'Color', [0.90, 0.60, 0], ...
%     'LineWidth', 2.5, 'DisplayName', 'Moose Creek LMWL');
% Little Susitna: Combine datasets
% LS_combined = [mLSMain; mLSTrib; mLSSub; mLSPeri; mLSSupra; mLSLake];
% LS_combined = LS_combined(~any(isnan(LS_combined(:, [9, 10])), 2), :);
% coeff_LS = polyfit(LS_combined(:, 9), LS_combined(:, 10), 1);
% plot(DO_range, polyval(coeff_LS, DO_range), '-', 'Color', [0, 0.60, 0.50], ...
%     'LineWidth', 1.5, 'DisplayName', 'Little Susitna LMWL');
% Gulkana: Combine datasets
% Gulkana_combined = [mGulkanaMain; mGulkanaTrib; mGulkanaSub; mGulkanaSupra];
% Gulkana_combined = Gulkana_combined(~any(isnan(Gulkana_combined(:, [9, 10])), 2), :);
% coeff_Gulkana = polyfit(Gulkana_combined(:, 9), Gulkana_combined(:, 10), 1);
% plot(DO_range, polyval(coeff_Gulkana, DO_range), '-', 'Color', [0.80, 0.40, 0], ...
%     'LineWidth', 1.5, 'DisplayName', 'Gulkana LMWL');
% Castner: Combine datasets
% Castner_combined = [mCastnerMain; mCastnerSub];
% Castner_combined = Castner_combined(~any(isnan(Castner_combined(:, [9, 10])), 2), :);
% coeff_Castner = polyfit(Castner_combined(:, 9), Castner_combined(:, 10), 1);
% plot(DO_range, polyval(coeff_Castner, DO_range), '-', 'Color', [0, 0.45, 0.70], ...
%     'LineWidth', 1.5, 'DisplayName', 'Castner LMWL');
% Canwell: Single dataset
% coeff_Canwell = polyfit(mCanwellData(:, 9), mCanwellData(:, 10), 1);
% plot(DO_range, polyval(coeff_Canwell, DO_range), '-', 'Color', [0.80, 0.60, 0.70], ...
%     'LineWidth', 1.5, 'DisplayName', 'Canwell LMWL');
% Plot average lines
plot(DO_range, polyval(avg_coeff_Group1, DO_range), ':', 'Color', 'r', ...
    'LineWidth', 2, 'DisplayName', 'Southcentral Avg LMWL');
plot(DO_range, polyval(avg_coeff_Group2, DO_range), ':', 'Color', 'g', ...
    'LineWidth', 2, 'DisplayName', 'Delta Range Avg LMWL');
plot(DO_range, polyval(avg_coeff_All, DO_range), '-.', 'Color', [0.2, 0.2, 0.2], ...
    'LineWidth', 2.5, 'DisplayName', 'Overall Avg LMWL');
% Alaska MWL
amwl_slope = 7.95;
amwl_intercept = 6.37;
DD_amwl = amwl_slope * DO_range + amwl_intercept;
plot(DO_range, DD_amwl, '--', 'Color','b', 'LineWidth', 2, 'DisplayName', 'Alaska MWL (Lachniet et al., 2016)');

xlabel('δ¹⁸O (‰)', 'FontSize', 20); 
ylabel('δ²H (‰)', 'FontSize', 20);  
ax = gca;
ax.FontSize = 15; 
legendHandle = legend('Location', 'northwest');
set(legendHandle, 'FontSize', 12);
set(legendHandle, 'ItemTokenSize', [20, 20]);
grid off;
hold off;
set(gcf, 'Position', [100, 100, 1000, 1000]);
folderName = 'U:/GoA plots/NewPlots';
fileName = 'LMWLs.jpg';
fullFilePath = fullfile(folderName, fileName);
saveas(gcf, fullFilePath, 'jpg');
%% Pring equations
% Initialize arrays to store slope and intercept values
coeff_Group1 = []; % For Knik, Matanuska, Little Susitna, Moose Creek
coeff_Group2 = []; % For Gulkana, Castner, Canwell

% Define a function to print regression equation
printEquation = @(coeff, name) fprintf('%s LMWL: δ²H = %.2f * δ¹⁸O + %.2f\n', name, coeff(1), coeff(2));

% Knik
coeff_Knik = polyfit(Knik_combined(:, 9), Knik_combined(:, 10), 1);
coeff_Group1 = [coeff_Group1; coeff_Knik];
printEquation(coeff_Knik, 'Knik');

% Matanuska
coeff_Matanuska = polyfit(Matanuska_combined(:, 9), Matanuska_combined(:, 10), 1);
coeff_Group1 = [coeff_Group1; coeff_Matanuska];
printEquation(coeff_Matanuska, 'Matanuska');

% % Moose Creek
% coeff_MC = polyfit(MC_combined(:, 9), MC_combined(:, 10), 1);
% coeff_Group1 = [coeff_Group1; coeff_MC];
% printEquation(coeff_MC, 'Moose Creek');

% Little Susitna
coeff_LS = polyfit(LS_combined(:, 9), LS_combined(:, 10), 1);
coeff_Group1 = [coeff_Group1; coeff_LS];
printEquation(coeff_LS, 'Little Susitna');

% Gulkana
coeff_Gulkana = polyfit(Gulkana_combined(:, 9), Gulkana_combined(:, 10), 1);
coeff_Group2 = [coeff_Group2; coeff_Gulkana];
printEquation(coeff_Gulkana, 'Gulkana');

% Castner
coeff_Castner = polyfit(Castner_combined(:, 9), Castner_combined(:, 10), 1);
coeff_Group2 = [coeff_Group2; coeff_Castner];
printEquation(coeff_Castner, 'Castner');

% Canwell
coeff_Canwell = polyfit(mCanwellData(:, 9), mCanwellData(:, 10), 1);
coeff_Group2 = [coeff_Group2; coeff_Canwell];
printEquation(coeff_Canwell, 'Canwell');

% Compute average lines
avg_coeff_Group1 = mean(coeff_Group1, 1); % Average slope and intercept for Group 1
avg_coeff_Group2 = mean(coeff_Group2, 1); % Average slope and intercept for Group 2
coeff_All = [coeff_Group1; coeff_Group2];
avg_coeff_All = mean(coeff_All, 1);

% Print average equations
fprintf('Group 1 Average LMWL: δ²H = %.2f * δ¹⁸O + %.2f\n', avg_coeff_Group1(1), avg_coeff_Group1(2));
fprintf('Group 2 Average LMWL: δ²H = %.2f * δ¹⁸O + %.2f\n', avg_coeff_Group2(1), avg_coeff_Group2(2));
fprintf('Overall Average LMWL: δ²H = %.2f * δ¹⁸O + %.2f\n', avg_coeff_All(1), avg_coeff_All(2));
fprintf('Alaska MWL: δ²H = %.2f * δ¹⁸O + %.2f\n', amwl_slope, amwl_intercept);

% Plot average lines
% plot(DO_range, polyval(avg_coeff_Group1, DO_range), '--', 'Color', [0.4, 0.4, 0.4], ...
    % 'LineWidth', 2, 'DisplayName', 'Southcentral Avg LMWL');
% plot(DO_range, polyval(avg_coeff_Group2, DO_range), '--', 'Color', [0.6, 0.2, 0.2], ...
    % 'LineWidth', 2, 'DisplayName', 'Alaska Range Avg LMWL');
% plot(DO_range, polyval(avg_coeff_All, DO_range), '--', 'Color', [0.2, 0.2, 0.2], ...
    % 'LineWidth', 2.5, 'DisplayName', 'Overall Avg LMWL');
%% Conc v time with discharge: USGS - done
columnNames = vColumnLabels;
folderName = 'U:/GoA plots/NewPlots';

dPosition = 10;  % Static position for the column

mKR3_2022 = mKR3Data(mKR3Year == 2022, :);
mKR3_Date22 = mKR3Dates(mKR3Year == 2022);
mKR3_2022(2, :) = [];
mKR3_Date22(2) = [];
% sTitle = ['Element ', sInput];    
figure;
set(gcf, 'Position', [600 100 400 400]);

% Plot data for KR3
y_KR3_2022 = mKR3_2022(:, dPosition);
y_KR3_2022 = fillmissing(y_KR3_2022, 'linear');
plot(mKR3_Date22, y_KR3_2022, 'o-', 'DisplayName', 'K3', 'Color', [0.35, 0.70, 0.90], ...
    'MarkerFaceColor', [0.35, 0.70, 0.90], 'MarkerEdgeColor', 'k');
hold on;

% Plot data for MR4
y_MR4_2022 = mMR4_2022(:, dPosition);
y_MR4_2022 = fillmissing(y_MR4_2022, 'linear');
plot(mMR4_Date22, y_MR4_2022, 'o-', 'DisplayName', 'M4', 'Color', [0.95, 0.90, 0.25], ...
    'MarkerFaceColor', [0.95, 0.90, 0.25], 'MarkerEdgeColor', 'k');

% Plot data for MC
y_MC_2022 = mMC_2022(:, dPosition);
y_MC_2022 = fillmissing(y_MC_2022, 'linear');
plot(mMC_Date22, y_MC_2022, 'o-', 'DisplayName', 'MC', 'Color', [0.90, 0.60, 0], ...
    'MarkerFaceColor', [0.90, 0.60, 0], 'MarkerEdgeColor', 'k');

% Plot data for LS2
y_LS2_2022 = mLS2_2022(:, dPosition);
y_LS2_2022 = fillmissing(y_LS2_2022, 'linear');
plot(mLS2_Date22, y_LS2_2022, 'o-', 'DisplayName', 'LS2', 'Color', [0, 0.60, 0.50], ...
    'MarkerFaceColor', [0, 0.60, 0.50], 'MarkerEdgeColor', 'k');

% Labels and Title
ylabel('δ²H (‰)', 'FontSize', 14);
title('δ²H', 'FontSize', 18);

% Adjust Legend
legendObj = legend('Location', 'southeast');
legendObj.FontSize = 12;  
ax = gca;
ax.XAxis.FontSize = 12;
ax.YColor = 'k';
ax.YAxis.FontSize = 12;
grid on;
% set(gca, 'YScale', 'log');
hold off;

% Save the figure
sInputChar = char(sInput);
formatFileName = "Hydrogen_iso.jpg";
fileName = sprintf(formatFileName, sInputChar);
fullFilePath = fullfile(folderName, fileName);
saveas(gcf, fullFilePath, 'jpg');


