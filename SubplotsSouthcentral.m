%% Plots
clear all, close all, clc %use AlaskaMasterdatabase2022_2023_spatial
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
%% Plot LS data
% Extract month from dates
columnNames = vColumnLabels;
folderName = 'U:/GoA plots/NewPlots';
years = year(mLSDates); % Extract years from dates

% Define the months for which to create lines
targetMonths = [5, 6, 7, 8, 9];
colors = {[0.8, 0.6, 0.7], [0.9, 0.6, 0], [0.95, 0.9, 0.25], [0, 0.6, 0.5], [0.35, 0.7, 0.9]};
monthLabels = {'May', 'June', 'July', 'August', 'September'};

% Define the dPosition values to iterate over
dPositions = [42, 29, 11, 60]; % e.g., As, Ca, HCO₃, Yb
customLabels = {'As', 'Ca', 'HCO_3', 'Yb'};

% Initialize array to store y-axis limits for each dPosition
yLimits = zeros(numel(dPositions), 2); % Each row stores [yMin, yMax] for a pair of plots

% Create a figure with 8 subplots (2 rows x 4 columns)
figure;
set(gcf, 'Position', [100, 100, 1400, 600]);

for yearIndex = 1:2 % 1 for 2022, 2 for 2023
    if yearIndex == 1
        currentYear = 2022;
    else
        currentYear = 2023;
    end
    
    for dpIndex = 1:numel(dPositions)
        dPosition = dPositions(dpIndex);
        sInput = columnNames{dPosition}; % Get the label for the current dPosition
        subplotIndex = (yearIndex - 1) * 4 + dpIndex; % Calculate subplot index for 2x4 layout
        subplot(2, 4, subplotIndex); % Create subplot
        
        % Filter data for the current year and dPosition
        currentData = mLSData(years == currentYear, :);
        currentSites = mLSSites(years == currentYear);
        months = month(mLSDates(years == currentYear));
        
        % Get unique and sorted site names
        uniqueSites = unique(currentSites);
        sortedSites = sort(uniqueSites);
        
        % Plot data for each month
        hold on;
        for i = 1:numel(targetMonths)
            xCoords = [];
            yCoords = [];
            for j = 1:numel(sortedSites)
                currentSite = sortedSites{j};
                monthData = currentData(strcmp(currentSites, currentSite) & months == targetMonths(i), dPosition);
                nonNaNIndices = ~isnan(monthData);
                xCoords = [xCoords, repmat(j, 1, sum(nonNaNIndices))];
                yCoords = [yCoords, monthData(nonNaNIndices)];
            end
            if ~isempty(xCoords)
                plot(xCoords, yCoords, 'o-', 'Color', colors{i}, 'MarkerFaceColor', colors{i}, 'MarkerEdgeColor', 'k', 'DisplayName', monthLabels{i});
            end
        end
        
        % Customize plot
        xticks(1:numel(sortedSites));
        xticklabels(sortedSites);
        title([sInput ' ' num2str(currentYear)], 'Interpreter', 'tex');
        if dpIndex == 1
            ylabel('Concentration (mg/L)');
        end
        
        % Add legend only in the first subplot (2022, As)
        if dpIndex == 1 && yearIndex == 1
            legendHandle = legend('Location', 'northeast');
            set(legendHandle, 'FontSize', 5.5);
            set(legendHandle, 'ItemTokenSize', [5, 5]);
        end
        
        % Store y-axis limits for comparison
        currentYLimits = ylim;
        if yearIndex == 1
            yLimits(dpIndex, :) = currentYLimits; % Store 2022 limits
        else
            yLimits(dpIndex, :) = [min(yLimits(dpIndex, 1), currentYLimits(1)), ...
                                   max(yLimits(dpIndex, 2), currentYLimits(2))]; % Update with 2023 limits
        end
    end
end

% Apply consistent y-axis limits across each column
for dpIndex = 1:numel(dPositions)
    for yearIndex = 1:2
        subplotIndex = (yearIndex - 1) * 4 + dpIndex;
        subplot(2, 4, subplotIndex);
        ylim(yLimits(dpIndex, :)); % Set matching y-limits
    end
end

% Save the figure
saveas(gcf, fullfile(folderName, 'HorizontalSubplotsLS.svg'), 'svg');


%% Plot Mat data
% Extract month from dates
columnNames = vColumnLabels;
folderName = 'U:/GoA plots/NewPlots';
years = year(mMatDates); % Extract years from dates

% Define the months for which to create lines
targetMonths = [5, 6, 7, 8, 9];
colors = {[0.8, 0.6, 0.7], [0.9, 0.6, 0], [0.95, 0.9, 0.25], [0, 0.6, 0.5], [0.35, 0.7, 0.9]};
monthLabels = {'May', 'June', 'July', 'August', 'September'};

% Define the dPosition values to iterate over
dPositions = [42, 29, 11, 60]; % e.g., As, Ca, HCO₃, Yb
customLabels = {'As', 'Ca', 'HCO_3', 'Yb'};

% Initialize array to store y-axis limits for each dPosition
yLimits = zeros(numel(dPositions), 2); % Each row stores [yMin, yMax] for each column

% Create a figure with 8 subplots (2 rows x 4 columns)
figure;
set(gcf, 'Position', [100, 100, 1400, 600]); 

for yearIndex = 1:2 % 1 for 2022, 2 for 2023
    if yearIndex == 1
        currentYear = 2022;
    else
        currentYear = 2023;
    end
    
    for dpIndex = 1:numel(dPositions)
        dPosition = dPositions(dpIndex);
        sInput = columnNames{dPosition}; % Get the label for the current dPosition
        subplotIndex = (yearIndex - 1) * 4 + dpIndex; % Calculate subplot index for 2x4 layout
        subplot(2, 4, subplotIndex); % Create subplot
        
        % Filter data for the current year and dPosition
        currentData = mMatData(years == currentYear, :);
        currentSites = mMatSites(years == currentYear);
        months = month(mMatDates(years == currentYear));
        
        % Get unique and sorted site names
        uniqueSites = unique(currentSites);
        sortedSites = sort(uniqueSites);
        
        % Plot data for each month
        hold on;
        for i = 1:numel(targetMonths)
            xCoords = [];
            yCoords = [];
            for j = 1:numel(sortedSites)
                currentSite = sortedSites{j};
                monthData = currentData(strcmp(currentSites, currentSite) & months == targetMonths(i), dPosition);
                nonNaNIndices = ~isnan(monthData);
                xCoords = [xCoords, repmat(j, 1, sum(nonNaNIndices))];
                yCoords = [yCoords, monthData(nonNaNIndices)];
            end
            if ~isempty(xCoords)
                plot(xCoords, yCoords, 'o-', 'Color', colors{i}, 'MarkerFaceColor', colors{i}, 'MarkerEdgeColor', 'k', 'DisplayName', monthLabels{i});
            end
        end
        
        % Customize plot
        xticks(1:numel(sortedSites));
        xticklabels(sortedSites);
        title([sInput ' ' num2str(currentYear)], 'Interpreter', 'tex');
        if dpIndex == 1
            ylabel('Concentration (mg/L)');
        end
        
        % Store y-axis limits for comparison
        currentYLimits = ylim;
        if yearIndex == 1
            yLimits(dpIndex, :) = currentYLimits; % Store 2022 limits
        else
            yLimits(dpIndex, :) = [min(yLimits(dpIndex, 1), currentYLimits(1)), ...
                                   max(yLimits(dpIndex, 2), currentYLimits(2))]; % Update with 2023 limits
        end
        
        % Add legend only in the first subplot (2022, As)
        if dpIndex == 1 && yearIndex == 1
            legendHandle = legend('Location', 'northeast');
            set(legendHandle, 'FontSize', 5.5);
            set(legendHandle, 'ItemTokenSize', [5, 5]);
        end
    end
end

% Apply consistent y-axis limits across each column
for dpIndex = 1:numel(dPositions)
    for yearIndex = 1:2
        subplotIndex = (yearIndex - 1) * 4 + dpIndex;
        subplot(2, 4, subplotIndex);
        ylim(yLimits(dpIndex, :)); % Set matching y-limits
    end
end

% Save the figure
saveas(gcf, fullfile(folderName, 'HorizontalSubplotsMat.svg'), 'svg');



%% Plot Knik
% Extract month from dates
columnNames = vColumnLabels;
folderName = 'U:/GoA plots/NewPlots';
years = year(mKnikDates); % Extract years from dates

% Define the months for which to create lines
targetMonths = [5, 6, 7, 8, 9];
colors = {[0.8, 0.6, 0.7], [0.9, 0.6, 0], [0.95, 0.9, 0.25], [0, 0.6, 0.5], [0.35, 0.7, 0.9]};
monthLabels = {'May', 'June', 'July', 'August', 'September'};

% Define the dPosition values to iterate over
dPositions = [42, 29, 11, 60]; % e.g., As, Ca, HCO₃, Yb

% Initialize array to store y-axis limits for each dPosition
yLimits = zeros(numel(dPositions), 2); % Each row stores [yMin, yMax] for each column

% Create a figure with 8 subplots (2 rows x 4 columns)
figure;
set(gcf, 'Position', [100, 100, 1400, 600]); 

for yearIndex = 1:2 % 1 for 2022, 2 for 2023
    if yearIndex == 1
        currentYear = 2022;
    else
        currentYear = 2023;
    end
    
    for dpIndex = 1:numel(dPositions)
        dPosition = dPositions(dpIndex);
        sInput = columnNames{dPosition}; % Get the label for the current dPosition
        subplotIndex = (yearIndex - 1) * 4 + dpIndex; % Calculate subplot index for 2x4 layout
        subplot(2, 4, subplotIndex); % Create subplot
        
        % Filter data for the current year and dPosition
        currentData = mKnikData(years == currentYear, :);
        currentSites = mKnikSites(years == currentYear);
        months = month(mKnikDates(years == currentYear));
        
        % Get unique and sorted site names
        uniqueSites = unique(currentSites);
        sortedSites = sort(uniqueSites);
        
        % Plot data for each month
        hold on;
        for i = 1:numel(targetMonths)
            xCoords = [];
            yCoords = [];
            for j = 1:numel(sortedSites)
                currentSite = sortedSites{j};
                monthData = currentData(strcmp(currentSites, currentSite) & months == targetMonths(i), dPosition);
                nonNaNIndices = ~isnan(monthData);
                xCoords = [xCoords, repmat(j, 1, sum(nonNaNIndices))];
                yCoords = [yCoords, monthData(nonNaNIndices)];
            end
            if ~isempty(xCoords)
                plot(xCoords, yCoords, 'o-', 'Color', colors{i}, 'MarkerFaceColor', colors{i}, 'MarkerEdgeColor', 'k', 'DisplayName', monthLabels{i});
            end
        end
        
        % Customize plot
        xticks(1:numel(sortedSites));
        xticklabels(sortedSites);
        title([sInput ' ' num2str(currentYear)], 'Interpreter', 'tex');
        if dpIndex == 1
            ylabel('Concentration (mg/L)');
        end
        
        % Store y-axis limits for comparison
        currentYLimits = ylim;
        if yearIndex == 1
            yLimits(dpIndex, :) = currentYLimits; % Store 2022 limits
        else
            yLimits(dpIndex, :) = [min(yLimits(dpIndex, 1), currentYLimits(1)), ...
                                   max(yLimits(dpIndex, 2), currentYLimits(2))]; % Update with 2023 limits
        end
        
        % Add legend only in the first subplot (2022, As)
        if dpIndex == 1 && yearIndex == 1
            legendHandle = legend('Location', 'northeast');
            set(legendHandle, 'FontSize', 5.5);
            set(legendHandle, 'ItemTokenSize', [5, 5]);
        end
    end
end

% Apply consistent y-axis limits across each column
for dpIndex = 1:numel(dPositions)
    for yearIndex = 1:2
        subplotIndex = (yearIndex - 1) * 4 + dpIndex;
        subplot(2, 4, subplotIndex);
        ylim(yLimits(dpIndex, :)); % Set matching y-limits
    end
end

% Save the figure
saveas(gcf, fullfile(folderName, 'HorizontalSubplotsKnik.svg'), 'svg');
