%% Spatial Plots
clear all, close all, clc
[sFile, sPath] = uigetfile('*.xlsx', 'Select Database File');
sFullPath = fullfile(sPath, sFile);

tData = readtable(sFullPath);
%% 
tData = removevars(tData, 'ICP_MS');
vColumnLabels = tData.Properties.VariableNames(20:end);
vColumnLabelsArray = cellstr(vColumnLabels);
mFullData = table2array(tData(1:345, 20:end)); %first value is the number of samples (each as a row); second value is the numeric values that will be plotted, must be in an order that can be iterated through (DONT INCLUDE STRINGS)

vSampleDates = datetime(convertStringsToChars(string(table2cell(tData(:, 13))))); % just date not datetime
vSampleLocations = string(table2cell(tData(:, 2))); %labels of each sample (KR4, MR3, etc)
vWatershed = string(table2cell(tData(:, 4))); % column that could be used to discriminate between watersheds each location is contained in
vEndmember = string(table2cell(tData(:, 6))); % column that defines types of endmembers 
%% 
% Knik River 1 elemental data and sample dates
mKR1Data = mFullData(vSampleLocations == "KR1", :);
mKR1Dates = vSampleDates(vSampleLocations == "KR1");

% Knik River 2 elemental data and sample dates
mKR2Data = mFullData(vSampleLocations == "KR2", :);
mKR2Dates = vSampleDates(vSampleLocations == "KR2");

% Knik River 3 elemental data and sample dates
mKR3Data = mFullData(vSampleLocations == "KR3", :);
mKR3Dates = vSampleDates(vSampleLocations == "KR3");

% Knik River 4 elemental data and sample dates
mKR4Data = mFullData(vSampleLocations == "KR4", :);
mKR4Dates = vSampleDates(vSampleLocations == "KR4");

% Matanuska River 0.5 elemental data and sample dates
mMRHalfData = mFullData(vSampleLocations == "MR0.5", :);
mMRHalfDates = vSampleDates(vSampleLocations == "MR0.5");

% Matanuska River 1 elemental data and sample dates
mMR1Data = mFullData(vSampleLocations == "MR1", :);
mMR1Dates = vSampleDates(vSampleLocations == "MR1");

% Matanuska River 2 elemental data and sample dates
mMR2Data = mFullData(vSampleLocations == "MR2", :);
mMR2Dates = vSampleDates(vSampleLocations == "MR2");

% Matanuska River 3 elemental data and sample dates
mMR3Data = mFullData(vSampleLocations == "MR3", :);
mMR3Dates = vSampleDates(vSampleLocations == "MR3");

% Matanuska River 4 elemental data and sample dates
mMR4Data = mFullData(vSampleLocations == "MR4", :);
mMR4Dates = vSampleDates(vSampleLocations == "MR4");

% Matanuska River 5 elemental data and sample dates
mMR5Data = mFullData(vSampleLocations == "MR5", :);
mMR5Dates = vSampleDates(vSampleLocations == "MR5");

% Little Susitna River 1 elemental data and samples dates
mLS1Data = mFullData(vSampleLocations == "LS1", :);
mLS1Dates = vSampleDates(vSampleLocations == "LS1");

% Little Susitna River 1.5 elemental data and samples dates
mLS15Data = mFullData(vSampleLocations == "LS1.5", :);
mLS15Dates = vSampleDates(vSampleLocations == "LS1.5");

% Little Susitna River 2 elemental data and samples dates
mLS2Data = mFullData(vSampleLocations == "LS2", :);
mLS2Dates = vSampleDates(vSampleLocations == "LS2");

% Little Susitna River 3 elemental data and samples dates
mLS3Data = mFullData(vSampleLocations == "LS3", :);
mLS3Dates = vSampleDates(vSampleLocations == "LS3");

% Little Susitna River 4 elemental data and samples dates
mLS4Data = mFullData(vSampleLocations == "LS4", :);
mLS4Dates = vSampleDates(vSampleLocations == "LS4");

% Moose Creek elemental data and samples dates
mMCData = mFullData(vSampleLocations == "MC", :);
mMCDates = vSampleDates(vSampleLocations == "MC");

% Gulkana 1 elemental data and samples dates
mG1Data = mFullData(vSampleLocations == "G1", :);
mG1Dates = vSampleDates(vSampleLocations == "G1");

% Gulkana 2 elemental data and samples dates
mG2Data = mFullData(vSampleLocations == "G2", :);
mG2Dates = vSampleDates(vSampleLocations == "G2");

% Gulkana 3 elemental data and samples dates
mG3Data = mFullData(vSampleLocations == "3", :);
mG3Dates = vSampleDates(vSampleLocations == "3");

% Castner 1 elemental data and samples dates
mCT1Data = mFullData(vSampleLocations == "CT1", :);
mCT1Dates = vSampleDates(vSampleLocations == "CT1");

% Castner 2 elemental data and samples dates
mCT2Data = mFullData(vSampleLocations == "CT2", :);
mCT2Dates = vSampleDates(vSampleLocations == "CT2");

% Canwell 1 elemental data and samples dates
mCW1Data = mFullData(vSampleLocations == "CW1", :);
mCW1Dates = vSampleDates(vSampleLocations == "CW1");

% Canwell 2 elemental data and samples dates
mCW2Data = mFullData(vSampleLocations == "CW2", :);
mCW2Dates = vSampleDates(vSampleLocations == "CW2");

%% Define watersheds
mMatData = mFullData(vWatershed == "Matanuska", :);
mMatDates = vSampleDates(vWatershed == "Matanuska");
mMatSites = vSampleLocations(vWatershed == "Matanuska");

mKnikData = mFullData(vWatershed == "Knik", :);
mKnikDates = vSampleDates(vWatershed == "Knik");
mKnikSites = vSampleLocations(vWatershed == "Knik");

mLSData = mFullData(vWatershed == "Little Susitna", :);
mLSDates = vSampleDates(vWatershed == "Little Susitna");
mLSSites = vSampleLocations(vWatershed == "Little Susitna");

mGulkanaData = mFullData(vWatershed == "Gulkana", :);
mGulkanaDates = vSampleDates(vWatershed == "Gulkana");
mGulkanaSites = vSampleLocations(vWatershed == "Gulkana");

mCanwellData = mFullData(vWatershed == "Canwell", :);
mCanwellDates = vSampleDates(vWatershed == "Canwell");
mCanwellSites = vSampleLocations(vWatershed == "Canwell");

mCastnerData = mFullData(vWatershed == "Castner", :);
mCastnerDates = vSampleDates(vWatershed == "Castner");
mCastnerSites = vSampleLocations(vWatershed == "Castner");

%% Spatial figures: Matanuska Watershed
% Extract month and year from dates
columnNames = vColumnLabels;
folderName = 'U:/GoA plots/NewPlots';
years = year(mMatDates); % Extract years from dates

% Separate data for 2022 and 2023
data2022 = mMatData(years == 2022, :);
data2023 = mMatData(years == 2023, :);

% Define the months for which to create lines
targetMonths = [5, 6, 7, 8, 9];
colors = ['r', 'g', 'b', 'c', 'm'];
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
                plot(xCoords, yCoords, 'o-', 'Color', colors(i), 'DisplayName', monthLabels{i});
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
            ylabel('Concentration (mg/kg)');
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
    formatFileName = "spatial_mat_%s.jpg";
    fileName = sprintf(formatFileName, sInputChar);
    fullFilePath = fullfile(folderName, fileName);
    saveas(gcf, fullFilePath, 'jpg');
end