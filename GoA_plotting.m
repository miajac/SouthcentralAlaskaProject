%% Plots
clear all, close all, clc
[sFile, sPath] = uigetfile('*.xlsx', 'Select Database File');
sFullPath = fullfile(sPath, sFile);

tData = readtable(sFullPath);

%% 
vColumnLabels = tData.Properties.VariableNames(20:end);
vColumnLabelsArray = cellstr(vColumnLabels);
mFullData = table2array(tData(:, 20:end));

vSampleDates = datetime(convertStringsToChars(string(table2cell(tData(:, 13)))));
vSampleLocations = string(table2cell(tData(:, 3))); 


% Knik River 1 elemental data and sample dates
mKR1Data = mFullData(vSampleLocations == "Knik River 1", :);
mKR1Dates = vSampleDates(vSampleLocations == "Knik River 1");

% Knik River 2 elemental data and sample dates
mKR2Data = mFullData(vSampleLocations == "Knik River 2", :);
mKR2Dates = vSampleDates(vSampleLocations == "Knik River 2");

% Knik River 3 elemental data and sample dates
mKR3Data = mFullData(vSampleLocations == "Knik River 3", :);
mKR3Dates = vSampleDates(vSampleLocations == "Knik River 3");

% Knik River 4 elemental data and sample dates
mKR4Data = mFullData(vSampleLocations == "Knik River 4", :);
mKR4Dates = vSampleDates(vSampleLocations == "Knik River 4");

% Matanuska River 0.5 elemental data and sample dates
mMRHalfData = mFullData(vSampleLocations == "Matanuska River 0.5", :);
mMRHalfDates = vSampleDates(vSampleLocations == "Matanuska River 0.5");

% Matanuska River 1 elemental data and sample dates
mMR1Data = mFullData(vSampleLocations == "Matanuska River 1", :);
mMR1Dates = vSampleDates(vSampleLocations == "Matanuska River 1");

% Matanuska River 2 elemental data and sample dates
mMR2Data = mFullData(vSampleLocations == "Matanuska River 2", :);
mMR2Dates = vSampleDates(vSampleLocations == "Matanuska River 2");

% Matanuska River 3 elemental data and sample dates
mMR3Data = mFullData(vSampleLocations == "Matanuska River 3", :);
mMR3Dates = vSampleDates(vSampleLocations == "Matanuska River 3");

% Matanuska River 4 elemental data and sample dates
mMR4Data = mFullData(vSampleLocations == "Matanuska River 4", :);
mMR4Dates = vSampleDates(vSampleLocations == "Matanuska River 4");

% Matanuska River 5 elemental data and sample dates
mMR5Data = mFullData(vSampleLocations == "Matanuska River 5", :);
mMR5Dates = vSampleDates(vSampleLocations == "Matanuska River 5");

% Little Susitna River 1 elemental data and samples dates
mLS1Data = mFullData(vSampleLocations == "Little Susitna River 1", :);
mLS1Dates = vSampleDates(vSampleLocations == "Little Susitna River 1");

% Little Susitna River 1.5 elemental data and samples dates
mLS15Data = mFullData(vSampleLocations == "Little Susitna River 1.5", :);
mLS15Dates = vSampleDates(vSampleLocations == "Little Susitna River 1.5");

% Little Susitna River 2 elemental data and samples dates
mLS2Data = mFullData(vSampleLocations == "Little Susitna River 2", :);
mLS2Dates = vSampleDates(vSampleLocations == "Little Susitna River 2");

% Little Susitna River 3 elemental data and samples dates
mLS3Data = mFullData(vSampleLocations == "Little Susitna River 3", :);
mLS3Dates = vSampleDates(vSampleLocations == "Little Susitna River 3");

% Little Susitna River 4 elemental data and samples dates
mLS4Data = mFullData(vSampleLocations == "Little Susitna River 4", :);
mLS4Dates = vSampleDates(vSampleLocations == "Little Susitna River 4");

% Moose Creek elemental data and samples dates
mMCData = mFullData(vSampleLocations == "Moose Creek", :);
mMCDates = vSampleDates(vSampleLocations == "Moose Creek");

% Gulkana 1 elemental data and samples dates
mG1Data = mFullData(vSampleLocations == "Gulkana 1", :);
mG1Dates = vSampleDates(vSampleLocations == "Gulkana 1");

% Gulkana 2 elemental data and samples dates
mG2Data = mFullData(vSampleLocations == "Gulkana 2", :);
mG2Dates = vSampleDates(vSampleLocations == "Gulkana 2");

% Gulkana 3 elemental data and samples dates
mG3Data = mFullData(vSampleLocations == "Gulkana 3", :);
mG3Dates = vSampleDates(vSampleLocations == "Gulkana 3");

% Castner 1 elemental data and samples dates
mCT1Data = mFullData(vSampleLocations == "Castner 1", :);
mCT1Dates = vSampleDates(vSampleLocations == "Castner 1");

% Castner 2 elemental data and samples dates
mCT2Data = mFullData(vSampleLocations == "Castner 2", :);
mCT2Dates = vSampleDates(vSampleLocations == "Castner 2");

% Canwell 1 elemental data and samples dates
mCW1Data = mFullData(vSampleLocations == "Canwell 1", :);
mCW1Dates = vSampleDates(vSampleLocations == "Canwell 1");

% Canwell 2 elemental data and samples dates
mCW2Data = mFullData(vSampleLocations == "Canwell 2", :);
mCW2Dates = vSampleDates(vSampleLocations == "Canwell 2");

%% Define watersheds
vWatershedLocations = string(table2cell(tData(:, 5))); 
mMatData = mFullData(vWatershedLocations == "Matanuska",:);
mMatDates = vSampleDates(vWatershedLocations == "Matanuska");
mMatSites = vSampleLocations(vWatershedLocations == "Matanuska");

mKnikData = mFullData(vWatershedLocations == "Knik",:);
mKnikDates = vSampleDates(vWatershedLocations == "Knik");
mKnikSites = vSampleLocations(vWatershedLocations == "Knik");

mLSData = mFullData(vWatershedLocations == "Little Susitna",:);
mLSDates = vSampleDates(vWatershedLocations == "Little Susitna");
mLSSites = vSampleLocations(vWatershedLocations == "Little Susitna");

mGulkanaData = mFullData(vWatershedLocations == "Gulkana",:);
mGulkanaDates = vSampleDates(vWatershedLocations == "Gulkana");
mGulkanaSites = vSampleLocations(vWatershedLocations == "Gulkana");

mCanwellData = mFullData(vWatershedLocations == "Canwell",:);
mCanwellDates = vSampleDates(vWatershedLocations == "Canwell");
mCanwellSites = vSampleLocations(vWatershedLocations == "Canwell");

mCastnerData = mFullData(vWatershedLocations == "Castner",:);
mCastnerDates = vSampleDates(vWatershedLocations == "Castner");
mCastnerSites = vSampleLocations(vWatershedLocations == "Castner");

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
%% Time-series figure: USGS
figure;
plot(mKR3Dates, mKR3Data(:,dPosition), 'o-', 'DisplayName', 'KR3','Color','k','MarkerFaceColor',[0.35,0.70,0.90]);
hold on;
plot(mMR4Dates, mMR4Data(:,dPosition), 'o-', 'DisplayName', 'MR4','Color','k','MarkerFaceColor',[0.95,0.90,0.25]);
plot(mMCDates, mMCData(:,dPosition), 'o-', 'DisplayName', 'Moose Creek','Color','k','MarkerFaceColor',[0.90,0.60,0]);
plot(mLS2Dates, mLS2Data(:,dPosition),'o-', 'DisplayName', 'LS2','Color','k','MarkerFaceColor',[0,0.60,0.50])
xlabel('Time');  % Label for x-axis
ylabel('Concentration in mg/L');  % Label for y-axis
title(sInput,'in mg/L');
legend('Location', 'Best');
grid on;
hold off;
%% Time-series figure: Matanuska Watershed
figure;
plot(mMRHalfDates, mMRHalfData(:,dPosition), 'o-', 'DisplayName', 'MR0.5');
hold on;
plot(mMR1Dates, mMR1Data(:,dPosition), 'o-', 'DisplayName', 'MR1');
plot(mMR2Dates, mMR2Data(:,dPosition), 'o-', 'DisplayName', 'MR2');
plot(mMR3Dates, mMR3Data(:,dPosition),'o-', 'DisplayName', 'MR3');
plot(mMR4Dates, mMR4Data(:,dPosition),'o-', 'DisplayName', 'MR4');
plot(mMR5Dates, mMR5Data(:,dPosition),'o-', 'DisplayName', 'MR5')
xlabel('Time');  % Label for x-axis
ylabel('Concentration in mg/L');  % Label for y-axis
title(sInput,'in mg/L');
legend('Location', 'Best');
grid on;
hold off;
%% Time-series figure: Knik Watershed
figure;
plot(mKR1Dates, mKR1Data(:,dPosition), 'o-', 'DisplayName', 'KR1');
hold on;
plot(mKR2Dates, mKR2Data(:,dPosition), 'o-', 'DisplayName', 'KR2');
plot(mKR3Dates, mKR3Data(:,dPosition), 'o-', 'DisplayName', 'KR3');
plot(mKR4Dates, mKR4Data(:,dPosition),'o-', 'DisplayName', 'KR4')
xlabel('Time');  % Label for x-axis
ylabel('Concentration in mg/L');  % Label for y-axis
title(sInput,'in mg/L');
legend('Location', 'Best');
grid on;
hold off;
%% Time-series figure: Little Susitna Watershed
figure;
plot(mLS1Dates, mLS1Data(:,dPosition), 'o-', 'DisplayName', 'LS1');
hold on;
plot(mLS15Dates, mLS15Data(:,dPosition), 'o-', 'DisplayName', 'LS1.5');
plot(mLS2Dates, mLS2Data(:,dPosition), 'o-', 'DisplayName', 'LS2');
plot(mLS3Dates, mLS3Data(:,dPosition),'o-', 'DisplayName', 'LS3');
plot(mLS4Dates, mLS4Data(:,dPosition),'o-', 'DisplayName', 'LS4')
xlabel('Time');  % Label for x-axis
ylabel('Concentration in mg/L');  % Label for y-axis
title(sInput,'in mg/L');
legend('Location', 'Best');
grid on;
hold off;
%% Time-series figure: Canwell Watershed
figure;
plot(mCW1Dates, mCW1Data(:,dPosition), 'o-', 'DisplayName', 'Canwell 1');
hold on;
plot(mCW2Dates, mCW2Data(:,dPosition), 'o-', 'DisplayName', 'Canwell 2')
xlabel('Time');  % Label for x-axis
ylabel('Concentration in mg/L');  % Label for y-axis
title(sInput,'in mg/L');
legend('Location', 'Best');
grid on;
hold off;
%% Time-series figure: Castner Watershed
figure;
plot(mCT1Dates, mCT1Data(:,dPosition), 'o-', 'DisplayName', 'Castner 1');
hold on;
plot(mCT2Dates, mCT2Data(:,dPosition), 'o-', 'DisplayName', 'Castner 2')
xlabel('Time');  % Label for x-axis
ylabel('Concentration in mg/L');  % Label for y-axis
title(sInput,'in mg/L');
legend('Location', 'Best');
grid on;
hold off;
%% Time-series figure: Gulkana Watershed
figure;
plot(mG1Dates, mG1Data(:,dPosition), 'o-', 'DisplayName', 'Gulkana 1');
hold on;
plot(mG2Dates, mG2Data(:,dPosition), 'o-', 'DisplayName', 'Gulkana 2');
plot(mG3Dates, mG3Data(:,dPosition), 'o-', 'DisplayName', 'Gulkana 3')
xlabel('Time');  % Label for x-axis
ylabel('Concentration in mg/L');  % Label for y-axis
title(sInput,'in mg/L');
legend('Location', 'Best');
grid on;
hold off;
%% Spatial figures: Matanuska Watershed
% Extract month from dates
months = month(mMatDates);

% Define the months for which to create lines
targetMonths = [5, 6, 7, 8, 9];
colors = ['r', 'g', 'b', 'c', 'm'];
monthLabels = {'May', 'June', 'July', 'August', 'September'};

% Get unique site names and sort them alphabetically
uniqueSites = unique(mMatSites);
sortedSites = sort(uniqueSites);

% Create a figure
figure;

for i = 1:numel(targetMonths)
    % Initialize variables to store plot coordinates
    xCoords = [];
    yCoords = [];
    
    % Iterate through each unique site
    for j = 1:numel(sortedSites)
        currentSite = sortedSites{j};
        monthData = mMatData(strcmp(mMatSites, currentSite) & months == targetMonths(i), dPosition);
        
        % Replace NaN values with zeros
        monthData(isnan(monthData)) = 0;
        
        % Store x and y coordinates
        xCoords = [xCoords, j];
        yCoords = [yCoords, monthData];
    end
    
    % Connect points from the same month with a line
    if ~isempty(xCoords)
        plot(xCoords, yCoords, 'o-', 'Color', colors(i), 'DisplayName', monthLabels{i});
        hold on;
    end
end

% Customize the plot
xlabel('Site Names');
ylabel('Concentration');
xticks(1:numel(sortedSites));
xticklabels(sortedSites);  % Use sorted site names

% Add legend manually with proper coloring
legendEntries = cell(numel(targetMonths), 1);
for i = 1:numel(targetMonths)
    legendEntries{i} = sprintf('%s', monthLabels{i});
end

% Display the legend with color coordination
legend(legendEntries, 'Location', 'eastoutside', 'TextColor', 'k', 'Color', 'w');

title('Concentration vs. Site Names for Specific Element');
hold off;
%% Spatial figures: Knik Watershed
% Extract month from dates
months = month(mKnikDates);

% Define the months for which to create lines
targetMonths = [5, 6, 7, 8, 9];
colors = ['r', 'g', 'b', 'c', 'm'];
monthLabels = {'May', 'June', 'July', 'August', 'September'};

% Get unique site names and sort them alphabetically
uniqueSites = unique(mKnikSites);
sortedSites = sort(uniqueSites);

% Create a figure
figure;

for i = 1:numel(targetMonths)
    % Initialize variables to store plot coordinates
    xCoords = [];
    yCoords = [];
    
    % Iterate through each unique site
    for j = 1:numel(sortedSites)
        currentSite = sortedSites{j};
        monthData = mKnikData(strcmp(mKnikSites, currentSite) & months == targetMonths(i), dPosition);
        
        % Replace NaN values with zeros
        monthData(isnan(monthData)) = 0;
        
        % Store x and y coordinates
        if ~isempty(monthData)
            xCoords = [xCoords, j];
            yCoords = [yCoords, mean(monthData)]; % Use mean to handle different lengths
        end
    end
    
    % Connect points from the same month with a line
    if ~isempty(xCoords)
        plot(xCoords, yCoords, 'o-', 'Color', colors(i), 'DisplayName', monthLabels{i});
        hold on;
    end
end

% Customize the plot
xlabel('Site Names');
ylabel('Mean Concentration'); % Adjusted ylabel
xticks(1:numel(sortedSites));
xticklabels(sortedSites);  % Use sorted site names

% Add legend manually with proper coloring
legendEntries = cell(numel(targetMonths), 1);
for i = 1:numel(targetMonths)
    legendEntries{i} = sprintf('%s', monthLabels{i});
end

% Display the legend with color coordination
legend(legendEntries, 'Location', 'eastoutside', 'TextColor', 'k', 'Color', 'w');
title(sInput,'in mg/L');
hold off;
%% Spatial figures: Little Susitna Watershed
% Extract month from dates
months = month(mLSDates);

% Define the months for which to create lines
targetMonths = [5, 6, 7, 8, 9];
colors = ['r', 'g', 'b', 'c', 'm'];
monthLabels = {'May', 'June', 'July', 'August', 'September'};

% Get unique site names and sort them alphabetically
uniqueSites = unique(mLSSites);
sortedSites = sort(uniqueSites);

% Create a figure
figure;

for i = 1:numel(targetMonths) % Initialize variables to store plot coordinates
    xCoords = [];
    yCoords = [];
    
    for j = 1:numel(sortedSites)  % Iterate through each unique site
        currentSite = sortedSites{j};
        monthData = mLSData(strcmp(mLSSites, currentSite) & months == targetMonths(i), dPosition);
        
        monthData(isnan(monthData)) = 0;% Replace NaN values with zeros
        
        if ~isempty(monthData) 
            xCoords = [xCoords, j];
            yCoords = [yCoords, mean(monthData)]; % Use mean to handle different lengths
        end
    end
    
    if ~isempty(xCoords) % Connect points from the same month with a line
        plot(xCoords, yCoords, 'o-', 'Color', colors(i), 'DisplayName', monthLabels{i});
        hold on;
    end
end

% Customize the plot
xlabel('Site Names');
ylabel('Mean Concentration'); 
xticks(1:numel(sortedSites));
xticklabels(sortedSites);  

legendEntries = cell(numel(targetMonths), 1);
for i = 1:numel(targetMonths)
    legendEntries{i} = sprintf('%s', monthLabels{i});
end

legend(legendEntries, 'Location', 'eastoutside', 'TextColor', 'k', 'Color', 'w');
title(sInput,'in mg/L');
hold off;
%% Spatial figures: Gulkana Watershed
months = month(mGulkanaDates);

targetMonths = [5, 6, 7, 9];
colors = ['r', 'g', 'b', 'c'];
monthLabels = {'May', 'June', 'July', 'September'};

uniqueSites = unique(mGulkanaSites);
sortedSites = sort(uniqueSites);

% Create a figure
figure;

for i = 1:numel(targetMonths) % Initialize variables to store plot coordinates
    xCoords = [];
    yCoords = [];
    
    for j = 1:numel(sortedSites)  % Iterate through each unique site
        currentSite = sortedSites{j};
        monthData = mGulkanaData(strcmp(mGulkanaSites, currentSite) & months == targetMonths(i), dPosition);
        
        monthData(isnan(monthData)) = 0;% Replace NaN values with zeros
        
        if ~isempty(monthData) 
            xCoords = [xCoords, j];
            yCoords = [yCoords, mean(monthData)]; % Use mean to handle different lengths
        end
    end
    
    if ~isempty(xCoords) % Connect points from the same month with a line
        plot(xCoords, yCoords, 'o-', 'Color', colors(i), 'DisplayName', monthLabels{i});
        hold on;
    end
end

% Customize the plot
xlabel('Site Names');
ylabel('Mean Concentration'); 
xticks(1:numel(sortedSites));
xticklabels(sortedSites);  

legendEntries = cell(numel(targetMonths), 1);
for i = 1:numel(targetMonths)
    legendEntries{i} = sprintf('%s', monthLabels{i});
end

legend(legendEntries, 'Location', 'eastoutside', 'TextColor', 'k', 'Color', 'w');
title(sInput,'in mg/L');
hold off;
%% Spatial figures: Canwell Watershed
months = month(mCanwellDates);

targetMonths = [5, 6, 7, 9];
colors = ['r', 'g', 'b', 'c'];
monthLabels = {'May', 'June', 'July', 'September'};

uniqueSites = unique(mCanwellSites);
sortedSites = sort(uniqueSites);

% Create a figure
figure;

for i = 1:numel(targetMonths) % Initialize variables to store plot coordinates
    xCoords = [];
    yCoords = [];
    
    for j = 1:numel(sortedSites)  % Iterate through each unique site
        currentSite = sortedSites{j};
        monthData = mCanwellData(strcmp(mCanwellSites, currentSite) & months == targetMonths(i), dPosition);
        
        monthData(isnan(monthData)) = 0;% Replace NaN values with zeros
        
        if ~isempty(monthData) 
            xCoords = [xCoords, j];
            yCoords = [yCoords, mean(monthData)]; % Use mean to handle different lengths
        end
    end
    
    if ~isempty(xCoords) 
        plot(xCoords, yCoords, 'o-', 'Color', colors(i), 'DisplayName', monthLabels{i});
        hold on;
    end
end

xlabel('Site Names');
ylabel('Mean Concentration'); 
xticks(1:numel(sortedSites));
xticklabels(sortedSites);  

legendEntries = cell(numel(targetMonths), 1);
for i = 1:numel(targetMonths)
    legendEntries{i} = sprintf('%s', monthLabels{i});
end

legend(legendEntries, 'Location', 'eastoutside', 'TextColor', 'k', 'Color', 'w');
title(sInput,'in mg/L');
hold off;
%% Spatial figures: Castner Watershed
months = month(mCastnerDates);

targetMonths = [5, 6, 7, 9];
colors = ['r', 'g', 'b', 'c'];
monthLabels = {'May', 'June', 'July', 'September'};

uniqueSites = unique(mCastnerSites);
sortedSites = sort(uniqueSites);

% Create a figure
figure;

for i = 1:numel(targetMonths) 
    xCoords = [];
    yCoords = [];
    
    for j = 1:numel(sortedSites)
        currentSite = sortedSites{j};
        monthData = mCastnerData(strcmp(mCastnerSites, currentSite) & months == targetMonths(i), dPosition);
        
        monthData(isnan(monthData)) = 0;
        
        if ~isempty(monthData) 
            xCoords = [xCoords, j];
            yCoords = [yCoords, mean(monthData)]; 
        end
    end
    
    if ~isempty(xCoords) 
        plot(xCoords, yCoords, 'o-', 'Color', colors(i), 'DisplayName', monthLabels{i});
        hold on;
    end
end

% Customize the plot
xlabel('Site Names');
ylabel('Mean Concentration'); 
xticks(1:numel(sortedSites));
xticklabels(sortedSites);  

legendEntries = cell(numel(targetMonths), 1);
for i = 1:numel(targetMonths)
    legendEntries{i} = sprintf('%s', monthLabels{i});
end

legend(legendEntries, 'Location', 'eastoutside', 'TextColor', 'k', 'Color', 'w');
title(sInput,'in mg/L');
hold off;
%% Stable isotope comparison: USGS
figure;
plot(mKR3Data(:,9), mKR3Data(:,10), 'o', 'DisplayName', 'KR3','Color','k','MarkerFaceColor',[0.35,0.70,0.90]);
hold on;
plot(mMR4Data(:,9), mMR4Data(:,10), 'o', 'DisplayName', 'MR4','Color','k','MarkerFaceColor',[0.95,0.90,0.25]);
plot(mMCData(:,9), mMCData(:,10), 'o', 'DisplayName', 'Moose Creek','Color','k','MarkerFaceColor',[0.90,0.60,0]);
plot(mLS2Data(:,9), mLS2Data(:,10),'o', 'DisplayName', 'LS2','Color','k','MarkerFaceColor',[0,0.60,0.50])
xlabel('Average DO'); 
ylabel('Average DD');  
title('USGS stable isotope comparison');
legend('Location', 'Best');
grid on;
hold off;
%% Stable isotope comparison: all
figure;
plot(mKR1Data(:,9), mKR1Data(:,10), 'o', 'DisplayName', 'KR1','Color', 'b');
hold on;
plot(mKR2Data(:,9), mKR2Data(:,10), 'o', 'DisplayName', 'KR2','Color', 'r');
plot(mKR3Data(:,9), mKR3Data(:,10), 'o', 'DisplayName', 'KR3','Color', 'g');
plot(mKR4Data(:,9), mKR4Data(:,10),'o', 'DisplayName', 'KR4','Color', 'y');
plot(mMRHalfData(:,9), mMRHalfData(:,10), '+', 'DisplayName', 'MR0.5','Color', 'b');
plot(mMR1Data(:,9), mMR1Data(:,10), '+', 'DisplayName', 'MR1','Color', 'r');
plot(mMR2Data(:,9), mMR2Data(:,10),'+', 'DisplayName', 'MR2','Color', 'g');
plot(mMR3Data(:,9), mMR3Data(:,10), '+', 'DisplayName', 'MR3','Color', 'y');
plot(mMR4Data(:,9), mMR4Data(:,10), '+', 'DisplayName', 'MR4','Color', 'm');
plot(mMR5Data(:,9), mMR5Data(:,10),'+', 'DisplayName', 'MR5','Color', 'k');
plot(mMCData(:,9), mMCData(:,10), '+', 'DisplayName', 'MC','Color', [0.75, 0.95, 0.0]);
plot(mLS1Data(:,9), mLS1Data(:,10), '*', 'DisplayName', 'LS1','Color', 'b');
plot(mLS15Data(:,9), mLS15Data(:,10),'*', 'DisplayName', 'LS1.5','Color', 'r');
plot(mLS2Data(:,9), mLS2Data(:,10), '*', 'DisplayName', 'LS2','Color', 'g');
plot(mLS3Data(:,9), mLS3Data(:,10),'*', 'DisplayName', 'LS3','Color', 'y');
plot(mLS4Data(:,9), mLS4Data(:,10), '*', 'DisplayName', 'LS4','Color', 'm');
plot(mG1Data(:,9), mG1Data(:,10),'x', 'DisplayName', 'G1','Color', 'b');
plot(mG2Data(:,9), mG2Data(:,10),'x', 'DisplayName', 'G2','Color', 'r');
plot(mG3Data(:,9), mG3Data(:,10),'x', 'DisplayName', 'G3','Color', 'g');
plot(mCT1Data(:,9), mCT1Data(:,10),'s', 'DisplayName', 'CT1','Color', 'b');
plot(mCT2Data(:,9), mCT2Data(:,10),'s', 'DisplayName', 'CT2','Color', 'r');
plot(mCW1Data(:,9), mCW1Data(:,10),'d', 'DisplayName', 'CW1','Color', 'b');
plot(mCW2Data(:,9), mCW2Data(:,10),'d', 'DisplayName', 'CW2','Color', 'r');
xlabel('Average DO'); 
ylabel('Average DD');  
title('Stable isotope comparison');
legend('Location', 'Best');
grid on;
hold off;
%% Stable Isotopes over time: USGS
figure;
plot(mKR3Dates, mKR3Data(:,9), 'o-', 'DisplayName', 'KR3 average DO','Color','k','MarkerFaceColor',[0.35,0.70,0.90]);
hold on;
plot(mMR4Dates, mMR4Data(:,9), 'o-', 'DisplayName', 'MR4 average DO','Color','k','MarkerFaceColor',[0.95,0.90,0.25]);
plot(mMCDates, mMCData(:,9), 'o-', 'DisplayName', 'Moose Creek average DO','Color','k','MarkerFaceColor',[0.90,0.60,0]);
plot(mLS2Dates, mLS2Data(:,9),'o-', 'DisplayName', 'LS2 average DO','Color','k','MarkerFaceColor',[0,0.60,0.50]);
plot(mKR3Dates, mKR3Data(:,10), '--o', 'DisplayName', 'KR3 average DD','Color','k','MarkerFaceColor',[0.35,0.70,0.90]);
plot(mMR4Dates, mMR4Data(:,10), '--o', 'DisplayName', 'MR4 average DD','Color','k','MarkerFaceColor',[0.95,0.90,0.25]);
plot(mMCDates, mMCData(:,10), '--o', 'DisplayName', 'Moose Creek average DD','Color','k','MarkerFaceColor',[0.90,0.60,0]);
plot(mLS2Dates, mLS2Data(:,10),'--o', 'DisplayName', 'LS2 average DD','Color','k','MarkerFaceColor',[0,0.60,0.50]);
xlabel('Time');  
ylabel('Stable isotope value');  
title('Stable isotope variations over time');
legend('Location', 'eastoutside');
grid on;
hold off;