%% spatial
% Load and Prepare Data
clear all, close all, clc 

% Load Excel file
[sFile, sPath] = uigetfile('*.xlsx', 'Select Database File');
sFullPath = fullfile(sPath, sFile);

tData = readtable(sFullPath);

% Extract relevant columns
vColumnLabels = tData.Properties.VariableNames(24:end);
vColumnLabelsArray = cellstr(vColumnLabels);
mFullData = table2array(tData(1:end, 24:end));

vSampleDates = datetime(convertStringsToChars(string(table2cell(tData(:, 14)))));
vSampleLocations = string(table2cell(tData(:, 2)));
vWatershed = string(table2cell(tData(:,5))); 
vSampleYear = table2array(tData(:,16)); 
vDistGl = table2array(tData(:,8)); % Distance from glacier

%% Define months and colors
targetMonths = [5, 6, 7, 8, 9];
colors = {[0.8, 0.6, 0.7], [0.9, 0.6, 0], [0.95, 0.9, 0.25], [0, 0.6, 0.5], [0.35, 0.7, 0.9]};
monthLabels = {'May', 'June', 'July', 'August', 'September'};

% Define the dPosition values to iterate over
dPositions = [29, 16, 3, 47];

% Define save folder
folderName = 'U:/GoA plots/NewPlots';

disp('Data successfully loaded and preprocessed.');
 
%% Matanuska Subplots 
columnNames = vColumnLabels;
folderName = 'U:/GoA plots/NewPlots';
years = year(vSampleDates(vWatershed == "Matanuska"));

figure;
set(gcf, 'Position', [100, 100, 600, 800]);

for dpIndex = 1:numel(dPositions)
    dPosition = dPositions(dpIndex);
    sInput = columnNames{dPosition};
    yLimits = []; 

    if contains(sInput, 'HCO3') 
        sInput = strrep(sInput, 'HCO3', 'HCO_3^-');
    end
    if contains(sInput, 'Ca') && ~contains(sInput, 'Ca^{2+}')
        sInput = strrep(sInput, 'Ca', 'Ca^{2+}');
    end

    for yearIndex = 1:2
        currentYear = 2022 + (yearIndex - 1);
        yearMask = (years == currentYear);

        currentData = mFullData(vWatershed == "Matanuska", :);
        currentDists = vDistGl(vWatershed == "Matanuska");
        months = month(vSampleDates(vWatershed == "Matanuska"));

        validMask = ~isnan(currentDists) & ~isnan(currentData(:, dPosition)) & yearMask;
        currentDists = currentDists(validMask);
        currentData = currentData(validMask, dPosition);
        months = months(validMask);

        [sortedDists, sortIdx] = sort(currentDists);
        sortedData = currentData(sortIdx);
        sortedMonths = months(sortIdx);

        subplotIndex = (dpIndex - 1) * 2 + yearIndex;
        subplot(4, 2, subplotIndex);
        hold on;

        for i = 1:numel(targetMonths)
            monthMask = sortedMonths == targetMonths(i);
            xCoords = sortedDists(monthMask);
            yCoords = sortedData(monthMask);

            if ~isempty(xCoords) && ~isempty(yCoords)
                plot(xCoords, yCoords, 'o-', 'Color', colors{i}, ...
                    'MarkerFaceColor', colors{i}, 'MarkerEdgeColor', 'k', ...
                    'DisplayName', monthLabels{i});
            end
        end

        if subplotIndex >= 7 % Only show x-axis label on the last two subplots
            xlabel('Distance from Glacier (km)');
        end

        if dPosition == 29 || dPosition == 47
            ylabel('Concentration (μg/L)');
        else
            ylabel('Concentration (mg/L)');
        end

        title([sInput ' ' num2str(currentYear)], 'Interpreter', 'tex');

        % if dpIndex == 1 && yearIndex == 1
        %     legend('Location', 'northeast');
        % end

        yLimits = [yLimits, ylim];
    end

    y_min = min(yLimits);
    y_max = max(yLimits);
    for yearIndex = 1:2
        subplotIndex = (dpIndex - 1) * 2 + yearIndex;
        subplot(4, 2, subplotIndex);
        ylim([y_min, y_max]);
    end
end

saveas(gcf, fullfile(folderName, 'spatial_mat_subplot.svg'), 'svg');

%% Knik Subplots 
columnNames = vColumnLabels;
folderName = 'U:/GoA plots/NewPlots';
years = year(vSampleDates(vWatershed == "Knik"));

figure;
set(gcf, 'Position', [100, 100, 600, 800]);

for dpIndex = 1:numel(dPositions)
    dPosition = dPositions(dpIndex);
    sInput = columnNames{dPosition};
    yLimits = []; 

    if contains(sInput, 'HCO3') 
        sInput = strrep(sInput, 'HCO3', 'HCO_3^-');
    end
    if contains(sInput, 'Ca') && ~contains(sInput, 'Ca^{2+}')
        sInput = strrep(sInput, 'Ca', 'Ca^{2+}');
    end

    for yearIndex = 1:2
        currentYear = 2022 + (yearIndex - 1);
        yearMask = (years == currentYear);

        currentData = mFullData(vWatershed == "Knik", :);
        currentDists = vDistGl(vWatershed == "Knik");
        months = month(vSampleDates(vWatershed == "Knik"));

        validMask = ~isnan(currentDists) & ~isnan(currentData(:, dPosition)) & yearMask;
        currentDists = currentDists(validMask);
        currentData = currentData(validMask, dPosition);
        months = months(validMask);

        [sortedDists, sortIdx] = sort(currentDists);
        sortedData = currentData(sortIdx);
        sortedMonths = months(sortIdx);

        subplotIndex = (dpIndex - 1) * 2 + yearIndex;
        subplot(4, 2, subplotIndex);
        hold on;

        for i = 1:numel(targetMonths)
            monthMask = sortedMonths == targetMonths(i);
            xCoords = sortedDists(monthMask);
            yCoords = sortedData(monthMask);

            if ~isempty(xCoords) && ~isempty(yCoords)
                plot(xCoords, yCoords, 'o-', 'Color', colors{i}, ...
                    'MarkerFaceColor', colors{i}, 'MarkerEdgeColor', 'k', ...
                    'DisplayName', monthLabels{i});
            end
        end

        if subplotIndex >= 7 % Only show x-axis label on the last two subplots
            xlabel('Distance from Glacier (km)');
        end

        if dPosition == 29 || dPosition == 47
            ylabel('Concentration (μg/L)');
        else
            ylabel('Concentration (mg/L)');
        end

        title([sInput ' ' num2str(currentYear)], 'Interpreter', 'tex');

        % if dpIndex == 1 && yearIndex == 1
        %     legend('Location', 'northeast');
        % end

        yLimits = [yLimits, ylim];
    end

    y_min = min(yLimits);
    y_max = max(yLimits);
    for yearIndex = 1:2
        subplotIndex = (dpIndex - 1) * 2 + yearIndex;
        subplot(4, 2, subplotIndex);
        ylim([y_min, y_max]);
    end
end

saveas(gcf, fullfile(folderName, 'spatial_knik_subplot.svg'), 'svg');

%% Little Susitna Subplots 
columnNames = vColumnLabels;
folderName = 'U:/GoA plots/NewPlots';
years = year(vSampleDates(vWatershed == "Little Susitna"));

figure;
set(gcf, 'Position', [100, 100, 600, 800]);

for dpIndex = 1:numel(dPositions)
    dPosition = dPositions(dpIndex);
    sInput = columnNames{dPosition};
    yLimits = []; 

    if contains(sInput, 'HCO3') 
        sInput = strrep(sInput, 'HCO3', 'HCO_3^-');
    end
    if contains(sInput, 'Ca') && ~contains(sInput, 'Ca^{2+}')
        sInput = strrep(sInput, 'Ca', 'Ca^{2+}');
    end

    for yearIndex = 1:2
        currentYear = 2022 + (yearIndex - 1);
        yearMask = (years == currentYear);

        currentData = mFullData(vWatershed == "Little Susitna", :);
        currentDists = vDistGl(vWatershed == "Little Susitna");
        months = month(vSampleDates(vWatershed == "Little Susitna"));

        validMask = ~isnan(currentDists) & ~isnan(currentData(:, dPosition)) & yearMask;
        currentDists = currentDists(validMask);
        currentData = currentData(validMask, dPosition);
        months = months(validMask);

        [sortedDists, sortIdx] = sort(currentDists);
        sortedData = currentData(sortIdx);
        sortedMonths = months(sortIdx);

        subplotIndex = (dpIndex - 1) * 2 + yearIndex;
        subplot(4, 2, subplotIndex);
        hold on;

        for i = 1:numel(targetMonths)
            monthMask = sortedMonths == targetMonths(i);
            xCoords = sortedDists(monthMask);
            yCoords = sortedData(monthMask);

            if ~isempty(xCoords) && ~isempty(yCoords)
                plot(xCoords, yCoords, 'o-', 'Color', colors{i}, ...
                    'MarkerFaceColor', colors{i}, 'MarkerEdgeColor', 'k', ...
                    'DisplayName', monthLabels{i});
            end
        end

        if subplotIndex >= 7 % Only show x-axis label on the last two subplots
            xlabel('Distance from Glacier (km)');
        end

        if dPosition == 29 || dPosition == 47
            ylabel('Concentration (μg/L)');
        else
            ylabel('Concentration (mg/L)');
        end

        title([sInput ' ' num2str(currentYear)], 'Interpreter', 'tex');

        %if dpIndex == 1 && yearIndex == 1
            %legend('Location', 'northeast');
        %end

        yLimits = [yLimits, ylim];
    end

    y_min = min(yLimits);
    y_max = max(yLimits);
    for yearIndex = 1:2
        subplotIndex = (dpIndex - 1) * 2 + yearIndex;
        subplot(4, 2, subplotIndex);
        ylim([y_min, y_max]);
    end
end

saveas(gcf, fullfile(folderName, 'spatial_ls_subplot.svg'), 'svg');

%% Knik legend seperate %% Knik Subplots 
columnNames = vColumnLabels;
folderName = 'U:/GoA plots/NewPlots';
years = year(vSampleDates(vWatershed == "Knik"));

figure;
set(gcf, 'Position', [100, 100, 600, 800]);

for dpIndex = 1:numel(dPositions)
    dPosition = dPositions(dpIndex);
    sInput = columnNames{dPosition};
    yLimits = []; 

    if contains(sInput, 'HCO3') 
        sInput = strrep(sInput, 'HCO3', 'HCO_3^-');
    end
    if contains(sInput, 'Ca') && ~contains(sInput, 'Ca^{+2}')
        sInput = strrep(sInput, 'Ca', 'Ca^{2+}');
    end

    for yearIndex = 1:2
        currentYear = 2022 + (yearIndex - 1);
        yearMask = (years == currentYear);

        currentData = mFullData(vWatershed == "Knik", :);
        currentDists = vDistGl(vWatershed == "Knik");
        months = month(vSampleDates(vWatershed == "Knik"));

        validMask = ~isnan(currentDists) & ~isnan(currentData(:, dPosition)) & yearMask;
        currentDists = currentDists(validMask);
        currentData = currentData(validMask, dPosition);
        months = months(validMask);

        [sortedDists, sortIdx] = sort(currentDists);
        sortedData = currentData(sortIdx);
        sortedMonths = months(sortIdx);

        subplotIndex = (dpIndex - 1) * 2 + yearIndex;
        subplot(4, 2, subplotIndex);
        hold on;

        for i = 1:numel(targetMonths)
            monthMask = sortedMonths == targetMonths(i);
            xCoords = sortedDists(monthMask);
            yCoords = sortedData(monthMask);

            if ~isempty(xCoords) && ~isempty(yCoords)
                plot(xCoords, yCoords, 'o-', 'Color', colors{i}, ...
                    'MarkerFaceColor', colors{i}, 'MarkerEdgeColor', 'k', ...
                    'DisplayName', monthLabels{i});
            end
        end

        if subplotIndex >= 7 % Only show x-axis label on the last two subplots
            xlabel('Distance from Glacier (km)');
        end

        if dPosition == 29 || dPosition == 47
            ylabel('Concentration (μg/L)');
        else
            ylabel('Concentration (mg/L)');
        end

        title([sInput ' ' num2str(currentYear)], 'Interpreter', 'tex');

        % Commented out the legend here
        % if dpIndex == 1 && yearIndex == 1
        %     legend('Location', 'northeast');
        % end

        yLimits = [yLimits, ylim];
    end

    y_min = min(yLimits);
    y_max = max(yLimits);
    for yearIndex = 1:2
        subplotIndex = (dpIndex - 1) * 2 + yearIndex;
        subplot(4, 2, subplotIndex);
        ylim([y_min, y_max]);
    end
end

saveas(gcf, fullfile(folderName, 'spatial_knik_subplot.svg'), 'svg');

%% Separate Legend Figure
legendFigure = figure;
set(legendFigure, 'Position', [100, 100, 300, 200]);
hold on;

% Dummy plots to create legend
legendHandles = gobjects(1, numel(targetMonths));
for i = 1:numel(targetMonths)
    legendHandles(i) = plot(nan, nan, 'o-', 'Color', colors{i}, ...
        'MarkerFaceColor', colors{i}, 'MarkerEdgeColor', 'k', ...
        'DisplayName', monthLabels{i});
end

legend(legendHandles, monthLabels, 'Location', 'northwest');
axis off; % Hide axes
saveas(legendFigure, fullfile(folderName, 'spatial_legend.svg'), 'svg');
