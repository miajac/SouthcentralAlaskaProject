%% CQ-log transform
clear all, close all, clc %use USGSGaugeDataset
[sFile, sPath] = uigetfile('*.xlsx', 'Select Database File');
sFullPath = fullfile(sPath, sFile);

tData = readtable(sFullPath);
%% 
vColumnLabels = tData.Properties.VariableNames(13:end);
vColumnLabelsArray = cellstr(vColumnLabels);
mFullData = table2array(tData(1:end, 13:end));
vSampleLocations = string(table2cell(tData(:, 3)));
vSampleDates = datetime(convertStringsToChars(string(table2cell(tData(:, 10)))));
vSampleMonths = month(vSampleDates);
%% Knik River 3 elemental data and sample dates
mKR3Data = mFullData(vSampleLocations == "KR3", :);
mMR4Data = mFullData(vSampleLocations == "MR4", :);
mLS2Data = mFullData(vSampleLocations == "LS2", :);
mMCData = mFullData(vSampleLocations == "MC", :);

%% Linear Regression and Power Law Form with Recession Lines (no unit corrections- expand if needed)
columnNames = vColumnLabels;
folderName = 'U:/GoA plots/NewPlots';
for dPosition = 1:numel(columnNames)
    sInput = columnNames{dPosition};
    sTitle = sInput;
    figure;
    set(gcf, 'Position', [100, 100, 800, 450]);
    hold on;
    
    % Define colors and markers for each site
    siteColors = {[0.35,0.70,0.90], [0.95,0.90,0.25], [0.90,0.60,0], [0,0.60,0.50]};
    siteMarkers = {'o', 'o', 'o', 'o'};
    siteLabels = {'K3', 'M4', 'MC', 'LS2'};
    siteData = {mKR3Data, mMR4Data, mMCData, mLS2Data};
    
    legendEntries = {}; % Store legend text
    
    for i = 1:numel(siteData)
        y_data = siteData{i}(:, dPosition);  % Log-transformed concentration
        x_data = siteData{i}(:, 1);  % Log-transformed discharge
        y_data = fillmissing(y_data, 'linear'); % Fill missing values
        
        % Plot data points
        plot(x_data, y_data, siteMarkers{i}, 'DisplayName', siteLabels{i}, ...
            'Color', siteColors{i}, 'MarkerFaceColor', siteColors{i}, 'MarkerEdgeColor', 'k');

        % Perform linear regression
        validIdx = ~isnan(x_data) & ~isnan(y_data); % Remove NaNs
        p = polyfit(x_data(validIdx), y_data(validIdx), 1); % Fit a line
        slope = p(1);
        intercept = p(2);
        r2 = 1 - sum((y_data(validIdx) - polyval(p, x_data(validIdx))).^2) / sum((y_data(validIdx) - mean(y_data(validIdx))).^2);

        % Convert to power law form C = (10^a) * Q^b
        a = 10^intercept;
        b = slope;

        % Generate regression line points
        x_fit = linspace(min(x_data), max(x_data), 100);
        y_fit = polyval(p, x_fit);
        plot(x_fit, y_fit, '-', 'Color', siteColors{i}, 'LineWidth', 1.5, 'HandleVisibility', 'off');

        % Store equation in legend entry
        legendEntries{end+1} = sprintf('%s: C = %.6f * Q^{%.2f}, R^2 = %.3f', siteLabels{i}, a, b, r2);
    end
    
    xlabel('Log discharge (m^{3}/s)'); 
    ylabel(sprintf('Log dissolved %s (log(mg/L))', sInput)); 
    %title(sTitle);
    legend(legendEntries, 'Location', 'eastoutside');
    grid on;
    hold off;
    
    % Save figure
    sInputChar = char(sInput);
    formatFileName = "conc_dis_log_usgs_power_%s.jpg";
    fileName = sprintf(formatFileName, sInputChar);
    fullFilePath = fullfile(folderName, fileName);
    saveas(gcf, fullFilePath, 'jpg');
end

%% Linear Regression and Power Law Form with Recession Lines (with units)
columnNames = vColumnLabels;
folderName = 'U:/GoA plots/NewPlots';
for dPosition = 1:numel(columnNames)
    sInput = columnNames{dPosition};
    sTitle = sInput;
    figure;
    set(gcf, 'Position', [100, 100, 800, 450]);
    hold on;
    
    % Define colors and markers for each site
    siteColors = {[0.35,0.70,0.90], [0.95,0.90,0.25], [0.90,0.60,0], [0,0.60,0.50]};
    siteMarkers = {'o', 'o', 'o', 'o'};
    siteLabels = {'K3', 'M4', 'MC', 'LS2'};
    siteData = {mKR3Data, mMR4Data, mMCData, mLS2Data};
    
    legendEntries = {}; % Store legend text
    
    for i = 1:numel(siteData)
        y_data = siteData{i}(:, dPosition);  % Log-transformed concentration
        x_data = siteData{i}(:, 1);  % Log-transformed discharge
        y_data = fillmissing(y_data, 'linear'); % Fill missing values
        
        % Plot data points
        plot(x_data, y_data, siteMarkers{i}, 'DisplayName', siteLabels{i}, ...
            'Color', siteColors{i}, 'MarkerFaceColor', siteColors{i}, 'MarkerEdgeColor', 'k');

        % Perform linear regression
        validIdx = ~isnan(x_data) & ~isnan(y_data); % Remove NaNs
        p = polyfit(x_data(validIdx), y_data(validIdx), 1); % Fit a line
        slope = p(1);
        intercept = p(2);
        r2 = 1 - sum((y_data(validIdx) - polyval(p, x_data(validIdx))).^2) / sum((y_data(validIdx) - mean(y_data(validIdx))).^2);

        % Convert to power law form C = (10^a) * Q^b
        a = 10^intercept;
        b = slope;

        % Generate regression line points
        x_fit = linspace(min(x_data), max(x_data), 100);
        y_fit = polyval(p, x_fit);
        plot(x_fit, y_fit, '-', 'Color', siteColors{i}, 'LineWidth', 1.5, 'HandleVisibility', 'off');

        % Store equation in legend entry
        legendEntries{end+1} = sprintf('%s: C = %.6f * Q^{%.2f}, R^2 = %.3f', siteLabels{i}, a, b, r2);
    end
    
    % Determine appropriate y-axis label
    if ismember(dPosition, [15,16,19]) || (dPosition >= 22 && dPosition <= 55)
        yUnit = 'μg/L';
    elseif dPosition == 8
        yUnit = 'mV';
    elseif dPosition == 7
        yUnit = 'pH';
    elseif dPosition == 6
        yUnit = 'μS/cm';
    elseif dPosition == 1
        yUnit = 'm^{3}/s';
    elseif dPosition == 2
        yUnit = '{\circ}C'; % Degree symbol formatted for MATLAB
    elseif dPosition == 3
        yUnit = 'mmHg';
    elseif dPosition == 4
        yUnit = '%';
    else
        yUnit = 'mg/L';
    end

    % Convert specific chemical species into proper notation
    switch sInput
        case 'F'
            sFormattedInput = 'F^{-}';
        case 'HCO3'
            sFormattedInput = 'HCO_{3}^{-}';
        case 'Cl'
            sFormattedInput = 'Cl^{-}';
        case 'NO3'
            sFormattedInput = 'NO_{3}^{-}';
        case 'SO4'
            sFormattedInput = 'SO_{4}^{2-}';
        case 'Na'
            sFormattedInput = 'Na^{+}';
        case 'Mg'
            sFormattedInput = 'Mg^{2+}';
        case 'K'
            sFormattedInput = 'K^{+}';
        case 'Ca'
            sFormattedInput = 'Ca^{2+}';
        otherwise
            sFormattedInput = sInput; % Keep original label if no special formatting needed
    end

    xlabel('Log discharge (m^{3}/s)'); 
    ylabel(sprintf('Log dissolved %s (%s)', sFormattedInput, yUnit), 'Interpreter', 'tex'); 

    legend(legendEntries, 'Location', 'eastoutside');
    grid on;
    hold off;
    
    % Save figure
    sInputChar = char(sInput);
    formatFileName = "conc_dis_log_usgs_power_unit_%s.jpg";
    fileName = sprintf(formatFileName, sInputChar);
    fullFilePath = fullfile(folderName, fileName);
    saveas(gcf, fullFilePath, 'jpg');
end
%% Seasonality CQ
% Site info
siteData = {mKR3Data, mMR4Data, mMCData, mLS2Data};
siteLabels = {'K3', 'M4', 'MC', 'LS2'};
siteIDs = {'KR3', 'MR4', 'MC', 'LS2'};
siteColors = {
    [0.35, 0.70, 0.90],  % KR3 = Blue
    [0.95, 0.90, 0.25],  % MR4 = Yellow
    [0.90, 0.60, 0.00],  % MC  = Orange
    [0.00, 0.60, 0.50]   % LS2 = Teal
};

% Month label and marker maps (May–Oct)
monthIDs = [5, 6, 7, 8, 9, 10];
monthMarkers = {'o', 's', '^', 'd', 'p', 'p'};
monthLabels = {'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct'};

folderName = 'U:/GoA plots/NewPlots';
regressionSummary = cell(0, 7);

%% Loop through each element
for dPosition = 1:numel(vColumnLabels)
    sInput = vColumnLabels{dPosition};
    figure;
    set(gcf, 'Position', [100, 100, 900, 500]);
    hold on;

    legendHandles = [];
    legendLabels = [];

    a_vals = nan(1, numel(siteData));
    b_vals = nan(1, numel(siteData));
    r2_vals = nan(1, numel(siteData));
    eqn_strings = cell(1, numel(siteData));  % Updated to cell array

    for i = 1:numel(siteData)
        currentSite = siteIDs{i};
        currentColor = siteColors{i};
        dataMatrix = siteData{i};
        siteIdx = vSampleLocations == currentSite;
        siteMonths = vSampleMonths(siteIdx);
        siteDischarge = dataMatrix(:, 1);
        siteConc = dataMatrix(:, dPosition);
        siteConc = fillmissing(siteConc, 'linear');

        for m = 1:length(monthIDs)
            monthVal = monthIDs(m);
            monthIdx = siteMonths == monthVal;
            if sum(monthIdx) < 2, continue; end

            x = siteDischarge(monthIdx);
            y = siteConc(monthIdx);
            marker = monthMarkers{m};

            h = scatter(x, y, 50, ...
                'Marker', marker, ...
                'MarkerEdgeColor', 'k', ...
                'MarkerFaceColor', currentColor, ...
                'DisplayName', sprintf('%s - %s', siteLabels{i}, monthLabels{m}));

            legendHandles(end+1) = h;
            legendLabels{end+1} = sprintf('%s - %s', siteLabels{i}, monthLabels{m});
        end

        % Regression line (in same site color)
        validIdx = ~isnan(siteDischarge) & ~isnan(siteConc);
        p = polyfit(siteDischarge(validIdx), siteConc(validIdx), 1);
        slope = p(1); intercept = p(2);
        r2 = 1 - sum((siteConc(validIdx) - polyval(p, siteDischarge(validIdx))).^2) / ...
                sum((siteConc(validIdx) - mean(siteConc(validIdx))).^2);
        a = 10^intercept;
        b = slope;

        a_vals(i) = a;
        b_vals(i) = b;
        r2_vals(i) = r2;
        eqn_strings{i} = sprintf('C = %.4f * Q^{%.2f}', a, b);

        x_fit = linspace(min(siteDischarge), max(siteDischarge), 100);
        y_fit = polyval(p, x_fit);
        h_line = plot(x_fit, y_fit, '-', 'Color', currentColor, 'LineWidth', 1.5);
        set(h_line, 'HandleVisibility', 'off');
    end

    % === Dynamic y-axis label ===
    if ismember(dPosition, [15,16,19]) || (dPosition >= 22 && dPosition <= 55)
        yUnit = 'μg/L';
    elseif dPosition == 8
        yUnit = 'mV';
    elseif dPosition == 7
        yUnit = 'pH';
    elseif dPosition == 6
        yUnit = 'μS/cm';
    elseif dPosition == 1
        yUnit = 'm^{3}/s';
    elseif dPosition == 2
        yUnit = '{\circ}C';
    elseif dPosition == 3
        yUnit = 'mmHg';
    elseif dPosition == 4
        yUnit = '%';
    else
        yUnit = 'mg/L';
    end

    % === Format chemical symbol for y-axis ===
    switch sInput
        case 'F';     sFormattedInput = 'F^{-}';
        case 'HCO3';  sFormattedInput = 'HCO_{3}^{-}';
        case 'Cl';    sFormattedInput = 'Cl^{-}';
        case 'NO3';   sFormattedInput = 'NO_{3}^{-}';
        case 'SO4';   sFormattedInput = 'SO_{4}^{2-}';
        case 'Na';    sFormattedInput = 'Na^{+}';
        case 'Mg';    sFormattedInput = 'Mg^{2+}';
        case 'K';     sFormattedInput = 'K^{+}';
        case 'Ca';    sFormattedInput = 'Ca^{2+}';
        otherwise;    sFormattedInput = sInput;
    end

    % === Axis labels and unified legend ===
    xlabel('Log discharge (m^{3}/s)');
    ylabel(sprintf('Log dissolved %s (%s)', sFormattedInput, yUnit), 'Interpreter', 'tex');

    lgd = legend(legendHandles, legendLabels, 'Location', 'westoutside');
    title(lgd, 'Color = Site, Shape = Month');

    grid on;
    hold off;

    % === Save figure ===
    sInputChar = char(sInput);
    fileName = sprintf('conc_dis_log_usgs_power_unit_%s_Seasonality.jpg', sInputChar);
    fullFilePath = fullfile(folderName, fileName);
    saveas(gcf, fullFilePath, 'jpg');

    for j = 1:numel(siteData)
    regressionSummary(end+1, :) = {
        sInput, ...
        siteIDs{j}, ...
        a_vals(j), ...
        b_vals(j), ...
        r2_vals(j), ...
        eqn_strings{j}, ...
        log10(a_vals(j))  
    };
    end
end

%% Write regression results to Excel (separate sheets by site)
% Define headers for export
header = {'Variable', 'SiteID', 'a', 'b', 'R_squared', 'Equation', 'Intercept'};

% Combine header and data
outputFile = fullfile(folderName, 'regression_summary_all_sites.xlsx');
writecell([header; regressionSummary], outputFile);
