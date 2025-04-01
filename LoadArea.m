%% Plots
clear all, close all, clc %use USGSGaugeDataset
[sFile, sPath] = uigetfile('*.xlsx', 'Select Database File');
sFullPath = fullfile(sPath, sFile);

tData = readtable(sFullPath);
%% 
vColumnLabels = tData.Properties.VariableNames(13:end);
vColumnLabelsArray = cellstr(vColumnLabels);
mFullData = table2array(tData(1:end, 13:end)); %first value is the number of samples (each as a row);second value is the numeric values that will be plotted, must be in an order that can be iterated through (DONT INCLUDE STRINGS)

vSampleDates = datetime(convertStringsToChars(string(table2cell(tData(:, 10))))); % just date not datetime
vSampleLocations = string(table2cell(tData(:, 3))); %labels of each sample (KR4, MR3, etc)
vWatershed = string(table2cell(tData(:,4))); % column that could be used to discriminate between watersheds each location is contained in 
vSampleYear = table2array(tData(:,12)); %column for year
vDischarge = table2array(tData(:,13));
%%
% Knik River 3 elemental data and sample dates
mKR3Data = mFullData(vSampleLocations == "KR3", :);
mKR3Dates = vSampleDates(vSampleLocations == "KR3");
mKR3Year = vSampleYear(vSampleLocations == "KR3");
mKR3_2022 = mKR3Data(mKR3Year == 2022, :);
mKR3_2023 = mKR3Data(mKR3Year == 2023, :);
mKR3_Date22 = mKR3Dates(mKR3Year == 2022);
mKR3_Date23 = mKR3Dates(mKR3Year == 2023);

% Matanuska River 4 elemental data and sample dates
mMR4Data = mFullData(vSampleLocations == "MR4", :);
mMR4Dates = vSampleDates(vSampleLocations == "MR4");
mMR4Year = vSampleYear(vSampleLocations == "MR4");
mMR4_2022 = mMR4Data(mMR4Year == 2022, :);
mMR4_2023 = mMR4Data(mMR4Year == 2023, :);
mMR4_Date22 = mMR4Dates(mMR4Year == 2022);
mMR4_Date23 = mMR4Dates(mMR4Year == 2023);

% Little Susitna River 2 elemental data and samples dates
mLS2Data = mFullData(vSampleLocations == "LS2", :);
mLS2Dates = vSampleDates(vSampleLocations == "LS2");
mLS2Year = vSampleYear(vSampleLocations == "LS2");
mLS2_2022 = mLS2Data(mLS2Year == 2022, :);
mLS2_2023 = mLS2Data(mLS2Year == 2023, :);
mLS2_Date22 = mLS2Dates(mLS2Year == 2022);
mLS2_Date23 = mLS2Dates(mLS2Year == 2023);

% Moose Creek elemental data and samples dates
mMCData = mFullData(vSampleLocations == "MC", :);
mMCDates = vSampleDates(vSampleLocations == "MC");
mMCYear = vSampleYear(vSampleLocations == "MC");
mMC_2022 = mMCData(mMCYear == 2022, :);
mMC_2023 = mMCData(mMCYear == 2023, :);
mMC_Date22 = mMCDates(mMCYear == 2022);
mMC_Date23 = mMCDates(mMCYear == 2023);

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
    title(sFormattedInput);
    legend('Location', 'none');
    legend('Position', [0.48 0.5 0.1 0.1], 'Location', 'none');
    grid on;
    xlim([datetime(2022,5,1), datetime(2023,10,03)]);
    hold off;
    sInputChar = char(sInput);
    formatFileName = ['load_area_by_time', sInputChar, 'jpg'];
    fileName = sprintf(formatFileName);
    fullFilePath = fullfile(folderName, fileName);
    saveas(gcf, fullFilePath, 'jpg');
end

%% save data to excel sheet
vSampleID = tData.SampleID;    
vSampleDate = tData.SampleDate;  
vBYULabID = tData.BYULabID;             

% Watershed areas (in km^2)
area_KR3 = 3156.7014;
area_MR4 = 5329.1646;
area_MC = 125.4186;
area_LS2 = 160.4475;

% Preallocate matrix for load/area data
[m, n] = size(mFullData);
mLoadAreaData = nan(m, n);

% Calculate load/area for each row
for i = 1:m
    loc = vSampleLocations(i);
    discharge = mFullData(i, 1); % assuming column 1 = discharge
    if isnan(discharge)
        continue
    end

    switch loc
        case "KR3"
            area = area_KR3;
        case "MR4"
            area = area_MR4;
        case "MC"
            area = area_MC;
        case "LS2"
            area = area_LS2;
        otherwise
            continue
    end

    mLoadAreaData(i, :) = mFullData(i, :) * 86.4 * discharge / area;
end

% Create final table
T_out = table(vSampleID, vSampleDate, vBYULabID, 'VariableNames', {'SampleID', 'SampleDate', 'BYULabID'});
T_values = array2table(mLoadAreaData, 'VariableNames', vColumnLabels);
T_final = [T_out T_values];

% Save to Excel
writetable(T_final, 'USGS_Load_Area_Data.xlsx');
disp('✅ Load/area values saved to USGS_Load_Area_Data.xlsx');
%%
% Separate discharge and element data
vDischargeOnly = mFullData(:,1);          % First column = discharge
mConcentrationOnly = mFullData(:,2:end);  % Remaining columns = concentrations

% Preallocate matrix for load/area results (same size as concentration-only part)
[m, n] = size(mConcentrationOnly);
mLoadAreaData = nan(m, n);

% Calculate load/area only for elemental concentrations
for i = 1:m
    loc = vSampleLocations(i);
    discharge = vDischargeOnly(i);
    if isnan(discharge)
        continue
    end

    switch loc
        case "KR3"
            area = area_KR3;
        case "MR4"
            area = area_MR4;
        case "MC"
            area = area_MC;
        case "LS2"
            area = area_LS2;
        otherwise
            continue
    end

    mLoadAreaData(i, :) = mConcentrationOnly(i, :) * 86.4 * discharge / area;
end

% Reattach the discharge column if you want to keep it in final table
T_out = table(vSampleID, vSampleDate, vBYULabID, vDischargeOnly, 'VariableNames', ...
              {'SampleID', 'SampleDate', 'BYULabID', 'Discharge_m3_s'});

T_values = array2table(mLoadAreaData, 'VariableNames', vColumnLabels(2:end)); % exclude discharge label
T_final = [T_out T_values];

% Save to Excel
writetable(T_final, 'USGS_Load_Area_Data_new.xlsx');
