% This script is used for converting raw pressure measurements from motes to Cp statistics

% fprintf('Note: calculating TurbIntensity with WSmovmean subtracted (line 156)\n')

clear

debug = false;
R = 287.05;
D = 37.2; % SN diameter, for eta calculation
fsamp = 12.5; % measurement sampling frequency

warning off MATLAB:table:RowsAddedExistingVars

%% Options
acqPeriodWindow = 10; %minutes
WSavgperiod = 600; %600;

% Data collection 2, May-July 2022
% xls_filename = '../Data/Seattle_deployment_2_May-July/May-July_Sensor_acquisition_times.xlsx';
% data_dir = '../Data/Seattle_deployment_2_May-July';
% KBFI_filename = '../Data/Seattle_deployment_2_May-July/KBFI_ASOS.txt';
% wind_mote = 'CM32';
% Cpstats_name = 'SN_Cpstats_May-July_Gumbel_10min_60minmovmean'; % name to save .mat and .csv

% Data collection 3, Nov-Dec 2022
% xls_filename = '../Data/Seattle_deployment_3_Nov-Dec/Nov-Dec_Sensor_acquisition_times.xlsx';
% data_dir = '../Data/Seattle_deployment_3_Nov-Dec';
% KBFI_filename = '../Data/Seattle_deployment_3_Nov-Dec/KBFI_ASOS.txt';
% wind_mote = 'WM28';
% Cpstats_name = 'SN_Cpstats_Nov-Dec_Gumbel_10min_60minmovmean'; % name to save .mat and .csv

% Data collection 4, Dec 2022-Jan 2023
% xls_filename = '../Data/Seattle_deployment_4_Dec-Jan/Dec-Jan_Sensor_acquisition_times.xlsx';
% data_dir = '../Data/Seattle_deployment_4_Dec-Jan';
% KBFI_filename = '../Data/Seattle_deployment_4_Dec-Jan/KBFI_ASOS.txt';
% wind_mote = 'WM28';
% ref_mote = 'WM26';
% Cpstats_name = 'SN_Cpstats_Dec-Jan_Gumbel_10min_60minmovmean'; % name to save .mat and .csv

% Data collection 5, Feb 2023
xls_filename = '../Data/Seattle_deployment_5_Feb/Feb_Sensor_acquisition_times.xlsx';
data_dir = '../Data/Seattle_deployment_5_Feb';
KBFI_filename = '../Data/Seattle_deployment_5_Feb/KBFI_ASOS.txt';
wind_mote = 'WM28';
Cpstats_name = 'SN_Cpstats_Feb_Gumbel_10min'; % name to save .mat and .csv

%% Input is in XLS file
% Motes listed in `z` sheet should be in the same order as in `Acq` sheet
input = readtable(xls_filename, 'Sheet', 'Acq'); % 1st sheet has columns: Mote | Sensor A mask | B mask | C mask | Start time | End time
motes = unique(input.mote, 'stable');  % 'stable' ensures no reordering
input.mask = [input.mask_A, input.mask_B, input.mask_C];

z = xlsread(xls_filename, 'z'); % 2nd sheet has z, calculated using dP_ss_calc.m

%% Import raw data
data = loadMoteData(data_dir, motes);
if length(motes)==1
    data = {data};
end
for i=1:length(motes)
    data{i}.z = z(i,:);
end
windData = loadMoteData(data_dir, wind_mote);
if exist('ref_mote','var')
    refData = loadMoteData(data_dir, ref_mote);
end
KBFI = readAWOSData(KBFI_filename);

%% Moving mean on wind data:
windData.WSmovmean = movmean(windData.WS, WSavgperiod);

%% Plot
% motePlot(figure(1), data, windData, true, true);
% xlim([tAcq(1,1), tAcq(end,2)]);


%% Calculate stats
% Initialize table:
initLeng = 10000;
FS_meas = table();
FS_meas.Mote = cell([initLeng,1]);
FS_meas.AcqID = NaN([initLeng,1]);
FS_meas.AcqStart = NaT([initLeng,1]);
FS_meas.AcqLength = NaN([initLeng,1]);
FS_meas.SensorMask = NaN([initLeng,3]);

% Windspeed magnitude statistics
FS_meas.WSavg = NaN([initLeng,1]);
FS_meas.WSstd = NaN([initLeng,1]);
FS_meas.WDiravg = NaN([initLeng,1]);
FS_meas.WDirstd = NaN([initLeng,1]);
FS_meas.TurbIntensity = NaN([initLeng,1]);

% Windspeed component statistics
FS_meas.Uavg = NaN([initLeng,1]);
FS_meas.TurbIntensity_x = NaN([initLeng,1]);
FS_meas.TurbIntensity_z = NaN([initLeng,1]);
FS_meas.Tu = NaN([initLeng,1]);
FS_meas.Tw = NaN([initLeng,1]);
FS_meas.Lux = NaN([initLeng,1]);
FS_meas.Lwx = NaN([initLeng,1]);
FS_meas.eta = NaN([initLeng,1]);

% Cp statistics
FS_meas.Cpmean = NaN([initLeng,3]);
FS_meas.Cprms = NaN([initLeng,3]);
FS_meas.Cpmax = NaN([initLeng,3]);
FS_meas.Cpmin = NaN([initLeng,3]);
FS_meas.Cp_skewness = NaN([initLeng,3]);
FS_meas.Cp_kurtosis = NaN([initLeng,3]);
FS_meas.dCpmean = NaN([initLeng,3]);
FS_meas.dCprms = NaN([initLeng,3]);
FS_meas.dCpmax = NaN([initLeng,3]);
FS_meas.dCpmin = NaN([initLeng,3]);
FS_meas.dCpmin_noEV = NaN([initLeng,3]); % just the min, no extreme value analysis
FS_meas.dCp_skewness = NaN([initLeng,3]);
FS_meas.dCp_kurtosis = NaN([initLeng,3]);

FS_meas.Cpmean_std = NaN([initLeng,3]); % quantification of the variation of the mean within the sample
FS_meas.Cprms_std = NaN([initLeng,3]); % quantification of the variation of the rms within the sample

table_idx = 1;
clear Cp dCp
for i=1:size(input,1)
    fprintf('i = %d\n', i);
    mote_idx = find(contains(motes,input.mote(i)));
    if i==1 || ~strcmp(input.mote{i}, input.mote{i-1})
        acq_idx = 1;
    else
        acq_idx = acq_idx + 1;
    end
    
    windows = input.t_start(i):minutes(acqPeriodWindow):input.t_end(i);
    
    for j=1:length(windows)-1
        fprintf('Calculating Cp time series: %s-%s...\n', datestr(windows(j), 'hh:MM'), datestr(windows(j+1), 'hh:MM'));
        acqPortion = extractPortion(data{mote_idx}, windows(j), windows(j+1));
        windPortion = extractPortion(windData, windows(j), windows(j+1));
        
        % Check there is enough wind data and acq data in this time range:
        % Requires >70% of pressure data and >70% of wind data
        if length(windPortion.t) > 0.7*60*acqPeriodWindow && length(acqPortion.t) > 0.7*fsamp*60*acqPeriodWindow           
            if exist('ref_mote', 'var')
                [t, Cp{i,j}, dCp{i,j}] = CptimeSeries(acqPortion, refData, 'mote', KBFI, windPortion, false);
            else
                [t, Cp{i,j}, dCp{i,j}] = CptimeSeries(acqPortion, KBFI, 'mote', KBFI, windPortion, false);
            end
            
            % Record quantities in table:
            FS_meas.Mote{table_idx} = input.mote{i};
            FS_meas.AcqID(table_idx) = acq_idx;
            FS_meas.AcqStart(table_idx) = windows(j);
            FS_meas.AcqLength(table_idx) = acqPeriodWindow;
            FS_meas.SensorMask(table_idx,:) = input.mask(i,:);
            
            % Windspeed magnitude statistics:
            FS_meas.WSavg(table_idx) = mean(windPortion.WS);
            FS_meas.WSstd(table_idx) = std(windPortion.WS); % - windPortion.WSmovmean
            FS_meas.TurbIntensity(table_idx) = FS_meas.WSstd(table_idx) / FS_meas.WSavg(table_idx);
            [FS_meas.WDiravg(table_idx), FS_meas.WDirstd(table_idx)] = windDirStats(windPortion.WDir);
            
            % Windspeed component statistics:
            u{i,j} = cosd(FS_meas.WDiravg(table_idx) - windPortion.WDir) .* windPortion.WS;
            w{i,j} = sind(FS_meas.WDiravg(table_idx) - windPortion.WDir) .* windPortion.WS;
            
            FS_meas.Uavg(table_idx) = mean(u{i,j});
            FS_meas.TurbIntensity_x(table_idx) = std(u{i,j}) / mean(u{i,j});
            FS_meas.TurbIntensity_z(table_idx) = std(w{i,j}) / mean(u{i,j});
            
            FS_meas.Tu(table_idx) = integralTimescaleSimple(u{i,j}, 1);
            FS_meas.Tw(table_idx) = integralTimescaleSimple(w{i,j}, 1);
            
            FS_meas.Lux(table_idx) = FS_meas.Uavg(table_idx) * FS_meas.Tu(table_idx);
            FS_meas.Lwx(table_idx) = FS_meas.Uavg(table_idx) * FS_meas.Tw(table_idx);
            FS_meas.eta(table_idx) = FS_meas.TurbIntensity_x(table_idx) * (FS_meas.Lux(table_idx) / D).^0.15;
            
            % Cp statistics:
            [Cp_stats, dCp_stats] = Cpstats(Cp{i,j}, dCp{i,j}, FS_meas.SensorMask(table_idx,:), fsamp);
            
            FS_meas.Cpmean(table_idx,:) = Cp_stats(:,1)';
            FS_meas.Cprms(table_idx,:) = Cp_stats(:,2)';
            FS_meas.Cpmax(table_idx,:) = Cp_stats(:,3)';
            FS_meas.Cpmin(table_idx,:) = Cp_stats(:,4)';
            FS_meas.Cp_skewness(table_idx,:) = Cp_stats(:,5)';
            FS_meas.Cp_kurtosis(table_idx,:) = Cp_stats(:,6)';
            
            FS_meas.dCpmean(table_idx,:) = dCp_stats(:,1)';
            FS_meas.dCprms(table_idx,:) = dCp_stats(:,2)';
            FS_meas.dCpmax(table_idx,:) = dCp_stats(:,3)';
            FS_meas.dCpmin(table_idx,:) = dCp_stats(:,4)';
            FS_meas.dCp_skewness(table_idx,:) = dCp_stats(:,5)';
            FS_meas.dCp_kurtosis(table_idx,:) = dCp_stats(:,6)';
            
            for k=1:3
                if FS_meas.SensorMask(table_idx,k)
                    FS_meas.dCpmin_noEV(table_idx,k) = min(dCp{i,j}(:,k), [], 1, 'omitnan');
                else
                    FS_meas.dCpmin_noEV(table_idx,k) = NaN;
                end
            end
            
            % Break the window into 10 and calc mean and std for each, then
            % take overall std to quantify variability:
            sub_idx = round(linspace(1, length(Cp{i,j}), 11));
            cur_subwindows_means = zeros([length(sub_idx)-1, 3]);
            cur_subwindows_stds = zeros([length(sub_idx)-1, 3]);
            for l=1:(length(sub_idx)-1)
                cur_subwindows_means(l,:) = mean(Cp{i,j}(sub_idx(l):sub_idx(l+1), :), 1, 'omitnan');
                cur_subwindows_stds(l,:) = std(Cp{i,j}(sub_idx(l):sub_idx(l+1), :), [], 1, 'omitnan');
            end
            FS_meas.Cpmean_std(table_idx,:) = std(cur_subwindows_means, [], 1, 'omitnan');
            FS_meas.Cprms_std(table_idx,:) = std(cur_subwindows_stds, [], 1, 'omitnan');
            for k=1:3
                if ~FS_meas.SensorMask(table_idx,k)
                    FS_meas.Cpmean_std(table_idx,k) = NaN;
                    FS_meas.Cprms_std(table_idx,k) = NaN;
                end
            end

            table_idx = table_idx+1;
        end
    end
end

% Truncate and save:
FS_meas = FS_meas(~isnan(FS_meas.WSavg), :);

save(sprintf('%s/Cpstats/%s.mat', data_dir, Cpstats_name), 'FS_meas');
writetable(FS_meas, sprintf('%s/Cpstats/%s.csv', data_dir, Cpstats_name));

return


%% Plot in python - use instead LES_FS_comparison.py
% The following is obsolete:
%{
%% Plot
colors =   [0 0 1; ...
            0 0.5 1; ...
            0.1 0.9 0.7; ...
            0.5 0.9 0; ...
            0.9 0.9 0.1; ...
            1 0.7 0.1; ...
            1 0.5 0.5; ...
            1 0 0; ...
            0.7 0 0.9];
colors_sensor = [0 0.5 0; 1 0 0; 0 0 1];
markers = {'+', 'd', 'o', 'x'};
labels = {'C_{p, mean}', 'C_{p, max}', 'C_{p, min}', 'C_{p, rms}'};
yRanges = [-1.5 0.5; ...
           -2 2; ... % possibly change
           -2 1; ...
           0 0.3];

if plot_dCp
    plot_idx = 3:4; % plot min and rms only
else
    plot_idx = [1, 3, 4]; % plot mean, min, rms
end

for i=1:length(motes)
    f = figure(i); clf
    if plot_dCp
        f.Position = [10+300*(i-1) 100 400 450];
    else
        f.Position = [10+300*(i-1) 100 400 600];
    end
end

% Loop through mote / acquisition combinations (i), statistics (l)
for i=1:size(input,1)
    mote_idx = find(contains(motes,input.mote(i)));
    figure(mote_idx)
    if i==1 || ~strcmp(input.mote{i}, input.mote{i-1})
        color_idx = 1;
    else
        color_idx = color_idx + 1;
    end

    for l=1:length(plot_idx)
        subplot(length(plot_idx),1,l); hold on

        if strcmp(plot_by, 'sensor')
            % Additional loop over sensors (k)
            for k=1:3
                if plot_dCp
                    y = dCp_stats(i,:,k,plot_idx(l));
                else
                    y = Cp_stats(i,:,k,plot_idx(l));
                end
                x = WDiravg(i,:);

                plot(x(:), y(:), markers{plot_idx(l)}, 'color', colors_sensor(k,:), 'displayname', sprintf('Sensor %d', k))
            end
        else % either turb_intensity of acq_day
            if plot_dCp
                y = dCp_stats(i,:,(input.mask(i,:)==1),plot_idx(l));
            else
                y = Cp_stats(i,:,(input.mask(i,:)==1),plot_idx(l));
            end
            x = repmat(WDiravg(i,:), 1, size(y,3));

            if strcmp(plot_by, 'acq_day')
                if isnat(input.t_start(i))
                    seriesName = '';
                else
                    seriesName = datestr(input.t_start(i), 'dd mmm');
                end
                plot(x(:), y(:), markers{plot_idx(l)}, 'color', colors(color_idx,:), 'displayname', seriesName)
            else
                c = repmat(TI(i,:), 1, size(y,3));
                scatter(x(:), y(:), 20, c(:), markers{plot_idx(l)})
                caxis([0.05 0.15])
                colormap('jet');
            end
        end

        if plot_dCp
            ylabel(['\Delta', labels{plot_idx(l)}])
            yRanges(3,:) = [-1.2 0.05];
        else
            ylabel(labels{plot_idx(l)})
        end

        xlim(xRange)
        ylim(yRanges(plot_idx(l),:));

        if l==length(plot_idx)
            xlabel('\theta [�]')
        end
        set(gca, 'fontsize', 13)
    end
end
%}