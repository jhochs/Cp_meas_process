% This script is used for converting raw pressure measurements from motes to Cp statistics

clear

debug = false;
R = 287.05;

warning off MATLAB:table:RowsAddedExistingVars

%% Options
% % Data collection 1, Apr-May 2023
% motes = {'CM4', 'CM34', 'CM30', 'CM42', 'CM44', 'CM18', 'CM25', 'CM36', ...
%     'CM26', 'CM22', 'CM38', 'CM20', 'CM27', 'CM19', 'CM14', 'CM37', 'CM33', ...
%     'CM10', 'CM11', 'CM5', 'CM45', 'CM17', 'CM39', 'CM32', 'CM12'};
% data_dir = '../Data/650Cal_deployment_1_Apr-May';
% KOAK_filename = [data_dir, '/KOAK_ASOS.txt'];
% wind_filename = [data_dir, '/650CalRoof.mat']; % assuming CR300 used for data acquisition
% Cpstats_name = '650Cal_Cpstats_Gumbel_10min'; % name to save .mat and .csv

% Data collection 2, May-Jun 2023
motes = {'CM27', 'CM32', 'CM42', 'CM44', 'CM18', 'CM17', 'CM36', 'CM38', 'CM20'}; % excluded: 'CM4', 'CM11', 'CM22'
data_dir = '../Data/650Cal_deployment_2_May-Jun';
KOAK_filename = [data_dir, '/KOAK_ASOS.txt'];
wind_filename = [data_dir, '/650CalRoof.mat']; % assuming CR300 used for data acquisition
Cpstats_name = '650Cal_Cpstats_Gumbel_10min'; % name to save .mat and .csv

acqPeriodWindow = 10; %minutes
fsamp = 12.5; % measurement sampling frequency

%% Import raw data
data = loadMoteData(data_dir, motes);
if length(motes)==1
    data = {data};
end

load(wind_filename);  % imports as 'wind' struct

KOAK = readAWOSData(KOAK_filename);

%% Process wind data:
if ~isdatetime(wind.t)
    wind.t = datetime(wind.t*1000,'ConvertFrom','epochtime','TicksPerSecond',1e3,'Format','dd-MMM-yyyy HH:mm:ss.SSS');
end
if ~isfloat(wind.WDir)
    wind.WDir = double(wind.WDir);
end

WSavgperiod = 600;
wind.WSmovmean = movmean(wind.WS, WSavgperiod);

jMax = 1500;

%% Calculate stats
% Initialize table:
initLeng = 10000;
FS_meas = table();
FS_meas.Mote = cell([initLeng,1]);
FS_meas.AcqStart = NaT([initLeng,1]);
FS_meas.AcqLength = zeros([initLeng,1]);

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
FS_meas.dCpmean = NaN([initLeng,3]);
FS_meas.dCprms = NaN([initLeng,3]);
FS_meas.dCpmax = NaN([initLeng,3]);
FS_meas.dCpmin = NaN([initLeng,3]);
FS_meas.dCpmin_noEV = NaN([initLeng,3]); % just the min, no extreme value analysis
FS_meas.dCp_skewness = NaN([initLeng,3]);
FS_meas.dCp_kurtosis = NaN([initLeng,3]);

table_idx = 1;
for i=1:length(motes)
    fprintf('Processing mote %s:\n', motes{i});
    
    splitTimes = splitDataTimes(data{i}, 60*acqPeriodWindow);
    if ~isempty(splitTimes)
        windows = splitTimes(1,1):minutes(acqPeriodWindow):splitTimes(1,2);
        for j=2:size(splitTimes,1)
            windows = [windows, splitTimes(j,1):minutes(acqPeriodWindow):splitTimes(j,2)];
        end
        
        for j=1:length(windows)-1
            if minutes(windows(j+1)-windows(j)) > acqPeriodWindow
                fprintf('Skipping non-acquisition window: %s - %s...\n', datestr(windows(j), 'mmm dd hh:MM'), datestr(windows(j+1), 'mmm dd hh:MM'));
                continue
            end
            fprintf('Calculating Cp time series: %s - %s... \n', datestr(windows(j), 'mmm dd hh:MM'), datestr(windows(j+1), 'mmm dd hh:MM'));
            acqPortion = extractPortion(data{i}, windows(j), windows(j+1));
            windPortion = extractPortion(wind, windows(j), windows(j+1));
            
            % Check there is enough wind data and acq data in this time range:
            % Requires >70% of pressure data and >70% of wind data
            % Also verify mean windspeed > 5 m/s
            if length(windPortion.t) > 0.7*60*acqPeriodWindow && length(acqPortion.t) > 0.7*fsamp*60*acqPeriodWindow ...
                    && mean(windPortion.WS)>5
                [t, ~, dCp{i,j}] = CptimeSeries(acqPortion, KOAK, 'mote', KOAK, windPortion, false);
                
                [~, dCp_stats] = Cpstats(NaN, dCp{i,j}, [1 1 1], fsamp);
                
                % Record quantities in table:
                FS_meas.Mote{table_idx} = motes{i};
                FS_meas.AcqStart(table_idx) = windows(j);
                FS_meas.AcqLength(table_idx) = acqPeriodWindow;
                
                % Windspeed magnitude statistics:
                FS_meas.WSavg(table_idx) = mean(windPortion.WS);
                FS_meas.WSstd(table_idx) = std(windPortion.WS);
                FS_meas.TurbIntensity(table_idx) = FS_meas.WSstd(table_idx) / FS_meas.WSavg(table_idx);
                [FS_meas.WDiravg(table_idx), FS_meas.WDirstd(table_idx)] = windDirStats(windPortion.WDir);
                
                % Windspeed component statistics:
                u = cosd(FS_meas.WDiravg(table_idx) - windPortion.WDir) .* windPortion.WS;
                w = sind(FS_meas.WDiravg(table_idx) - windPortion.WDir) .* windPortion.WS;
                
                FS_meas.Uavg(table_idx) = mean(u);
                FS_meas.TurbIntensity_x(table_idx) = std(u) / mean(u);
                FS_meas.TurbIntensity_z(table_idx) = std(w) / mean(u);
                
                FS_meas.Tu(table_idx) = integralTimescaleSimple(u, 1);
                FS_meas.Tw(table_idx) = integralTimescaleSimple(w, 1);
                
                FS_meas.Lux(table_idx) = FS_meas.Uavg(table_idx) * FS_meas.Tu(table_idx);
                FS_meas.Lwx(table_idx) = FS_meas.Uavg(table_idx) * FS_meas.Tw(table_idx);
                
                % Cp statistics:
                FS_meas.dCpmean(table_idx,:) = dCp_stats(:,1)';
                FS_meas.dCprms(table_idx,:) = dCp_stats(:,2)';
                FS_meas.dCpmax(table_idx,:) = dCp_stats(:,3)';
                FS_meas.dCpmin(table_idx,:) = dCp_stats(:,4)';
                FS_meas.dCp_skewness(table_idx,:) = dCp_stats(:,5)';
                FS_meas.dCp_kurtosis(table_idx,:) = dCp_stats(:,6)';
                
                for k=1:3
                    if sum(isnan(dCp{i,j}(:,k))) < 0.2*size(dCp{i,j},1)
                        FS_meas.dCpmin_noEV(table_idx,k) = nanmin(dCp{i,j}(:,k));
                    end
                end
                
                table_idx = table_idx+1;
            end
        end
    end
end

%% Truncate and save table
FS_meas = FS_meas(~isnan(FS_meas.WSavg), :);

fprintf('Note: not saving resulting table\n');
return

save(sprintf('%s/Cpstats/%s.mat', data_dir, Cpstats_name), 'FS_meas');
writetable(FS_meas, sprintf('%s/Cpstats/%s.csv', data_dir, Cpstats_name));





