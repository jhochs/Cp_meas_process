function data = loadMoteData(dir, motes, calib_data)
%{
This function looks in dir for a .mat containing mote measurement data. To import 
the data from the AWS DynamoDB table or the .txt in the first  place, see 
PlotDynamoDBData.m (for cellular/wifi motes) or PlotSDData.m (for SD motes). 
INPUTS: 
dir        =  directory where data is stored
motes      =  a single mote ID or a cell array of mote IDs, e.g. {'cm22', 'wm1'}
calib_data =  (optional) whether should look for '_CALIB' suffix
%}

if ~exist('calib_data','var')
    % whether data is calib data is not specified, default to false
    calib_data = false;
end

if dir(end) == '/'
    dir = dir(1:end-1);
end

if ischar(motes)
    motes = {motes};
end

for i=1:length(motes)
    % Disconnected (SD) motes write calib data to a separate file:
    if calib_data
        temp = struct2cell(load(sprintf('%s/%s_CALIB.mat', dir, upper(motes{i}))));
    else % standard acquisition data
        try
            % Cellular/data mote format:
            temp = struct2cell(load(sprintf('%s/%s_data.mat', dir, lower(motes{i}))));
        catch
            % Disconnected (SD) mote format:
            temp = struct2cell(load(sprintf('%s/%s.mat', dir, upper(motes{i}))));
        end
    end
        
    data{i} = temp{1};
    if ~isdatetime(data{i}.t)
        data{i}.t = datetime(data{i}.t*1000,'ConvertFrom','epochtime','TicksPerSecond',1e3,'Format','dd-MMM-yyyy HH:mm:ss.SSS');
    end
    
    % Determine if reference data, which have only 2 pressure measurements:
    if ~isfield(data{i},'Pc') && isfield(data{i},'Pa')
        data{i}.Pavg = nanmean([data{i}.Pa; data{i}.Pb], 1);
        data{i}.Tavg = nanmean([data{i}.Ta; data{i}.Tb], 1);
    end
end

% If only one, don't need to store in cell array:
if length(data) == 1
    data = data{1};
end

end

