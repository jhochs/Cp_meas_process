function station = readAWOSData(filename)
%Read RAW ASOS data, exported from https://mesonet.agron.iastate.edu/request/download.phtml?network=CA_ASOS
%This function assumes the following are selected to be exported:
% 1. Air Temperature [deg C]
% 2. Dew Point [deg C]
% 3. Wind Speed [knots]
% 4. Altimeter [inches]
% 5. Cloud Coverage Level 1
% 6. Cloud Coverage Level 2
% 7. Cloud Coverage Level 3

fid = fopen(filename);
data = textscan(fid, '%q%q%q%q%q%q%q%q%q', 'Delimiter', ',', 'HeaderLines', 1);
fclose(fid);
station.t = datetime(data{2});
station.T = data{3};
station.Td = data{4};
station.WS = data{5};
station.P = data{6};
for i=1:length(data{7})
    if strcmp(data{7}{i}, 'M')
        station.skyCover(i) = NaN;
    elseif strcmp(data{7}{i}, 'CLR')
        station.skyCover(i) = 0;
    elseif ismember('OVC', {data{7}{i}, data{8}{i}, data{9}{i}})
        station.skyCover(i) = 1;
    elseif ismember('BKN', {data{7}{i}, data{8}{i}, data{9}{i}})
        station.skyCover(i) = 6/8;
    elseif ismember('SCT', {data{7}{i}, data{8}{i}, data{9}{i}})
        station.skyCover(i) = 3.5/8;
    elseif ismember('FEW', {data{7}{i}, data{8}{i}, data{9}{i}})
        station.skyCover(i) = 1.5/8;
    else
        station.skyCover(i) = NaN;
    end
end

%Remove missing data:
station.WS(contains(station.WS, 'M')) = {'NaN'};
station.P(contains(station.P, 'M')) = {'NaN'};
station.T(contains(station.T, 'M')) = {'NaN'};
station.Td(contains(station.Td, 'M')) = {'NaN'};

station.WS = str2double(station.WS);
station.P = str2double(station.P);
station.T = str2double(station.T);
station.Td = str2double(station.Td);

%Convert units:
station.WS = station.WS*0.5144; % knots to m/s conversion
station.P = station.P*3386; % inHg to Pa conversion
station.T = station.T + 273.1; %degC to K conversion
station.Td = station.Td + 273.1; %degC to K conversion

%Fill missing T and P:
station.P = fillmissing(station.P, 'nearest');
station.T = fillmissing(station.T, 'nearest');
station.Td = fillmissing(station.Td, 'nearest');

%Remove duplicates:
[station.t, idx, ~] = unique(station.t);
station.P = station.P(idx);
station.T = station.T(idx);
station.Td = station.Td(idx);
station.WS = station.WS(idx);
station.skyCover = station.skyCover(idx);

% Make row vectors:
station.t = station.t(:)';
station.T = station.T(:)';
station.Td = station.Td(:)';
station.WS = station.WS(:)';
station.P = station.P(:)';
station.skyCover = station.skyCover(:)';

end