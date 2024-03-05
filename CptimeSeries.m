function [t, Cp, dCp] = CptimeSeries(Pdata, PrefData, TrefData, TdData, windData, constWS)
%{
This function calculates Cp statistics from pressure time series
INPUTS:
Pdata       =   struct that contains the measurements of pressure on building, 
                must also contain z data (zero-wind offset) if Cp measurements 
                are desired. Not required for Cp' (dCp) measurements
PrefData    =   struct that contains barometric reference pressure
TrefData    =   either the struct that contains reference temperature data, or 
                'mote' to just use the temperatures contained in Pdata
TdData      =   struct that contains dewpoint data
windData    =   struct that contains wind measurements
constWS     =   true to use a constant windspeed, false to use the moving mean
OUTPUTS:
t           =   time
Cp          =   Cp time series
dCp         =   Cp' (mean-subtracted Cp) time series
%}

g = 9.81;
R = 287.05;
fs = 12.5; % sensor sampling frequency
dCp_movmean_period = 600 * fs; % period of moving mean to subtract to get Delta(Cp)

% Extract pressure data:
Pdata = removePressureOutliers(Pdata);
t = Pdata.t;

% Movmean, and interpolate wind data:
if constWS
    WSinterp = mean(windData.WS);
else
    WSinterp = interp1(windData.t, windData.WSmovmean, t);
end

% Interpolate reference pressure data:
if isfield(PrefData, 'Pavg')
    % Data source is a reference mote:
    Prefinterp = interp1(PrefData.t, PrefData.Pavg, t);
else
    % Data source is AWOS station:
    if length(PrefData.t)==1
        % only one measurement, can't interpolate
        Prefinterp = PrefData.P * ones(size(t));
    else
        Prefinterp = interp1(PrefData.t, PrefData.P, t);
    end
end

% Interpolate reference temperature data:
if ischar(TrefData) && strcmp(TrefData, 'mote')
    Trefinterp = [Pdata.Ta; Pdata.Tb; Pdata.Tc];
else
    if length(TrefData.t)==1
        % only one measurement, can't interpolate
        Trefinterp = TrefData.T * ones([3,length(t)]);
    else
        Trefinterp = repmat(interp1(TrefData.t, TrefData.T, t), [3,1]);
    end
end

% Interpolate dewpoint data:
if length(TdData.t)==1
    % only one measurement, can't interpolate
    Tdinterp = TdData.Td * ones(size(t));
else
    Tdinterp = interp1(TdData.t, TdData.Td, t);
end

% Convert to kelvin if in celsius
if mean(Trefinterp, 'all')<50
    Trefinterp = Trefinterp + 273.15;
end
if mean(Tdinterp, 'all')<50
    Tdinterp = Tdinterp + 273.15;
end

% Calculate density and dynamic pressure at every t
thetav = virtualPotentialTemp(Trefinterp, Tdinterp, Prefinterp);
rho = Prefinterp ./ (R * thetav);
q = 0.5 * rho .* WSinterp.^2; % will be more or less const. if constWS = true

dCp(:,1) = (Pdata.Pa - movmean(Pdata.Pa, dCp_movmean_period, 'omitnan')) ./ q(1,:);
dCp(:,2) = (Pdata.Pb - movmean(Pdata.Pb, dCp_movmean_period, 'omitnan')) ./ q(2,:);
dCp(:,3) = (Pdata.Pc - movmean(Pdata.Pc, dCp_movmean_period, 'omitnan')) ./ q(3,:);

if isfield(Pdata, 'z')
    for k=1:3
        dP_ss(k) = mean(Prefinterp .* (1 - exp(-g*Pdata.z(k)./(R*thetav(k,:)))), 'all', 'omitnan');
    end
    Cp(:,1) = (Pdata.Pa - Prefinterp + dP_ss(1)) ./ q(1,:);
    Cp(:,2) = (Pdata.Pb - Prefinterp + dP_ss(2)) ./ q(2,:);
    Cp(:,3) = (Pdata.Pc - Prefinterp + dP_ss(3)) ./ q(3,:);
else
    Cp = NaN * ones(size(dCp));
end




