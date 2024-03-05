function [Cp_max, Cp_min] = cookMayneGumbel(Cp, fsamp, T_target)
% Calculate Cppeak using Cook & Mayne method with Gumbel distrib.
% INPUTS:
% Cp        =  Cp time series
% fsamp     =  sampling freq of measurements
% T_target  =  target overall time series length (s)
% 
% OUTPUTS:
% Cp_max      =  Cook & Mayne method derived Cp_max
% Cp_min      =  Cook & Mayne method derived Cp_min

% If no T_target is given, assume T_target = T_actual
T_actual = length(Cp)/fsamp;
if ~exist('T_target','var')
    T_target = T_actual;
end

% 1. Divide time series into windows:
% window_edges = 1:round(window_dt*fsamp):length(Cp);
window_edges = round(linspace(1, length(Cp), 17));

% 2. Calculate the peaks in each window:
for i=1:length(window_edges)-1
    peaks_max(i) = max(Cp(window_edges(i):window_edges(i+1)), [], 'omitnan');
    peaks_min(i) = min(Cp(window_edges(i):window_edges(i+1)), [], 'omitnan');
end
peaks_max = peaks_max(~isnan(peaks_max));
peaks_min = peaks_min(~isnan(peaks_min));

% 3. Use the Gumbel equations (see Kasperski 2003)
alpha = T_target / T_actual;
alpha = max([1, alpha]); % can't be less than 1

% Max:
Vc = std(peaks_max) / abs(mean(peaks_max));
c_ad = 1 + 0.636 * Vc;
Cp_max = mean(peaks_max) * (c_ad + sqrt(6)/pi * log(alpha) * Vc);
% Min:
Vc = std(peaks_min) / abs(mean(peaks_min));
c_ad = 1 + 0.636 * Vc;
Cp_min = mean(peaks_min) * (c_ad + sqrt(6)/pi * log(alpha) * Vc);

