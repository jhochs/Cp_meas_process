function Tx = integralTimescaleSimple(x, fs, method)
% This function computes the integral timescale both by exponential
% fitting and integration of the curve. 
% 
% INPUTS:
% x          : QoI for which to calculate the timescale
% fs         : sampling frequency
% method     : 'integral' or 'fit'

% OUTPUTS:
% Tx         : mean timescale

if ~exist('method', 'var')
    % Default to integral method:
    method = 'integral';
end

intThresh = 0.02;
minLength = 100;

% Force column vector:
x = x(:);

% Fill missing values:
x = fillmissing(x, 'linear', 'EndValues', 'none');

% Truncate end values:
x(isnan(x)) = [];

N = length(x);

if N<minLength
%     error('Time series is too short to find integral time scale');
    Tx = NaN;
    return
end

[R, tau] = xcov(x, 'normalized');
R = R(N:end);
tau = tau(N:end)/fs;

if strcmp(method, 'fit')
    % Fit exponential to autocorrelation
    expFit = fittype('exp(-x/a)', 'dependent',{'y'}, 'independent',{'x'}, 'coefficients',{'a'});
    coeff = fit(tau', R, expFit, 'StartPoint', 1, 'Lower', 0);
    Tx = coeff.a;
elseif strcmp(method, 'integral')
    % Integrate autocorrelation
    zeroIdx = find(R - intThresh < 0);
    if isempty(zeroIdx)
        Tx = NaN;
    else
        zeroIdx = zeroIdx(1);
        Tx = trapz(tau(1:zeroIdx), R(1:zeroIdx));
    end
else
    error('Method ''%s'' not recognized, method should be either ''integral'' or ''fit''');
end

