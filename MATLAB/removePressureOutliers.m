function data = removePressureOutliers(data, avgTime, threshFac)

if ~exist('avgTime','var')
    % averaging time is not specified, default to 1 minute
    avgTime = minutes(1);
end
if ~exist('threshFac','var')
    % threshold is not specified, default to 20
    % this is a very high threshold since we don't want to unfairly exclude
    %  any peak measurements
    threshFac = 20;
end

if isfield(data, 'P')
    data.P = filloutliers(data.P, NaN, 'movmedian', avgTime, 'SamplePoints', data.t, 'ThresholdFactor', threshFac);
end
if isfield(data, 'Pa')
    data.Pa = filloutliers(data.Pa, NaN, 'movmedian', avgTime, 'SamplePoints', data.t, 'ThresholdFactor', threshFac);
end
if isfield(data, 'Pb')
    data.Pb = filloutliers(data.Pb, NaN, 'movmedian', avgTime, 'SamplePoints', data.t, 'ThresholdFactor', threshFac);
end
if isfield(data, 'Pc')
    data.Pc = filloutliers(data.Pc, NaN, 'movmedian', avgTime, 'SamplePoints', data.t, 'ThresholdFactor', threshFac);
end
