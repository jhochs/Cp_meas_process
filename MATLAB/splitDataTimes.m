function splitTimes = splitDataTimes(data, minTime)

% Split a struct containing data into discrete structs for each acquisition
% period

gapStartIdx = findMeasuringGaps(data);
if ~isempty(gapStartIdx)
    % First portion is up until the first measuring gap:
    splitTimes(1,:) = [data.t(1), data.t(gapStartIdx(1))];
    
    % If more than two portions, populate out the rest:
    i = 0; % if the following statement doesn't execute, this value is used on line 20
    if length(gapStartIdx)>1
        for i=1:length(gapStartIdx)-1
            splitTimes(i+1,:) = [data.t(gapStartIdx(i)+1), data.t(gapStartIdx(i+1))];
        end
    end
    
    % Last portion is after the last measuring gap:
    splitTimes(length(gapStartIdx)+1,:) = [data.t(gapStartIdx(i+1)+1), data.t(end)];
else
    splitTimes = [data.t(1), data.t(end)];
end

% Remove short periods:
splitTimes(diff(splitTimes,1,2) < seconds(minTime), :) = [];

end