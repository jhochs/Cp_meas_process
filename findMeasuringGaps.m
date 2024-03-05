function gapStartIdx = findMeasuringGaps(mote)

diffs = diff(mote.t);
mask = diffs > minutes(19);
gapStartIdx = find(mask);
% mask(end+1) = false;
% gapStarts = mote.t(mask); % must pad with false because diffs has dim N-1

end