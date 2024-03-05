function portion = extractPortion(varargin)
%{
For all measurement data stored in structs with field 't' as the time measurement,
this function can be used to extract a subset of that struct's measurements.
INPUTS:
1 (data)        =   struct of measurement data
2 (t0)          =   start time with which to filter
3 (t1)          =   end time with which to filter
4 (ornearest)   =   (optional) 'ornearest' flag, specifies to use the nearest 
                    measurement if none lie between the specified range. Useful
                    for AWOS measurements, for example, which only record one
                    measurement every 5 mins.
%}

if nargin<3
    error('Input 1: data struct, input 2: t0, input 3: t1');
end

data = varargin{1};
t0 = varargin{2};
t1 = varargin{3};

if ~isdatetime(data.t)
    data.t = datetime(data.t*1000,'ConvertFrom','epochtime','TicksPerSecond',1e3,'Format','dd-MMM-yyyy HH:mm:ss.SSS');
end
if ischar(t0)
    t0 = datetime(t0);
end
if ischar(t1)
    t1 = datetime(t1);
end

mask = isbetween(data.t, t0, t1);
if sum(mask)==0 && nargin>3 && strcmp(varargin{4}, 'ornearest')
    [~, mask] = min(abs(data.t - t0));
end

flds = fieldnames(data);
z_mask = strcmp(flds,'z');
if sum(z_mask)==1
    portion.z = data.z;
    flds(z_mask) = []; % not necessary for the field 'z' since this is const.
end
for i=1:length(flds)
    if ~strcmpi('diffs', flds{i})
        % check whether is a column or row array:
        if size(data.(flds{i}), 1)>size(data.(flds{i}), 2)
            portion.(flds{i}) = data.(flds{i})(mask,:);
        else
            portion.(flds{i}) = data.(flds{i})(:,mask);
        end
    end
end
