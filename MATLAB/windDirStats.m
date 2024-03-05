function [dir_mean, dir_std] = windDirStats(angles)
% Calculate mean and standard deviation of wind direction from a vector of
% angles in degrees

if isempty(angles)
    dir_mean = NaN;
    dir_std = NaN;
else
    angles = angles(:); % force into a column vector
    
    unit_vec_x = cosd(angles);
    unit_vec_y = sind(angles);
    xm = mean(unit_vec_x);
    ym = mean(unit_vec_y);
    
    dir_mean = wrapTo360(atan2d(ym, xm));
    
    from_mean = rad2deg(abs(angdiff(deg2rad(dir_mean*ones(size(angles))), deg2rad(angles))));
    dir_std = std(from_mean);
end
