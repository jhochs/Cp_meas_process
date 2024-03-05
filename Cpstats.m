function [Cp_stats, dCp_stats] = Cpstats(Cp, dCp, mask, CM_fsamp)

% Output columns are Cp_:
% 1. Mean
% 2. RMS
% 3. Max
% 4. Min
% 5. Skew
% 6. Kurt

Cp_stats = NaN([3,6]);
dCp_stats = NaN([3,6]);

for k=1:3
    if mask(k) % if sensor mask is 0, skip this one
        if sum(isnan(Cp)) ~= numel(Cp)
            Cp_stats(k,1) = nanmean(Cp(:,k));
            Cp_stats(k,2) = nanstd(Cp(:,k));
            [Cp_stats(k,3), Cp_stats(k,4)] = cookMayneGumbel(Cp(:,k), CM_fsamp, 160*60);
            Cp_stats(k,5) = skewness(Cp(:,k));
            Cp_stats(k,6) = kurtosis(Cp(:,k));
        end
        
        dCp_stats(k,1) = 0;
        dCp_stats(k,2) = nanstd(dCp(:,k));
        [dCp_stats(k,3), dCp_stats(k,4)] = cookMayneGumbel(dCp(:,k), CM_fsamp, 160*60);
        dCp_stats(k,5) = skewness(dCp(:,k));
        dCp_stats(k,6) = kurtosis(dCp(:,k));
    end
end