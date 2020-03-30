function [out, boot_vals] = boot_mean_wani(vals, nboot, varargin)

% function [out, boot_vals] = boot_mean_wani(vals, nboot, varargin)
%
% boot for mean
%
% vals: values you want to bootstrap
% nboot: number of bootstraps (e.g., 10000)

% *optional input:
%
% 'noverbose'

doverbose = true;

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            % functional commands
            case {'noverbose'}
                doverbose = false;
        end
    end
end

boot_vals = bootstrp(nboot, @mean, vals);
out.bootmean = mean(boot_vals);
out.bootste = std(boot_vals);
out.bootZ = out.bootmean./out.bootste;
out.bootP = 2 * (1 - normcdf(abs(out.bootZ)));
out.ci95 = [prctile(boot_vals, 2.5); prctile(boot_vals, 97.5)];

if doverbose
    l = [mean(boot_vals); out.bootste; out.bootZ; out.bootP];
    fprintf('\nTest results: mean = %1.4f, sem = %1.4f, z = %1.4f, p = %1.5f', l);
end

end
