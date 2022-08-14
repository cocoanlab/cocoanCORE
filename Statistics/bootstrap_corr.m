function [bootP boot_vals rngset] = bootstrap_corr(orig_corr, bootnum, varargin)

no_plot = false;
rngset = rng;

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            case {'only_p'}
                no_plot = true;
            case {'rng'}
                rng(varargin{i+1});
                rngset = rng;
        end
    end
end

orig_corr_z = reformat_r_new(orig_corr, 'r2z');
boot_vals = bootstrp(bootnum, @mean, orig_corr_z);
bootmean = mean(boot_vals);
bootste = std(boot_vals);
bootZ = bootmean./bootste;
bootP = 2 * (1 - normcdf(abs(bootZ)));
fprintf('Niter: %d.  Corr: %.4f\nPermutation - Mean: %.4f,  SD: %.4f,  P-value: %d\n', ...
    bootnum, mean(orig_corr), reformat_r_new(bootmean, 'z2r'), reformat_r_new(bootste, 'z2r'), bootP);

if ~no_plot
    bootci95 = [prctile(boot_vals, 2.5); prctile(boot_vals, 97.5)];
    z2r_bootci95 = reformat_r_new(bootci95, 'z2r');
    z2r_boot_vals = reformat_r_new(boot_vals, 'z2r');
    histogram(z2r_boot_vals);
    line([z2r_bootci95(1) z2r_bootci95(1)], get(gca, 'ylim'), 'color', 'r');
    line([z2r_bootci95(2) z2r_bootci95(2)], get(gca, 'ylim'), 'color', 'r');
end

end