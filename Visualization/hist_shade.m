function handle = hist_shade(w, varargin)

% Draw a histogram with shades of the upper and lower percentile areas.
%
% Usage:
% -------------------------------------------------------------------------
% h = hist_shade(w, varargin)
%
% Inputs:
% -------------------------------------------------------------------------
% w              the target vector (suitable for bootstrapped or permutated
%                values
%
% Optional inputs: Enter keyword followed by variable with values
% 'color'        'color', [.33 .66 1] (default: [.7 .7 .7]
% 'bins'         the number of bins (default: 100)
% 'shade_range'  lower and upper percentile, e.g., 'shade_range', [0 95]
%                default: [2.5 97.5]
%
% Example:
% w = normrnd(100,5,10000,1);
% h = hist_shade(w, 'color', [.33 .66 1], 'bins', 50, 'shade_range', [0 95]);

handle.gcf = create_figure('hist_shade');

set(gcf, 'position', [16   425   436   280]);

colors = [.7 .7 .7];
bins = 100;
lpct = 2.5;
hpct = 97.5;

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            % functional commands
            case {'color', 'colors'}
                colors = varargin{i+1}; 
            case {'bins'}
                bins = varargin{i+1};
            case {'shade_range'}
                lpct = varargin{i+1}(1);
                hpct = varargin{i+1}(2);
        end
    end
end

c2 = colors -.3;
c2(c2<0) = 0;

[handle.h, xx] = hist(w, bins);
    
handle.hbar = bar(xx, handle.h);
set(handle.hbar, 'FaceColor', colors, 'EdgeColor', colors);
    
wh1 = xx > prctile(w, hpct);
wh2 = xx < prctile(w, lpct);

h = handle.h;
h(~(wh1 | wh2)) = 0;
    
handle.hbar2 = bar(xx, h);
set(handle.hbar2, 'FaceColor', c2, 'EdgeColor',c2);
    
set(gca, 'fontsize', 20, 'linewidth', 2, 'tickdir', 'out', 'ticklength', [.02 .02]);

end