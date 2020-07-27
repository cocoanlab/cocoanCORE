function [h_line, h_patch] = wani_plot_shading(xaxis, mean, error, varargin)

% usage: [h_line, h_patch] = wani_plot_shading(xaxis, mean, error, varargin)
%
% input: xaxis, mean, error
%
% optional_inputs: 'color', 'color_shade', 'alpha', 'linewidth'

color = [0.3333, 0.6588, 1.0000]; % default color  

use_color_shade = false;
do_alpha = false;
alpha = 1;
linew = 2;
lines = '-';
upline = [];
lowline = [];

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            case {'color', 'color_line'} 
                color = varargin{i+1};
            case {'color_shade'} 
                use_color_shade = true;
                color2 = varargin{i+1};
            case {'alpha'} 
                do_alpha = true;
                alpha = varargin{i+1}; % degree-corrected SBM
            case {'linewidth'}
                linew = varargin{i+1};
            case {'linestyle'}
                lines = varargin{i+1};
            case {'upperline'}
                upline = varargin{i+1};
            case {'lowerline'}
                lowline = varargin{i+1};
        end
    end
end

if ~isempty(upline)
    upperline = upline; 
else
    upperline = mean + error;
end

if ~isempty(lowline)
    lowerline = lowline; 
else
    lowerline = mean - error;
end

xdata = [xaxis fliplr(xaxis) xaxis(1)];
ydata = [upperline fliplr(lowerline) upperline(1)];

if ~use_color_shade && ~do_alpha
    color2 = color +.3;
    color2(color2 > 1) = 1;
elseif ~use_color_shade && do_alpha
    color2 = color;
end

h_patch = patch(xdata,ydata,'y','linestyle', 'none', 'FaceColor', color2, 'faceAlpha', alpha);

hold on;

h_line = plot(xaxis, mean, lines, 'linewidth', linew, 'color', color, 'MarkerSize', 4, 'MarkerFaceColor', color);
% h_line = plot(xaxis, mean, '-', 'linewidth', linew, 'color', color);

end