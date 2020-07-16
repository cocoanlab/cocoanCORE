function out = plot_specificity_box(y1, y2, varargin)

% This function draws two violin plots and multiple lines for comparing two
% groups of same sample size. Each line connects two dependent samples, and
% color of the line indicates which sample is greater than the other.
% 
%
% :Usage:
% ::
%     plot_specificity_box(y1, y2, varargin)
%
%
% :Input:
% ::
%   - y1                 The first group to be compared
%   - y2                 The second group to be compared
%                        (N.B. y1 and y2 have same dimensions)
%
%
% :Optional Input:
% :: 
%   - 'colors', 'color'  [2 x 3] matrix for indiciating two colors of
%                        violin plot y1 (first row) and y2 (second row).
%                        (default: first row is red, second row is yellow)
%   - 'linecolors',      [2 x 3] matrix for indiciating colors of lines
%      or 'linecolor'    connecting two violin plots.
%                        First row indicates line color for y1 > y2,
%                        and the second row indicates the line color for
%                        y2 > y1.
%                        (default: first row is red, second row is blue)
%
%
% :Output:
% ::   
%
%
% :Example:
% ::
%   y1 = randi(100, 20, 1);
%   y2 = randi(100, 20, 1);
%   out = plot_specificity_box(y1, y2);
%
%   
% ..
%     Author and copyright information:
%
%     Copyright (C) Mar 2020  Choong-Wan Woo & Jae-Joong Lee
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
%
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
%
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
% ..


cols = [0.8353    0.2431    0.3098
    0.9922    0.6824    0.3804];
lcols = [0.8902    0.1020    0.1098
    0.1961    0.5333    0.7412];

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            % functional commands
            case {'colors', 'color'}
                cols = varargin{i+1};
            case {'linecolors', 'linecolor'}
                lcols = varargin{i+1};
        end
    end
end

data = [y1 y2];
boxplot_wani_2016(data, 'color', cols, 'linewidth', 2, 'boxlinewidth', 1, 'mediancolor', 'k', 'violin');

xdot{1} = ones(size(data,1),1)*1+.32;
xdot{2} = ones(size(data,1),1)*2-.32;

wh = data(:,1) > data(:,2);

out.h1 = line([xdot{1}(wh) xdot{2}(wh)]', data(wh,1:2)', 'color', lcols(1,:), 'linewidth', 1.5);
out.h2 = line([xdot{1}(~wh) xdot{2}(~wh)]', data(~wh,1:2)', 'color', lcols(2,:), 'linewidth', 1.5);

scatter(xdot{1}, data(:,1), 20, cols(1,:), 'filled'); 
scatter(xdot{2}, data(:,2), 20, cols(2,:), 'filled'); 

set(gcf, 'position', [1   518   204   187]);
set(gca, 'fontsize', 18, 'xlim', [0.5 2.5], 'linewidth', 1.5, 'ticklength', [.03 .03], 'xticklabel', '');

end