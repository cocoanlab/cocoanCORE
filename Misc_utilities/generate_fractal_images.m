function generate_fractal_images(varargin)

% This function generates fractal images.
%
% Reference: "Generation of fractal patterns for probing the visual memory"
%            Miyashita et al., 1992, Neuroscience Research
%
%
% :Usage:
% ::
%     generate_fractal_images
%
%
% :Input:
% ::
%
%
% :Optional Input:
% ::
%   - range_suppos       range of number of superposition. ex) [1 4], ...
%                        (default: [2 4])
%   - range_initvert     range of initial number of edges. ex) [1 4], ...
%                        (default: [3 5])
%   - range_recur        range of depth of recursion. ex) [1 4], ...
%                        (default: [3 5])
%   - range_cols         range of colormap.
%                        (default: 12 qualitative colors from Colorbrewer2)
%   - def_coef           coefficient of edge deflection.
%                        extent of deflection is defined as
%                        [randn * def_coef(1) + def_coef(2)].
%                        (default: [1 0])
%   - back_cols          background color.
%                        (default: [50 50 50] ./ 255)
%
%
% :Output:
% ::
%
%
% :Example:
% ::
%   % generate fractal images.
%   generate_fractal_images('range_suppos', [3 3], 'range_initvert', [5 10], 'range_recur', [5 7], 'def_coef', [0.2 0.1])
%
%
% ..
%     Author and copyright information:
%
%     Copyright (C) Jun 2020  Jae-Joong Lee
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

range_suppos = [2 4];
range_initvert = [3 5];
range_recur = [3 5];
range_cols = [166,206,227
    31,120,180
    178,223,138
    51,160,44
    251,154,153
    227,26,28
    253,191,111
    255,127,0
    202,178,214
    106,61,154
    255,255,153
    177,89,40] ./ 255; % 12 qualitative colors from Colorbrewer2
def_coef = [1 0];
back_cols = [50 50 50] ./ 255;

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            case {'range_suppos'}
                range_suppos = varargin{i+1};
            case {'range_initvert'}
                range_initvert = varargin{i+1};
            case {'range_recur'}
                range_recur = varargin{i+1};
            case {'range_cols'}
                range_cols = varargin{i+1};
            case {'def_coef'}
                def_coef = varargin{i+1};
            case {'back_cols'}
                back_cols = varargin{i+1};
        end
    end
end


n_suppos = randi(range_suppos); % number of superposition

figure;

for suppos_i = 1:n_suppos
    
    n_initvert = randi(range_initvert); % number of initial edges
    n_recur = randi(range_recur); % number of recursion
    GA = randn * def_coef(1) + def_coef(2);
    cols = range_cols(randi(size(range_cols,1)), :);
    
    pgon = nsidedpoly(n_initvert, 'Center', [0 0], 'Radius', 1);
    pgon = patch(pgon.Vertices(:,1), pgon.Vertices(:,2), cols);
    pgon.EdgeAlpha = 0;
    
    for recur_i = 1:n_recur
        vert_list = pgon.Vertices;
        delete(pgon);
        def_vert_list = deflect(vert_list, GA);
        pgon = patch(def_vert_list(:,1), def_vert_list(:,2), cols);
        pgon.EdgeAlpha = 0;
    end
    
end

xymax = max(abs([get(gca, 'xlim') get(gca, 'ylim')]));
set(gca, 'xlim', [-xymax xymax], 'ylim', [-xymax xymax]);
set(gcf, 'name', 'Fractal', 'color', 'white', 'position', [669   368   604   571]);
rectangle('Position', [-xymax -xymax 2*xymax 2*xymax], 'FaceColor', back_cols, 'EdgeColor', back_cols)
h = get(gca, 'Children');
set(gca, 'Children', h([2:end, 1]));
axis off;

end

function def_vert_list = deflect(vert_list, GA)

n_vert = size(vert_list,1);

def_vert_list = zeros(n_vert * 2, 2);
def_vert_list(2:2:end,:) = vert_list;

mx = (vert_list([1:n_vert], 1) + vert_list([n_vert 1:n_vert-1], 1)) / 2;
my = (vert_list([1:n_vert], 2) + vert_list([n_vert 1:n_vert-1], 2)) / 2;
dx = vert_list([1:n_vert], 1) - vert_list([n_vert 1:n_vert-1], 1);
dy = vert_list([1:n_vert], 2) - vert_list([n_vert 1:n_vert-1], 2);
theta = atan2(dy, dx);
def_vert_list(1:2:end, 1) = mx + GA .* sin(theta);
def_vert_list(1:2:end, 2) = my - GA .* cos(theta);

end
