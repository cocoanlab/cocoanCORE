function generate_fractal_images

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
%   - TBD                TBD
%
%
% :Output:
% ::
%
%
% :Example:
% ::
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

n_suppos = randi([2 4]); % number of superposition
cols_all = jet(128);
cols_all = [166,206,227
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
    177,89,40] ./ 255;

figure;

for suppos_i = 1:n_suppos
    
    n_edge = randi([3 6]); % number of edges
    n_recur = randi([3 5]); % number of recursion
    GA = randn;
    cols = cols_all(randi(size(cols_all,1)), :);
    
    pgon = nsidedpoly(n_edge, 'Center', [0 0], 'Radius', 1);
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
