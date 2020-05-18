function addblob_contour_out(fh, varargin)

% Move contour object behind Surface object forward.
%
% :Usage:
% ::
%     addblob_contour_out(fh, varargin)
%
% ..
%     Author and copyright information:
%
%     Copyright (C) 2020  Jae-Joong Lee
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
%
% :Inputs:
%
%   **fh:**
%        Figure handle (usually, the current handle obtained from gcf)
%
% :Optional Inputs:
%
%
%   **'threshold':**
%        threshold for contour line. Large threshold makes sharp line,
%        and small threshold makes smooth and loose line.
%        (range: 0~1; default: 0.8)
%

contour_thresh = 0.8;
for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            case 'threshold'
                contour_thresh = varargin{i+1};
        end
    end
end

h = get(fh);
n_axes = numel(h.Children);
for ax_i = 1:n_axes
    cont_objs = findobj(h.Children(ax_i).Children, 'Type', 'Contour');
    surf_objs = findobj(h.Children(ax_i).Children, 'Type', 'Surface');
    for surf_i = 1:numel(surf_objs)
        surf_objs(surf_i).ZData = zeros(size(surf_objs(surf_i).ZData));
    end
    for cont_i = 1:numel(cont_objs)
        cont_max_Z = max(cont_objs(1).ZData(:));
        cont_objs(cont_i).LevelList = cont_max_Z * contour_thresh;
    end
end

end