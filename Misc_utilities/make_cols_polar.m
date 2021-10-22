function out = make_cols_polar(cols, ratio, varargin)

% This function polarize colormap scheme based on the specified center and
% polarization ratio.
%
%
% :Usage:
% ::
%     make_cols_polar(cols, wh_center, ratio)
%
%
% :Input:
% ::
%   - cols               color scheme. (N x 3)
%   - ratio              polarization ratio.
%                        e.g., correlation matrix A ranges from [-0.4 0.8].
%                              -> to re-balance colormap by zero-centering,
%                                 ratio would be [0.4 0.8] or [1/3 2/3];
%
%
% :Optional Input:
% ::
%   - wh_center          index of color to be used as the center.
%                        (default is center; N/2 for even, (N-1)/2 for odd)
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
%     Copyright (C) Oct 2021  Jae-Joong Lee
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

n = size(cols, 1);
wh_center = floor(n/2);

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            case {'wh_center'}
                wh_center = varargin{i+1};
        end
    end
end

orig_ratio = [wh_center n-wh_center] ./ n;
ratio = ratio ./ sum(ratio);
if ratio(1) > orig_ratio(1) % lower part should be inflated
    lower_factor = (ratio(1) / ratio(2)) / (orig_ratio(1) / orig_ratio(2));
    col_idx = [round(linspace(1,wh_center,wh_center*lower_factor)) wh_center+1:n];
elseif ratio(2) > orig_ratio(2) % upper part should be inflated
    upper_factor = (ratio(2) / ratio(1)) / (orig_ratio(2) / orig_ratio(1));
    col_idx = [1:wh_center round(linspace(wh_center+1,n,wh_center*upper_factor))];
end

out = cols(col_idx, :);


end