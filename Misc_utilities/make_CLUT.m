function make_CLUT(cols, clutpath)

% This function makes MRICROGL color lookup table (CLUT) file.
%
%
% :Usage:
% ::
%     make_CLUT(cols, clutpath)
%
%
% :Input:
% ::
%   - cols               color scheme. (N x 3)
%   - clutpath           path to save CLUT file
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
%     Copyright (C) Oct 2019  Jae-Joong Lee
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


n_col = size(cols, 1);
if n_col > 255
    cols = interp1(1:n_col, cols, linspace(1, n_col, 256));
    n_col = 256;
elseif n_col == 1
    cols = repmat(cols, 2, 1);
    n_col = 2;
end

if isfloat(cols) && max(cols(:)) <= 1 % given that color was specified as a scale of 0-to-1
    cols = round(cols .* 255);
    cols(cols > 255) = 255;
end
nodes = round(linspace(0, 255, n_col));

fid = fopen(clutpath, 'w');
fprintf(fid, '[FLT]\n');
fprintf(fid, 'min=0\n');
fprintf(fid, 'max=0\n');
fprintf(fid, '[INT]\n');
fprintf(fid, 'numnodes=%d\n', n_col);
fprintf(fid, '[BYT]\n');
for i = 1:n_col
    fprintf(fid, 'nodeintensity%d=%d\n', i-1, nodes(i));
end
fprintf(fid, '[RGBA255]\n');
for i = 1:n_col
    fprintf(fid, 'nodergba%d=%d|%d|%d|%d\n', i-1, cols(i,:), 255);
end
fclose(fid);

end