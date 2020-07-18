function h = color_check(colors)

% This simple function show the colors
%
% :Usage:
% ::
%     h = color_check(colors)
%
% :Input:
% ::
%   - colors             color scheme. (N x 3)
%
% :Output:
% ::
%   - h                  figure handle
%
% :Example:
% ::
%
% colors = [0.9098    0.4902    0.4471;
%     0.3059    0.7098    0.9020];
% h = color_check(colors);
%
% % if you want to make this colors into CLUT for mricrogl, use make_CLUT.m
% clutpath{1} = '/Applications/MRIcroGL.app/Contents/Resources/lut/semic_stim.clut';
% clutpath{2} = '/Applications/MRIcroGL.app/Contents/Resources/lut/semic_cue.clut';
% for i = 1:2
%     make_CLUT(colors(i,:), clutpath{i})
% end
%
%
% ..
%     Author and copyright information:
%
%     Copyright (C) Oct 2020  Choong-Wan Woo
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


h = figure;
for i = 1:size(colors,1)
    patch(i+[0 0 1 1], [0 1 1 0],colors(i,:),'EdgeColor','none');
    hold on;
end

axis off;
set(gcf, 'color', 'w');