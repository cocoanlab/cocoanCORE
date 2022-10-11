function circos_multilayer(A, varargin)

% This function draws circos plot with multiple layers.
%
%
% :Usage:
% ::
%     circos_multilayer(A, varargin)
%
%
% :Input:
% ::
%   - A                  Adjacency matrix
%
%
% :Optional Input:
%
%   - group              group assignment of nodes. ex) [1,1,1,2,2,2,3,4], ...
%   - group_color        RGB color values for each group. [groups X 3]
%   - rotate             rotate circos plot counter-clockwise. (degree)
%   - length_ratio       relative length ratio between
%                        [region, interval between region, interval between groups]. (default: [1,0,1])
%   - patch_edge_alpha   alpha value for edge of region patch. (default: 0.5)
%   - patch_edge_color   color value for edge of region patch. (default: [0.5 0.5 0.5])
%   - patch_size_coef    size for region path. (default: 0.05)
%   - patch_resolution   number of sides for region patch. (default: 10)
%   - add_layer          add additional layers, specified by following key-value pairs.
%                        (e.g., 'add_layer', {'layer', deg_cent, 'color', deg_cent_cols})
%                        the values followed by 'layer' will be presented
%                        with the colormap values following 'color'.
%                        for example, if you want to add 'degree centrality'
%                        layer in the circos plot and define color with 5-level,
%                        put degree centrality value in 'deg_cent'
%                        and 5 RGB values [5 X 3] in 'deg_cent_cols'.
%                        (Note: all values for layer should be normalized,
%                        ranging from 0 to 1)
%   - region_names       labels for regions (or nodes)
%   - region_names_size  size of labels for regions
%   - laterality         laterality index for circos plot (usually for cortex)
%                        -1: Left, 1: Right, 0: No laterality
%   - radiological       laterality display in radiological convention. (default: neurological)
%   - alpha_fun          function handle of mapping rule for alpha value of connections.
%                        if conn_alpha is specified, this option will not be used.                         
%                        (default: ((abs(x) - min(abs(x))) ./ (max(abs(x)) - min(abs(x)))).^4.5)
%   - width_fun          function handle of mapping rule for width value of connections.
%                        if conn_width is specified, this option will not be used.
%                        (default: (abs(x) - min(abs(x))) ./ (max(abs(x)) - min(abs(x))) * 2.25 + 0.25)
%   - conn_color         color value of connections. [connections X 3] or [regions X regions X 3]
%                        (default: [255,0,0]./255 for positive w, [10,150,255]./255 for negative w)
%   - conn_alpha         alpha value of connections. [connections X 1] or [regions X regions]
%                        (default: based on alpha_fun)
%   - conn_width         width value of connections. [connections X 1] or [regions X regions]
%                        (default: based on width_fun)
%   - conn_order         order of drawing connections. [connections X 1] or [regions X regions]
%                        (default: order from find() function )
%   - each_patch_color   specify color of each patches(only for single patch
%                        ([number of nodes X 3])
%   - text_color         specify color of each color in 'region_names'
%                        ([number of nodes X 3] or [number of nodes X 'name of color in string'(e.g. 'red')])
% :Output:
% ::
%
%
% :Example:
% ::
%   % generating undirected and weighted adjacency matrix
%   A = randn(100,100); % to make realistic network, use other fx.
%   A = A + A';
%   A = A ./ max(abs(A(:)));
%   A(1:length(A)+1:end) = 0;
%   % defining groups of nodes
%   A_group = randi(5,100,1);
%   A_group_cols = [0    0.4470    0.7410
%       0.8500    0.3250    0.0980
%       0.9290    0.6940    0.1250
%       0.4940    0.1840    0.5560
%       0.4660    0.6740    0.1880];
%   % add 'weighted degree centrality' layer
%   A_pos_deg_cent = sum(A .* double(A>0), 2) ./ (size(A, 1) - 1);
%   A_pos_deg_cent_cols = cell2mat(arrayfun(@(x, y) linspace(x, y, 10), ...
%       [179,0,0], [254,240,17], 'UniformOutput', false)')' ./ 255;
%   A_neg_deg_cent = -sum(A .* double(A<0), 2) ./ (size(A, 1) - 1);
%   A_neg_deg_cent_cols = cell2mat(arrayfun(@(x, y) linspace(x, y, 10), ...
%       [8,104,172], [240,290,232], 'UniformOutput', false)')' ./ 255;
%   % region names
%   A_names = cellstr(strcat(repmat('node', 100, 1), num2str((1:100)')));
%   % laterality
%   A_lat = [repmat([-1; 1], 45, 1); zeros(10,1)];
%   % draw circos plot
%   circos_multilayer(A, 'group', A_group, 'group_color', A_group_cols, ...
%       'add_layer', {'layer', A_pos_deg_cent, 'color', A_pos_deg_cent_cols, ...
%       'layer', A_neg_deg_cent, 'color', A_neg_deg_cent_cols}, ...
%       'region_names', A_names, 'laterality', A_lat);
%
%   %% note history
%   2020/04/13 Jungwoo: add 'each_patch_color', 'text_color' option & 'Interpretor', 'None' for 'region_names' option
%
%     Author and copyright information:
%
%     Copyright (C) Jan 2019  Choong-Wan Woo & Jae-Joong Lee
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

group = ones(size(A,1), 1);
gcols = [0.5 0.5 0.5];
rotate_angle = 0;
length_region = 1;
interval_region = 0;
interval_group = 1;
patch_edge_alpha = 0.5;
patch_edge_color = [0.5 0.5 0.5];
patch_size_coef = 0.05;
patch_resolution = 10;
add_layer = {};
do_region_label = false;
region_names_size = 6;
laterality = false;
radiological = false;
alpha_fun = @(x) ((abs(x) - min(abs(x))) ./ (max(abs(x)) - min(abs(x)))).^4.5;
width_fun = @(x) (abs(x) - min(abs(x))) ./ (max(abs(x)) - min(abs(x))) * 2.25 + 0.25;
sep_pos_neg = false;
layer = {};

default_col_names = { ...
    'degree', ...
    'lesion1', ...
    'lesion2', ...
    'lesion3', ...
    'clcoef', ...
    };
default_col = { ...
    [255,255,178
    254,204,92
    253,141,60
    240,59,32
    189,0,38]./255, ...
    repmat(linspace(1, 0, 10)', 1, 3), ...
    [247,252,253
    224,236,244
    191,211,230
    158,188,218
    140,150,198
    140,107,177
    136,65,157
    129,15,124
    77,0,75]./255, ...
    [255,247,243
    253,224,221
    252,197,192
    250,159,181
    247,104,161
    221,52,151
    174,1,126
    122,1,119
    73,0,106]./255, ...
    [237,248,251
    178,226,226
    102,194,164
    44,162,95
    0,109,44]./255, ...
    };

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            case {'group'}
                group = varargin{i+1};
            case {'group_color'}
                gcols = varargin{i+1};
            case {'rotate'}
                rotate_angle = varargin{i+1};
            case {'length_ratio'}
                length_region = varargin{i+1}(1);
                interval_region = varargin{i+1}(2);
                interval_group = varargin{i+1}(3);
            case {'patch_edge_alpha'}
                patch_edge_alpha = varargin{i+1};
            case {'patch_edge_color'}
                patch_edge_color = varargin{i+1};
            case {'patch_size_coef'}
                patch_size_coef = varargin{i+1};
            case {'patch_resolution'}
                patch_resolution = varargin{i+1};
            case {'add_layer'}
                add_layer = varargin{i+1};
            case {'region_names'}
                do_region_label = true;
                region_names = varargin{i+1};
            case {'region_names_size'}
                region_names_size = varargin{i+1};
            case {'laterality'}
                laterality = true;
                lat_index = varargin{i+1};
            case {'radiological'}
                radiological = true;
            case {'alpha_fun'}
                alpha_fun = varargin{i+1};
            case {'width_fun'}
                width_fun = varargin{i+1};
            case {'conn_color'}
                conn_color = varargin{i+1};
                sz = size(conn_color);
                if numel(sz) == 3 && sz(1) == sz(2)
                    conn_color = conn_color(repmat(triu(true(sz(1),sz(2)), 1), 1,1,3));
                    conn_color = reshape(conn_color, [], 3);
                end
            case {'conn_alpha'}
                conn_alpha = varargin{i+1};
                sz = size(conn_alpha);
                if sz(1) == sz(2) && all(sz(1:2) > 1)
                    conn_alpha = conn_alpha(triu(true(sz(1),sz(2)), 1));
                end
            case {'conn_width'}
                conn_width = varargin{i+1};
                sz = size(conn_width);
                if sz(1) == sz(2) && all(sz(1:2) > 1)
                    conn_width = conn_width(triu(true(sz(1),sz(2)), 1));
                end
            case {'conn_order'}
                conn_order = varargin{i+1};
                sz = size(conn_order);
                if sz(1) == sz(2) && all(sz(1:2) > 1)
                    conn_order = conn_order(triu(true(sz(1),sz(2)), 1));
                end
            case {'each_patch_color'}
                patch_cmap = varargin{i+1};
            case {'text_color'}
                text_color = varargin{i+1};
        end
    end
end

j = 0;
for i = 1:length(add_layer)
    if ischar(add_layer{i})
        switch add_layer{i}
            case {'layer'}
                if ~isempty(add_layer{i+1})
                    j = j + 1;
                    layer{j} = add_layer{i+1};
                    if max(layer{j}) > 1 || min(layer{j}) < 0
                        error('Values of each layer should be between 0 and 1.');
                    end
                end
            case {'color'}
                if ~isempty(add_layer{i-1})
                    if ~ischar(add_layer{i+1})
                        layer_color{j} = add_layer{i+1};
                    elseif ischar(add_layer{i+1})
                        layer_color{j} = default_col{contains(default_col_names, add_layer{i+1})};
                    end
                end
        end
    end
end




%% Calculating theta

[~,~,group] = unique(group);
gcols = gcols(unique(group),:);

if laterality
    
    [ugroup_L,~,group_L] = unique(group(lat_index == -1));
    [ugroup_M,~,group_M] = unique(group(lat_index == 0));
    [ugroup_R,~,group_R] = unique(group(lat_index == 1));
    
    if ~radiological % Right is Right, Left is Left
        group_new = group;
        group_new(lat_index == 1) = group_R;
        group_new(lat_index == 0) = group_M + max(group_R);
        group_new(lat_index == -1) = group_L + max(group_R) + max(group_M);
        group = group_new;
        gcols = gcols([ugroup_R; ugroup_M; ugroup_L],:);
    elseif radiological % Right is Left, Left is Right
        group_new = group;
        group_new(lat_index == -1) = group_L;
        group_new(lat_index == 0) = group_M + max(group_L);
        group_new(lat_index == 1) = group_R + max(group_L) + max(group_M);
        group = group_new;
        gcols = gcols([ugroup_L; ugroup_M; ugroup_R],:);
    end
    
end

[group_val, group_idx] = sort(group, 'ascend');
if laterality
    if ~radiological % Right is Right, Left is Left
        wh_mirror = group_val > max(group_R) + max(group_M);
    elseif radiological % Right is Left, Left is Right
        wh_mirror = group_val > max(group_L) + max(group_M);
    end
    group_val(wh_mirror) = flipud(group_val(wh_mirror));
    group_idx(wh_mirror) = flipud(group_idx(wh_mirror));
end

orig_to_new_idx_mat = double(triu(logical(A),1));
orig_to_new_idx_mat(orig_to_new_idx_mat ~= 0) = 1:sum(sum(orig_to_new_idx_mat ~= 0));
orig_to_new_idx_mat = orig_to_new_idx_mat + orig_to_new_idx_mat';
orig_to_new_idx_mat = orig_to_new_idx_mat(group_idx, group_idx);
[~, ~, orig_to_new_idx] = find(triu(orig_to_new_idx_mat,1));

A = A(group_idx, group_idx);


for i = 1:numel(layer)
    layer{i} = layer{i}(group_idx);
end
if exist('region_names', 'var')
    region_names = region_names(group_idx);
end


N_node = size(A, 1);
N_group = numel(unique(group));
unit_theta = (2*pi) / (N_node * (length_region + interval_region) + N_group * interval_group);

wh_interval = find(diff([group_val]) ~= 0); % find where group index differs = find where interval is located

j = [0 length_region];
for i = 1:N_node
    
    range_theta{i} = -(unit_theta * j) + pi/2 + deg2rad(rotate_angle);
    range_theta{i} = linspace(range_theta{i}(1), range_theta{i}(2), patch_resolution); 
    j = j + length_region + interval_region;
    if ismember(i, wh_interval)
        j = j + interval_group; % interval
    end
    
end

%% Draw circular sectors

for i = 1:N_node
    
    ref_line = [cos(range_theta{i})', sin(range_theta{i})'];
    
    %% Layer
    for j = 1:length(layer)
        
        bottom_line = ref_line .* (1 + (j-1)*patch_size_coef);
        top_line = ref_line .* (1 + j*patch_size_coef);
        if mod(j, 2) == 1
            top_line = flipud(top_line);
        elseif mod(j, 2) == 0
            bottom_line = flipud(bottom_line);
        end
        patch_vec = [bottom_line; top_line];
        patch_color = layer_color{j}(sum(layer{j}(i) >= linspace(0, 1+eps, size(layer_color{j},1) + 1)), :);
        patch([patch_vec(:,1)], [patch_vec(:,2)], patch_color, 'linewidth', 0.5, 'edgecolor', patch_edge_color, 'edgealpha', patch_edge_alpha);
        
    end
    
    %% Group color
    if isempty(j); j = 1;
    elseif ~isempty(j); j = j + 1;
    end
    bottom_line = ref_line .* (1 + (j-1)*patch_size_coef);
    top_line = ref_line .* (1 + j*patch_size_coef);
    if mod(j, 2) == 1
        top_line = flipud(top_line);
    elseif mod(j, 2) == 0
        bottom_line = flipud(bottom_line);
    end
    patch_vec = [bottom_line; top_line];
    if exist('patch_cmap', 'var'), patch_color = patch_cmap(i, :);
    else, patch_color = gcols(group_val(i),:);end
    patch([patch_vec(:,1)], [patch_vec(:,2)], patch_color, 'linewidth', 0.5, 'edgecolor', patch_edge_color, 'edgealpha', patch_edge_alpha);
    
    %% ROI text
    if do_region_label
        text_line = mean(ref_line) .* (1 + j*patch_size_coef);
        text_rotate = rad2deg(mean(range_theta{i}));
        if text_rotate < -90 || text_rotate > 90
            text_rotate = text_rotate + 180;
            % h = text(text_line(1), text_line(2), [region_names{i} '- '], 'HorizontalAlignment', 'Right', 'Fontsize', region_names_size, 'Rotation', text_rotate);
            if exist('text_color', 'var'), text_col = text_color(i, :);
                h = text(text_line(1), text_line(2), [region_names{i} '- '], 'HorizontalAlignment', 'Right', 'Fontsize', region_names_size, ...
                    'Rotation', text_rotate, 'Interpreter', 'none', 'Color', text_col);
            else
                h = text(text_line(1), text_line(2), [region_names{i} '- '], 'HorizontalAlignment', 'Right', 'Fontsize', region_names_size, ...
                    'Rotation', text_rotate, 'Interpreter', 'none');
            end
        else
            % h = text(text_line(1), text_line(2), [' -' region_names{i}], 'HorizontalAlignment', 'Left', 'Fontsize', region_names_size, 'Rotation', text_rotate);
            if exist('text_color', 'var'), text_col = text_color(i, :);
                h = text(text_line(1), text_line(2), [' -' region_names{i}], 'HorizontalAlignment', 'Left', 'Fontsize', region_names_size, ...
                    'Rotation', text_rotate, 'Interpreter', 'none', 'Color', text_col);
            else
                h = text(text_line(1), text_line(2), [' -' region_names{i}], 'HorizontalAlignment', 'Left', 'Fontsize', region_names_size, ...
                    'Rotation', text_rotate, 'Interpreter', 'none');
            end
        end
    end
    
end


%% Draw connections

[row,col,w] = find(triu(A,1));

if ~exist('conn_color', 'var')
    pos_edge_color = [255,0,0]./255;
    neg_edge_color = [10,150,255]./255;
    conn_color = zeros(numel(w), 3);
    conn_color(w>0, :) = repmat(pos_edge_color, sum(w>0), 1);
    conn_color(w<0, :) = repmat(neg_edge_color, sum(w<0), 1);
else
    conn_color = conn_color(orig_to_new_idx, :);
end

if ~exist('conn_alpha', 'var')
    conn_alpha = alpha_fun(w);
else
    conn_alpha = conn_alpha(orig_to_new_idx);
end

if ~exist('conn_width', 'var')
    conn_width = width_fun(w);
else
    conn_width = conn_width(orig_to_new_idx);
end

if ~exist('conn_order', 'var')
    conn_order = 1:numel(w);
else
    conn_order = conn_order(orig_to_new_idx);
end
if iscolumn(conn_order); conn_order = conn_order'; end

for i = conn_order
    
    u = [cos(mean(range_theta{row(i)})), sin(mean(range_theta{row(i)}))];
    v = [cos(mean(range_theta{col(i)})), sin(mean(range_theta{col(i)}))];
    
    x0 = -(u(2)-v(2))/(u(1)*v(2)-u(2)*v(1));
    y0 =  (u(1)-v(1))/(u(1)*v(2)-u(2)*v(1));
    r  = sqrt(x0^2 + y0^2 - 1);
    thetaLim(1) = atan2(u(2)-y0,u(1)-x0);
    thetaLim(2) = atan2(v(2)-y0,v(1)-x0);
    
    if u(1) >= 0 && v(1) >= 0
        % ensure the arc is within the unit disk
        theta = [linspace(max(thetaLim),pi,50),...
            linspace(-pi,min(thetaLim),50)].';
    else
        theta = linspace(thetaLim(1),thetaLim(2)).';
    end
    
    if abs(u*v') <= 1+eps && abs(u*v') >= 1-eps  % u and v are almost the same
        line(...
            [u(1) v(1)],...
            [u(2) v(2)],...
            'LineWidth', conn_width(i),...
            'PickableParts','none', 'color', [conn_color(i,:) conn_alpha(i)]);
    else
        line(...
            r*cos(theta)+x0,...
            r*sin(theta)+y0,...
            'LineWidth', conn_width(i),...
            'PickableParts','none', 'color', [conn_color(i,:) conn_alpha(i)]);
    end
    
end

axis off;
set(gcf, 'color', 'w');

end