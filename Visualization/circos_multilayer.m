function out = circos_multilayer(A, varargin)

rotate_angle = 0;
add_layer = {};
do_region_label = false;
% pos_edge_color = [227,26,28]./255;
pos_edge_color = [255,0,0]./255;
% neg_edge_color = [43,131,186]./255;
neg_edge_color = [10,150,255]./255;
region_names_size = 6;
laterality = false;
radiological = false;
sep_pos_neg = false;
dot_node = 10;
dot_interval = 3;
patch_size_coef = 0.05;
layer = {};
% alpha_fun = @(x) (((abs(x) - min(abs(x))) ./ (max(abs(x)) - min(abs(x))))).^4.5;
% width_fun = @(x) (abs(x) - min(abs(x))) ./ (max(abs(x)) - min(abs(x))) * 2.25 + 0.25;
alpha_fun = @(x) (x - min(x)) ./ (max(x) - min(x)) * 0.9 + 0.1;
width_fun = @(x) (abs(x) - min(abs(x))) ./ (max(abs(x)) - min(abs(x))) * 2 + 1;

default_col_names = { ...
    'degree', ...
    'lesion', ...
    'clcoef', ...
    };
default_col = { ...
    [255,255,178
    254,204,92
    253,141,60
    240,59,32
    189,0,38]./255, ...
    repmat(linspace(1, 0, 10)', 1, 3), ...
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
            case {'sep_pos_neg'}
                sep_pos_neg = true;
            case {'alpha_fun'}
                alpha_fun = varargin{i+1};
            case {'width_fun'}
                width_fun = varargin{i+1};
        end
    end
end

j = 0;
for i = 1:length(add_layer)
    if ischar(add_layer{i})
        switch add_layer{i}
            case {'layer'}
                j = j + 1;
                layer{j} = add_layer{i+1};
                if max(layer{j}) ~= min(layer{j})
                    layer{j} = (layer{j} - min(layer{j})) ./ (max(layer{j}) - min(layer{j}));
                end
            case {'color'}
                if ~ischar(add_layer{i+1})
                    layer_color{j} = add_layer{i+1};
                elseif ischar(add_layer{i+1})
                    layer_color{j} = default_col{contains(default_col_names, add_layer{i+1})};
                end
        end
    end
end




%% Calculating theta

% A = reformat_r_new(w, 'reconstruct'); %DELETE IT!!!

if laterality
    
    before_N_group = max(group);
    group(lat_index == 0) = before_N_group + 1;
    if ~radiological % Right is Right, Left is Left
        group(lat_index == -1) = before_N_group*2 + 2 - group(lat_index == -1);
    elseif radiological % Right is Left, Left is Right
        group(lat_index == 1) = before_N_group*2 + 2 - group(lat_index == 1);
    end
    gcols = [gcols; gcols(end,:); flipud(gcols)];

end

N_node = size(A, 1);
N_group = numel(unique(group));
unit_theta = (2*pi) / (N_node * dot_node + N_group * dot_interval);

[group_val, group_idx] = sort(group, 'ascend');
if laterality
    for i = (before_N_group+1):before_N_group*2+1
        wh_mirror = group_val == i;
        group_val(wh_mirror) = flipud(group_val(wh_mirror));
        group_idx(wh_mirror) = flipud(group_idx(wh_mirror));
    end
end
A = A(group_idx, group_idx);

for i = 1:numel(layer)
    layer{i} = layer{i}(group_idx);
end
if exist('region_names', 'var')
    region_names = region_names(group_idx);
end
    
wh_interval = find(diff([group_val]) == 1); % find where group index differs = find where interval is located

j = 0:(dot_node-1);
% j = [0:(dot_node-1)] + dot_interval;
for i = 1:N_node
    
    range_theta{i} = -(unit_theta * j) + pi/2 + deg2rad(rotate_angle);
    j = j + dot_node;
    if ismember(i, wh_interval)
        j = j + dot_interval; % interval
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
        patch([patch_vec(:,1)], [patch_vec(:,2)], patch_color, 'linewidth', 0.5, 'edgecolor', [.5 .5 .5], 'edgealpha', .5);
        
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
    patch_color = gcols(group_val(i),:);
    patch([patch_vec(:,1)], [patch_vec(:,2)], patch_color, 'linewidth', 0.5, 'edgecolor', [.5 .5 .5], 'edgealpha', .5);
    
    %% ROI text
    if do_region_label
        text_line = mean(ref_line) .* (1 + j*patch_size_coef);
        text_rotate = rad2deg(mean(range_theta{i}));
        if text_rotate < -90
            text_rotate = text_rotate + 180;
            h = text(text_line(1), text_line(2), [region_names{i} '- '], 'HorizontalAlignment', 'Right', 'Fontsize', region_names_size, 'Rotation', text_rotate);
        else
            h = text(text_line(1), text_line(2), [' -' region_names{i}], 'HorizontalAlignment', 'Left', 'Fontsize', region_names_size, 'Rotation', text_rotate);
        end
    end
    
end


%% Draw connections

[row,col,w] = find(triu(A,1));

% alpha_w = (w - min(w)) ./ (max(w) - min(w)) * 0.9 + 0.1;
alpha_w = alpha_fun(w);
% width_w = (abs(w) - min(abs(w))) ./ (max(abs(w)) - min(abs(w))) * 2 + 1;
width_w = width_fun(w);

for i = 1:numel(w)
    
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
    
    if ~sep_pos_neg
        edge_color = pos_edge_color;
        line(...
            r*cos(theta)+x0,...
            r*sin(theta)+y0,...
            'LineWidth', width_w(i),...
            'PickableParts','none', 'color', [edge_color alpha_w(i)]);
    elseif sep_pos_neg
        if w(i) >= 0; edge_color = pos_edge_color;
        elseif w(i) < 0; edge_color = neg_edge_color;
        end
        line(...
            r*cos(theta)+x0,...
            r*sin(theta)+y0,...
            'LineWidth', width_w(i),...
            'PickableParts','none', 'color', [edge_color alpha_w(i)]);
    end

end

axis off
% set(gcf, 'color', 'w', 'position', [1   444   990   920]);

end