function out = circos_multilayer(A, varargin)

rotate_angle = 0;
lesion_proportions = zeros(size(A,1), 1);
edge_color = [215,25,28]./255;
% edge_color = [0 0 0]./255;
degree_color = [255,255,178
    254,204,92
    253,141,60
    240,59,32
    189,0,38]./255;
laterality = false;
radiological = false;
dot_node = 10;
dot_interval = 3;

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            case {'group'}
                group = varargin{i+1};
            case {'group_color'}
                gcols = varargin{i+1};
            case {'rotate'}
                rotate_angle = varargin{i+1};
            case {'lesion'}
                lesion_proportions = varargin{i+1};
                lesion_proportions = (lesion_proportions - min(lesion_proportions)) ./ (max(lesion_proportions) - min(lesion_proportions));
            case {'degree'}
                degree = varargin{i+1};
                degree = (degree - min(degree)) ./ (max(degree) - min(degree));
            case {'region_names'}
                region_names = varargin{i+1};
            case {'laterality'}
                laterality = true;
                lat_index = varargin{i+1};
            case {'radiological'}
                radiological = true;
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

lesion_proportions = lesion_proportions(group_idx);
degree = degree(group_idx);
region_names = region_names(group_idx);
    
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
    
    bottom_line = ref_line;
    top_line = flipud(ref_line) .* 1.05;
    patch_vec = [bottom_line; top_line];
    patch_color = gcols(group_val(i),:);
    patch([patch_vec(:,1)], [patch_vec(:,2)], patch_color, 'linewidth', 0.5, 'edgecolor', [.5 .5 .5], 'edgealpha', .5);
    
    bottom_line = flipud(ref_line) .* 1.05;
    top_line = ref_line .* 1.1;
    patch_vec = [bottom_line; top_line];
    patch_color = repmat(1-lesion_proportions(i), 1, 3);
    patch([patch_vec(:,1)], [patch_vec(:,2)], patch_color, 'linewidth', 0.5, 'edgecolor', [.5 .5 .5], 'edgealpha', .5);
    
    bottom_line = ref_line .* 1.1;
    top_line = flipud(ref_line) .* 1.15;
    patch_vec = [bottom_line; top_line];
    patch_color = degree_color(sum(degree(i) >= [0 0.2 0.4 0.6 0.8]), :);
    patch([patch_vec(:,1)], [patch_vec(:,2)], patch_color, 'linewidth', 0.5, 'edgecolor', [.5 .5 .5], 'edgealpha', .5);
    
    text_line = mean(ref_line) .* 1.16;
    text_rotate = rad2deg(mean(range_theta{i}));
    if text_rotate < -90
        text_rotate = text_rotate + 180;
        h = text(text_line(1), text_line(2), [region_names{i} '- '], 'HorizontalAlignment', 'Right', 'Fontsize', 8, 'Rotation', text_rotate);
    else
        h = text(text_line(1), text_line(2), [' -' region_names{i}], 'HorizontalAlignment', 'Left', 'Fontsize', 8, 'Rotation', text_rotate);
    end
end


%% Draw connections

[row,col,w] = find(triu(A,1));

[~, sorted_idx] = sort(w, 'descend');
alpha_w = zeros(size(w)) + 0.25;
alpha_w(sorted_idx(1:10)) = 1;
alpha_w(sorted_idx(11:50)) = 0.5;

% prctile_A = prctile(A(A~=0), [0 90 99]);
% alpha_vals = [0.1 0.2 1];
% alpha_A = zeros(size(A));
% for i = 1:numel(prctile_A)
%     alpha_A(A >= prctile_A(i)) = alpha_vals(i);
% end

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
    
    line(...
        r*cos(theta)+x0,...
        r*sin(theta)+y0,...
        'LineWidth', 2,...
        'PickableParts','none', 'color', [edge_color alpha_w(i)]);

end


end