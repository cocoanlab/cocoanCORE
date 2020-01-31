function out = circos_wani(A, varargin)

% This function draws circos plot.
%
% :Usage:
% ::
%
%    out = circos_wani(a, varargin)
%
% :Features:
%
%  - can draw ...
%
% :Inputs:
%
%   **a:**
%        adjacency matrix (weighted or not)
%

% Find non-zero values of s and their indices

laterality = false;
pos_edge_color = [215,25,28]./255;
neg_edge_color = [43,131,186]./255;
weighted_node_size = false;
do_node_alpha = false;
draw_node_top  = false;
draw_circle = false;
draw_node_edge = false;
rotate_angle = 0;
edge_alpha = 0.5;
node_size = 150;
draw_circle_width = 8;

% default (no group)
group = 1:size(A,1);
g_order = group;
gcols = repmat([0 0 0], size(A,1),1);

if isequal(unique(A), [0;1]) % binary adjacency matrix
    pos_edge_color = [0 0 0];
end

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            % functional commands
            case {'laterality'}
                laterality = true;
                lat_index = varargin{i+1};
            case {'group'}
                group = varargin{i+1};
                if size(group,1) < size(group,2)
                    group = group';
                end
%                 g_order = 1:numel(unique(group)); % default
                g_order = unique(group); % modified for more general purpose (2020.01.25, J.J.)
            case {'group_color'}
                gcols = varargin{i+1};
                gcols_edge = gcols - .2;
                gcols_edge(gcols_edge<0) = 0;
            case {'weighted_node_size'}
                weighted_node_size = true;
            case {'node_alpha'}
                do_node_alpha = true;
                node_alpha = varargin{i+1};
            case {'edge_alpha'}
                edge_alpha = varargin{i+1};
            case {'draw_node_top'}
                draw_node_top = true;
            case {'draw_circle'}
                draw_circle = true;
            case {'draw_node_edge'}
                draw_node_edge = true;
            case {'rotate'}
                rotate_angle = varargin{i+1};
            case {'pos_edge_color'}
                pos_edge_color = varargin{i+1};
            case {'neg_edge_color'}
                neg_edge_color = varargin{i+1};
            case {'edge_color'} % same color for all edges
                pos_edge_color = varargin{i+1};
                neg_edge_color = varargin{i+1};
            case {'node_color'} % same color for all nodes
                gcols = varargin{i+1};
                if size(gcols,1) == 1
                    gcols = repmat(varargin{i+1}, size(A,1),1);
                end
            case {'node_size'}
                node_size = varargin{i+1};
            case {'node_edge_color'}
                gcols_edge = varargin{i+1};
                if size(gcols_edge, 1) == 1
                    gcols_edge = repmat(gcols_edge, size(A,1),1);
                end
            case {'draw_circle_width'}
                draw_circle_width = varargin;
        end
    end
end

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            case {'group_order'}
                g_order = varargin{i+1};
        end
    end
end

if laterality
    [lat_val, lat_sort] = sort(lat_index, 'descend');
    
    group = group(lat_sort);
    A = A(lat_sort, lat_sort);
    lat_index = lat_index(lat_sort);
    
    g_idx_total = [];
    for k = [1 -1]
        temp = group .* double(lat_val==k);
        for i = 1:numel(g_order)
            g_idx = find(temp==g_order(i));
            g_idx_total = [g_idx_total; g_idx];
        end
    end
else
    g_idx_total = [];
    for i = 1:numel(g_order)
        temp = group;
        g_idx = find(temp==g_order(i));
        g_idx_total = [g_idx_total; g_idx];
    end
end

group = group(g_idx_total);
A = A(g_idx_total, g_idx_total);
if laterality, lat_index = lat_index(g_idx_total); end

[row,col,w] = find(A);
sumA = sum(A);
norm_factor = prctile(abs(w), 80); 
% norm_factor2 = min(abs(sumA(sumA~=0)));

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            % functional commands
            case {'norm_factor'}
                norm_factor = varargin{i+1}; % you can use your customized normalizing factor
        end
    end
end

if laterality
    t1 = linspace(-pi/2, pi/2, sum(lat_index>0)+4)'; % theta for each node
    t1([1 2 sum(lat_index>0)+3 sum(lat_index>0)+4]) = [];
    t2 = linspace(-pi/2, pi/2, sum(lat_index<0)+4)'+pi; % theta for each node
    t2([1 2 sum(lat_index<0)+3 sum(lat_index<0)+4]) = [];
    t = [flipud(t1); t2];
else
    t = linspace(-pi, pi,length(A) + 1)'; % theta for each node
end

t = t - deg2rad(rotate_angle);

if draw_circle
    tt = linspace(-pi, pi,1000)'; % theta for each node
    plot(cos(tt), sin(tt), 'color', [.8 .8 .8], 'linewidth', draw_circle_width);
    hold on;
end

for i = 1:length(A) %find(sumA == 0)
    % scatter(cos(t(i)),sin(t(i)), 50, gcols(group(i),:),  'filled', 'MarkerFaceAlpha', .5);
    scatter(cos(t(i)),sin(t(i)), node_size, gcols(group(i),:),  'filled');
    hold on;
end

for i = find(sumA ~= 0)
    % norm_factor3 = (abs(sumA(i))-norm_factor2)./(90*(norm_factor-norm_factor2))+1;
    
    if weighted_node_size
        % scatter(cos(t(i)).*norm_factor3,sin(t(i)).*norm_factor3, 150*abs(sumA(i))./norm_factor, gcols(group(i),:), 'filled');
        if draw_node_edge
            scatter(cos(t(i)), sin(t(i)), node_size*abs(sumA(i))./norm_factor, gcols(group(i),:), 'filled', 'MarkerEdgeColor', gcols_edge(group(i),:), 'LineWidth', 1);
        else
            scatter(cos(t(i)), sin(t(i)), node_size*abs(sumA(i))./norm_factor, gcols(group(i),:), 'filled');
        end
    else
        if draw_node_edge
            scatter(cos(t(i)), sin(t(i)), node_size, gcols(group(i),:), 'filled', 'MarkerEdgeColor', gcols_edge(group(i),:), 'LineWidth', 1);
        else
            scatter(cos(t(i)), sin(t(i)), node_size, gcols(group(i),:), 'filled');
        end
    end
    
    if do_node_alpha
        scatter(cos(t(i)),sin(t(i)), 300, 'w', ...gcols(group(i),:), ...
            'filled', 'MarkerFaceAlpha', node_alpha);
    end
    hold on;
end

% Calculate line widths based on values of s (stored in v).
minLineWidth  = 0.5;
lineWidthCoef = 5;
lineWidth = w./norm_factor;

if sum(lineWidth) == numel(lineWidth) % all lines are the same width.
    lineWidth = repmat(minLineWidth,numel(lineWidth),1);
else % lines of variable width.
    lineWidth = lineWidthCoef*lineWidth; % plus and minus
end

for i = 1:length(w)
    if row(i) ~= col(i)
        if abs(row(i) - col(i)) - length(A)/2 == 0
            % points are diametric, so draw a straight line
            u = [cos(t(row(i)));sin(t(row(i)))];
            v = [cos(t(col(i)));sin(t(col(i)))];
            if lineWidth(i) > 0
                line(...
                    [u(1);v(1)],...
                    [u(2);v(2)],...
                    'LineWidth', lineWidth(i),...
                    'PickableParts','none', 'color', [pos_edge_color edge_alpha]);
            else
                line(...
                    [u(1);v(1)],...
                    [u(2);v(2)],...
                    'LineWidth', -lineWidth(i),...
                    'PickableParts','none', 'color', [neg_edge_color edge_alpha]);
            end
        else % points are not diametric, so draw an arc
            u  = [cos(t(row(i)));sin(t(row(i)))];
            v  = [cos(t(col(i)));sin(t(col(i)))];
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
            
            %this.Node(row(i)).Connection(end+1) =
            if lineWidth(i) > 0
                line(...
                    r*cos(theta)+x0,...
                    r*sin(theta)+y0,...
                    'LineWidth', lineWidth(i),...
                    'PickableParts','none', 'color', [pos_edge_color edge_alpha]);
            else
                line(...
                    r*cos(theta)+x0,...
                    r*sin(theta)+y0,...
                    'LineWidth', -lineWidth(i),...
                    'PickableParts','none', 'color', [neg_edge_color edge_alpha]);
            end
        end
    end
end

if draw_node_top
    
    for i = find(sumA ~= 0)
        % norm_factor3 = (abs(sumA(i))-norm_factor2)./(90*(norm_factor-norm_factor2))+1;
        
        if weighted_node_size
            % scatter(cos(t(i)).*norm_factor3,sin(t(i)).*norm_factor3, 150*abs(sumA(i))./norm_factor, gcols(group(i),:), 'filled');
            if draw_node_edge
                scatter(cos(t(i)), sin(t(i)), node_size*abs(sumA(i))./norm_factor, gcols(group(i),:), 'filled', 'MarkerEdgeColor', gcols_edge(group(i),:), 'LineWidth', 1);
            else
                scatter(cos(t(i)), sin(t(i)), node_size*abs(sumA(i))./norm_factor, gcols(group(i),:), 'filled');
            end
        else
            if draw_node_edge
                scatter(cos(t(i)), sin(t(i)), node_size, gcols(group(i),:), 'filled', 'MarkerEdgeColor', gcols_edge(group(i),:), 'LineWidth', 1);
            else
                scatter(cos(t(i)), sin(t(i)), node_size, gcols(group(i),:), 'filled');
            end
        end
        
        if do_node_alpha
            scatter(cos(t(i)),sin(t(i)), 300, 'w', ...gcols(group(i),:), ...
                'filled', 'MarkerFaceAlpha', node_alpha);
        end
        hold on;
    end
    
    
end

set(gca, 'xlim', [-1 1], 'ylim', [-1 1]); % [2020. 01. 31] added by J.J.: sometimes xlim starts from -2.
set(gcf, 'color', 'w', 'position', [1   444   562   511]);
axis off

out.norm_factor = norm_factor;

end