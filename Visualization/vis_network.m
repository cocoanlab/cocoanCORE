function h = vis_network(W, varargin)

% Visualize network
%
% :Usage:
% ::
%     h = vis_network(W, varargin)
%
% ..
%     Author and copyright information:
%
%     Copyright (C) 2018  Wani Woo
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
%   **W:**
%        adjacency matrix (can be weighted)
%
% :Optional Inputs:
%
%   **'weighted':**
%        draw weighted network (positive: red, negative: blue)
%        weights are normalized between 1 and 4
%
%   **'degree':**
%        Make the node size proportional to degree (binary) (normalized
%        between 50-450)
%
%   **'weighted_absolute_deg':**
%        Make the node size proportional to the sum of absolue weights 
%        (normalized between 50-450)
%
%   **'weighted_raw_deg':**
%        Make the node size proportional to the sum of raw weights
%        (normalized between 50-450)
% 
%   **'group':**
%        Use different colors for different group members - this needs to
%        be used along with 'groupcolor' option
%        e.g., 'group', [1 1 1 2 2 3 3 3 4 4]
%
%   **'groupcolor'
%        group colors, color matrix (its row should be # of groups)
%        e.g., 'groupcolor', [0.9098    0.4902    0.4471; 0.7451    0.6039 0.2000; 0.4275    0.6863    0.2039; 0.3373    0.7373    0.5922]
%  
%   **'highlight'
%        you can use this option when you want to highlight some nodes. You
%        have to provide a vector of logical values (true or false) to
%        specify which nodes you want to highlight
%        e.g., 'highlight', [1 1 0 0 0 1 1 0 0 0]
%
%   **'label'
%       Add each node name using 'text' function (ln: 297)
%       e.g., 'label', {'1','2','3','4','5'}
%
% :Outputs:
%
%
%   **h:**
%        graphic handles
%
% :Examples:
% ::
% % create some fake W
% W = rand(10,10);
% W = W-0.5;
% W = reformat_r_new(W, 'symmetric_avg', 'remove_diag');
% W = W.*double(abs(W)>.1);
%
% grouping = [1 1 1 2 2 3 3 3 4 4];
% g_cols = [0.9098    0.4902    0.4471
%     0.7451    0.6039    0.2000
%     0.4275    0.6863    0.2039
%     0.3373    0.7373    0.5922];
% labelname = {'1','2','3','4','5','6','7','8','9','10'}; 
%
% h = vis_network(W, 'weighted', 'degree', 'group', grouping, 'groupcolor', g_cols)
% h = vis_network(W, 'weighted', 'degree', 'group', grouping, 'groupcolor', g_cols,'label',labelname)
%
% ..
%    Programmers' notes:
%    Created 8/28/18 by wani woo
%
%    08/29/18 : J.J.   - change to use only upper triangle because the input
%                      should not be duplicate matrix.
%    01/30/20 : Suhwan - add label option 
%    05/30/20 : Suhwan - add gravity option 
% ..
%


doweight = 0;
dogravity = false; 
do3d = 0;
dodegree = 0; % change size of the node using degree
do_w_degree = 0;
do_absw_degree = 0;
dogroup = 0;
do_manual_node_size = 0;
g_cols = {};
wh_hl = true(size(W,1),1);
noline = 0;
ln_pos_color = [0.8431    0.0980    0.1098];
ln_neg_color = [0.1686    0.5137    0.7294];
ln_width = 1.5;
do_arrow = false;
arrow_len = 10;
do_label = false;
node_edge_color = 'w';
xy = [];
xystart = [];
do_trans = false;
node_lwidth = 2;

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            % functional commands
            case {'weighted', 'weight'}
                doweight = 1;
            case {'gravity'}
                dogravity = true;
            case {'3d'} % not yet implemented
                do3d = 1;
                cl = varargin{i+1};
            case {'degree'}
                dodegree = 1;
            case {'weighted_absolute_deg'}
                do_absw_degree = 1;
            case {'weighted_raw_deg'}
                do_w_degree = 1;
            case {'group'}
                dogroup = 1;
                grouping = varargin{i+1};
                if size(grouping,1) == 1
                    grouping = grouping';
                end
            case {'groupcolor'}
                g_cols = varargin{i+1};
            case {'highlight'}
                wh_hl = varargin{i+1};
            case {'manual_node_size'}
                do_manual_node_size = 1;
                node_size = varargin{i+1};
            case {'noline'}
                noline = 1;
            case {'line_pos_color'}
                ln_pos_color = varargin{i+1};
            case {'line_neg_color'}
                ln_neg_color = varargin{i+1};
            case {'linewidth'}
                ln_width = varargin{i+1};
            case {'arrow'}
                do_arrow = true;
            case {'arrow_length'}
                arrow_len = varargin{i+1};
            case {'label', 'labels'}
                do_label = 1;
                labelname = varargin{i+1};
            case {'node_edge_color'}
                node_edge_color = varargin{i+1};
            case {'xy'}
                xy = varargin{i+1};
            case {'xystart'}
                xystart = varargin{i+1};
                if isa(xystart, 'region')
                    xy_C = cell2mat({xystart.center}');
                    xy_D = repmat(xy_C, 1, 1, size(xy_C,1));
                    xy_D = xy_D - permute(xy_D, [3 2 1]);
                    xy_D = squeeze(sum(xy_D .^ 2, 2) .^ 0.5);
                    xystart = cmdscale(xy_D,2);
                end
            case {'trans'}
                do_trans = true;
        end
    end
end

W(isnan(W)) = 0;
if issymmetric(W)
    W = reformat_r_new(W, 'upper_triangle'); 
end

[i, j, s] = find(W);

new_s = (abs(s) - min(abs(s)))./(max(abs(s))-min(abs(s))); % 0-1 normalization
new_s = new_s.*4+1; % 1~4
trans_val = ((new_s-min(new_s))./(max(new_s)-min(new_s)))*0.9+0.1; % 0.1-1 normalization

if doweight
    G = graph(i,j,s);
else
    G = graph(i,j);
end

% add the empty nodes at the end of the network
G = addnode(G, size(W,1)-max(max(i), max(j)));

hh = figure;

if doweight
    if ~isempty(xystart)
        h_graph = plot(G, 'Layout', 'force', 'WeightEffect', 'inverse', 'UseGravity', dogravity, 'XStart', xystart(:,1), 'YStart', xystart(:,2));
    elseif isempty(xystart)
        h_graph = plot(G, 'Layout', 'force', 'WeightEffect', 'inverse', 'UseGravity', dogravity);
    end
else
    if ~isempty(xystart)
        h_graph = plot(G, 'Layout', 'force', 'UseGravity', dogravity, 'XStart', xystart(:,1), 'YStart', xystart(:,2));
    elseif isempty(xystart)
        h_graph = plot(G, 'Layout', 'force', 'UseGravity', dogravity);
    end
end
X = [h_graph.XData',h_graph.YData'];

if ~isempty(xy)
    X = xy;
end

close(hh);

if do3d
    for i = 1:numel(cl), xyz(i,:) = cl(i).mm_center; end
end

if ~noline
    if do_arrow
        
        [neg_i,neg_j] = find(W < 0);
        
        for k = 1:numel(neg_i)

            x = [X(neg_i(k),1) X(neg_j(k),1)];
            y = [X(neg_i(k),2) X(neg_j(k),2)];
            
            [newx, newy] = prop_reduce(x,y, .7);

            h.edge_neg(k) = arrow([newx(1) newy(1)], [newx(2) newy(2)], 'length', arrow_len);
            if doweight
                set(h.edge_neg(k), 'edgecolor', ln_neg_color, 'facecolor', ln_neg_color,  'linewidth', new_s(sum([i==neg_i(k) j==neg_j(k)],2)==2));
            else
                set(h.edge_neg(k), 'edgecolor', ln_neg_color, 'facecolor', ln_neg_color, 'linewidth', ln_width);
            end
            hold on;
        end
        
        [pos_i,pos_j] = find(W > 0);
        
        for k = 1:numel(pos_i)
            
            x = [X(pos_i(k),1) X(pos_j(k),1)];
            y = [X(pos_i(k),2) X(pos_j(k),2)];
            
            [newx, newy] = prop_reduce(x,y, .7);
            
            h.edge_pos(k) = arrow([newx(1) newy(1)], [newx(2) newy(2)], 'length', arrow_len);            
            if doweight
                set(h.edge_pos(k), 'edgecolor', ln_pos_color, 'facecolor', ln_pos_color, 'linewidth', new_s(sum([i==pos_i(k) j==pos_j(k)],2)==2));
            else
                set(h.edge_pos(k), 'edgecolor', ln_pos_color, 'facecolor', ln_pos_color, 'linewidth', ln_width);
            end
            hold on;
        end
    else
        [neg_i,neg_j] = find(W < 0);
        
        for k = 1:numel(neg_i)
            h.edge_neg(k) = line([X(neg_i(k),1) X(neg_j(k),1)], [X(neg_i(k),2) X(neg_j(k),2)]);
            if doweight
                e_idx = sum([i==neg_i(k) j==neg_j(k)],2)==2;
                if do_trans
                    set(h.edge_neg(k), 'color', [ln_neg_color trans_val(e_idx)], 'linewidth', new_s(e_idx));
                else
                    set(h.edge_neg(k), 'color', ln_neg_color, 'linewidth', new_s(e_idx));
                end
            else
                if do_trans
                    set(h.edge_neg(k), 'color', [ln_neg_color .5], 'linewidth', ln_width);
                else
                    set(h.edge_neg(k), 'color', ln_neg_color, 'linewidth', ln_width);
                end
            end
            hold on;
        end
        
        [pos_i,pos_j] = find(W > 0);
        
        for k = 1:numel(pos_i)
            h.edge_pos(k) = line([X(pos_i(k),1) X(pos_j(k),1)], [X(pos_i(k),2) X(pos_j(k),2)]);
            if doweight
                e_idx = sum([i==pos_i(k) j==pos_j(k)],2)==2;
                if do_trans
                    set(h.edge_pos(k), 'color', [ln_pos_color trans_val(e_idx)], 'linewidth', new_s(e_idx));
                else
                    set(h.edge_pos(k), 'color', ln_pos_color, 'linewidth', new_s(e_idx));
                end
            else
                if do_trans
                    set(h.edge_pos(k), 'color', [ln_pos_color .5], 'linewidth', ln_width);
                else
                    set(h.edge_pos(k), 'color', ln_pos_color, 'linewidth', ln_width);
                end
            end
            hold on;
        end
    end
end

W = reformat_r_new(W, 'symmetric_sum');
if dodegree
    d = sum(W~=0);
    if all(d==d(1)), d = repmat(150, size(d));
    else, d = ((d-min(d))./(max(d)-min(d))*4+.5)*100; % 50-450
    end
elseif do_absw_degree
    d = sum(abs(W));
    if all(d==d(1)), d = repmat(150, size(d));
    else, d = ((d-min(d))./(max(d)-min(d))*4+.5)*100; % 50-450
    end
elseif do_w_degree
    d = sum(W);
    if all(d==d(1)), d = repmat(150, size(d));
    else, d = ((d-min(d))./(max(d)-min(d))*4+.5)*100; % 50-450
    end
elseif do_manual_node_size
    d = node_size;
else
    d = repmat(150, size(X,1), 1);
end

% for node_i = 1:size(X,1)
%     if wh_hl(node_i)
%         if dogroup
%             h.node{node_i} = scatter(X(node_i,1), X(node_i,2), d(node_i), 'filled', 'MarkerFaceColor', g_cols(grouping(node_i),:), 'MarkerEdgeColor', 'w' , 'LineWidth', 1.5);
%             hold on;
%         else
%             h.node{node_i} = scatter(X(node_i,1), X(node_i,2), d(node_i), 'filled', 'MarkerFaceColor', [.5 .5 .5], 'MarkerEdgeColor', 'w' , 'LineWidth', 1.5);
%             hold on;
%         end
%     else
%         h.node{node_i} = scatter(X(node_i,1), X(node_i,2), 50, 'filled', 'MarkerFaceColor', [.5 .5 .5], 'MarkerEdgeColor', 'w' , 'LineWidth', 1.5);
%         hold on;
%     end
% end

h.node_xy = NaN(size(X,1),2);

if dogroup
    for g_i = unique(grouping)'
        h.node_hl{g_i} = scatter(X(wh_hl & grouping==g_i,1), X(wh_hl & grouping==g_i,2), d(wh_hl & grouping==g_i), 'filled', 'MarkerFaceColor', g_cols(g_i,:), 'MarkerEdgeColor', node_edge_color , 'LineWidth', node_lwidth);
        h.node_xy(wh_hl & grouping==g_i,1) = h.node_hl{g_i}.XData;
        h.node_xy(wh_hl & grouping==g_i,2) = h.node_hl{g_i}.YData;
        hold on;
    end
    if any(wh_hl == 0)
        h.node_nohl = scatter(X(~wh_hl,1), X(~wh_hl,2), 50, 'filled', 'MarkerFaceColor', [.5 .5 .5], 'MarkerEdgeColor', node_edge_color , 'LineWidth', node_lwidth);
        h.node_xy(~wh_hl,1) = h.node_nohl.XData;
        h.node_xy(~wh_hl,2) = h.node_nohl.YData;
        hold on;
    end
else
    h.node_hl = scatter(X(wh_hl,1), X(wh_hl,2), d(wh_hl), 'filled', 'MarkerFaceColor', [.5 .5 .5], 'MarkerEdgeColor', node_edge_color , 'LineWidth', node_lwidth);
    h.node_xy(wh_hl,1) = h.node_hl.XData;
    h.node_xy(wh_hl,2) = h.node_hl.YData;
    hold on;
    if any(wh_hl == 0)
        h.node_nohl = scatter(X(~wh_hl,1), X(~wh_hl,2), 50, 'filled', 'MarkerFaceColor', [.5 .5 .5], 'MarkerEdgeColor', node_edge_color , 'LineWidth', node_lwidth);
        h.node_xy(~wh_hl,1) = h.node_nohl.XData;
        h.node_xy(~wh_hl,2) = h.node_nohl.YData;
        hold on;
    end
end

if do_label
    label_space = ((max(X(:,1))-min(X(:,1))))*0.02; 
    text(X(:,1)+label_space , X(:,2),labelname);
end

axis off;
set(gcf, 'position', [360   235   570   463], 'color', 'w');

end

function [newx, newy] = prop_reduce(x,y, prop)

xstep = (x(2)-x(1)).*((1-prop)/2);
newx = [x(1)+xstep x(2)-xstep];

ystep = (y(2)-y(1)).*((1-prop)/2);
newy = [y(1)+ystep y(2)-ystep];

end