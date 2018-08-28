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
% 
% h = vis_network(W, 'weighted', 'degree', 'group', grouping, 'groupcolor', g_cols)
%
% ..
%    Programmers' notes:
%    Created 8/28/18 by wani woo
% ..
%


doweight = 0;
do3d = 0;
dodegree = 0; % change size of the node using degree
do_w_degree = 0;
do_absw_degree = 0;
dogroup = 0;
g_cols = {};
wh_hl = true(size(W,1),1);

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            % functional commands
            case {'weighted', 'weight'}
                doweight = 1;
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
            case {'groupcolor'}
                g_cols = varargin{i+1};
            case {'highlight'}
                wh_hl = varargin{i+1};
        end
    end
end

W(isnan(W)) = 0;
[i, j, s] = find(W);

new_s = (abs(s) - min(abs(s)))./(max(abs(s))-min(abs(s))); % 0-1 normalization
new_s = new_s.*4+1; % 1~4

G = graph(i,j);
hh = figure;
h_graph = plot(G,'Layout','force');
X = [h_graph.XData',h_graph.YData'];
close(hh);

if do3d
    for i = 1:numel(cl), xyz(i,:) = cl(i).mm_center; end
end

[neg_i,neg_j] = find(W < 0);

for k = 1:numel(neg_i)
    h.edge_neg(k) = line([X(neg_i(k),1) X(neg_j(k),1)], [X(neg_i(k),2) X(neg_j(k),2)]);
    if doweight
        set(h.edge_neg(k), 'color', [0.1686    0.5137    0.7294], 'linewidth', new_s(k));
    else
        set(h.edge_neg(k), 'color', [0.1686    0.5137    0.7294], 'linewidth', 1.5);
    end
    hold on;
end

[pos_i,pos_j] = find(W > 0);

for k = 1:numel(pos_i)
    h.edge_pos(k) = line([X(pos_i(k),1) X(pos_j(k),1)], [X(pos_i(k),2) X(pos_j(k),2)]);
    if doweight
        set(h.edge_pos(k), 'color', [0.8431    0.0980    0.1098], 'linewidth', new_s(k));
    else
        set(h.edge_pos(k), 'color', [0.8431    0.0980    0.1098], 'linewidth', 1);
    end
    hold on;
end

if dodegree
    d = sum(W~=0);
    d = ((d-min(d))./(max(d)-min(d))*4+.5)*100; % 50~450
elseif do_absw_degree
    d = sum(abs(W));
    d = ((d-min(d))./(max(d)-min(d))*4+.5)*100; % 50~450
elseif do_w_degree
    d = sum(W);
    d = ((d-min(d))./(max(d)-min(d))*4+.5)*100; % 50~450
else
    d = repmat(150, size(X,1), 1);
end

for node_i = 1:size(X,1)
    if wh_hl(node_i)
        if dogroup
            for i = 1:numel(grouping)
                h.node = scatter(X(node_i,1), X(node_i,2), d(node_i), 'filled', 'MarkerFaceColor', g_cols(grouping(node_i),:), 'MarkerEdgeColor', 'w' , 'LineWidth', 1.5);
            end
        else
            h.node = scatter(X(node_i,1), X(node_i,2), d(node_i), 'filled', 'MarkerFaceColor', [.5 .5 .5], 'MarkerEdgeColor', 'w' , 'LineWidth', 1.5);
        end
    else
        h.node = scatter(X(node_i,1), X(node_i,2), 50, 'filled', 'MarkerFaceColor', [.5 .5 .5], 'MarkerEdgeColor', 'w' , 'LineWidth', 1.5);
    end
end

axis off;
set(gcf, 'position', [360   235   570   463], 'color', 'w');

end