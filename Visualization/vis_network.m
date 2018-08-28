function h = vis_network(W, varargin)

% under construction

doweight = 0;
do3d = 0;
dodegree = 0; % change size of the node using degree
dogroup = 0;
g_cols = {};
wh_hl = [];

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            % functional commands
            case {'weighted'}
                doweight = 1;
            case {'3d'}
                do3d = 1;
                cl = varargin{i+1};
            case {'degree'}
                dodegree = 1;
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
[i,j] = find(W);

G = graph(i,j);
h = figure;
h_graph = plot(G,'Layout','force');
X = [h_graph.XData',h_graph.YData'];
close(h);

if do3d
    for i = 1:numel(cl), xyz(i,:) = cl(i).mm_center; end
end

[neg_i,neg_j] = find(W < 0);

for k = 1:numel(neg_i)
    h.edge_neg(k) = line([X(neg_i(k),1) X(neg_j(k),1)], [X(neg_i(k),2) X(neg_j(k),2)]);
    set(h.edge_neg(k), 'color', [0.1686    0.5137    0.7294], 'linewidth', 1);
    hold on;
end

[pos_i,pos_j] = find(W > 0);

for k = 1:numel(pos_i)
    h.edge_pos(k) = line([X(pos_i(k),1) X(pos_j(k),1)], [X(pos_i(k),2) X(pos_j(k),2)]);
    set(h.edge_pos(k), 'color', [0.8431    0.0980    0.1098], 'linewidth', 1);
    hold on;
end

h.node = scatter(X(:,1), X(:,2), 200, 'filled', 'MarkerFaceColor', [1 0 0], 'MarkerEdgeColor', 'w' , 'LineWidth', 1.5);

axis off;
set(gcf, 'position', [360   235   570   463], 'color', 'w');

end