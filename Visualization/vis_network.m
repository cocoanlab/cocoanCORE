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

A(isnan(W)) = 0;
A = sparse(W);

% A(isnan(A)) = 0;
% A = sparse(A);
X = fruchterman_reingold_force_directed_layout(A);

if do3d
    for i = 1:numel(cl), xyz(i,:) = cl(i).mm_center; end
end

[i,j] = find(W < 0);

for k = 1:numel(i)
    hl(k) = line([X(i(k),1) X(j(k),1)], [X(i(k),2) X(j(k),2)]);
    set(hl(k), 'color', [0    0.2333    0.8627], 'linewidth', 1);
    hold on;
end

[i,j] = find(W > 0);

for k = 1:numel(i)
    hl(k) = line([X(i(k),1) X(j(k),1)], [X(i(k),2) X(j(k),2)]);
    set(hl(k), 'color', [0.8608    0.2020         0], 'linewidth', 1);
    hold on;
end

h(i) = scatter(X(:,1), X(:,2), 200, 'filled', 'MarkerFaceColor', [1 0 0], 'MarkerEdgeColor', 'w' , 'LineWidth', 1.5);

axis off;
set(gcf, 'position', [360   235   570   463], 'color', 'w');

end