function out = glass_brain_network(r, varargin)

% This function draws glass brain with nodes and edges for the whole brain.
%
% :Usage:
% ::
%
%    h = glass_brain_network(r, varargin)
%
% :Features:
%
%  - can draw glass brain. You can adjust cortex and cerebellum's alpha level
%  - can draw nodes as 3d sphere or just circles
%  - can display group membership using different colors and radii.
%  - can draw weighted edges. 
%  - You can choose to draw positive and negative edges only.
%
% :Inputs:
%
%   **r:**
%        region object 
%
% :Optional Inputs:
%
%   **sphere**
%        draw nodes as 3d sphere. No need more inputs. 
%   **radius**
%        Node radius, default is 2. If you want to assign different radius
%        to different groups, you can make a vector of radii [#group x 1]. 
%           e.g., 'radius', 3
%   **colors**
%        If you want to use one color for all regions, you can put one
%        1 x 3 color vector. If you want to use different colors for
%        different groups, you can make a vector of color [group# x 3]. 
%           e.g., 'colors', [.6 .3 .2]
%   **group**
%        Group membership for the regions. The length of the group vector 
%        should be same with the number of regions.
%           e.g., 'group', group_indices
%   **edge_weights**
%        edge weights, if you want to draw edges. This can be multiple weights
%        in a cell array
%           e.g., 'edge_weights', w
%   **edge_alpha**
%        edge alpha. This can be cell array.
%           e.g., 'edge_alpha', .3
%   **hl_node_edge**
%        highlight the only nodes that have edges. This needs three more 
%        inputs. The first one is the radius for highlighted regions. The
%        second one is the radius for un-highlighted regions. The third one
%        is the color for un-highlighted regions. 
%           e.g., 'hl_node_edge', 2.5, 1, [.7 .7 .7]
%   **positive_only**
%        If you want to draw positive edges only, you can use this option.
%        No need more input.
%   **negative_only**
%        If you want to draw negative edges only, you can use this option.
%        No need more input.
%   **pos_edge_color**
%        If you want to use a specific color for positive edges, use this
%        option. This can be cell array with multiple colors for different 
%        weight sets. 
%           e.g., 'pos_edge_color', [1 0 0]
%   **neg_edge_color**
%        If you want to use a specific color for negative edges, use this
%        option. This can be cell array with multiple colors for different 
%        weight sets. 
%           e.g., 'neg_edge_color', [0 0 1]
%   **cortex_alpha**
%        You can adjust the alpha value for the cortex using this option.
%           e.g., 'cortex_alpha', .03 (default)
%   **cerebellum_alpha**
%        You can adjust the alpha value for the cortex using this option.
%           e.g., 'cerebellum_alpha', .1 (default)
%
% :Example:
%    
%    r = region(FanBS, 'unique_mask_values');
%    h = glass_brain_network(r, 'sphere', 'group', cluster_idx, 'colors', cols);
% 
%    h = glass_brain_network(r, 'group', cluster_idx, 'colors', cols, 'edge_weights', w, 'hl_node_edge', 2.5, 1, [.7 .7 .7], 'positive_only');

% default
dosphere = false;
group_radius = 2;
group_cols = [.7 .7 .7];
color_unhl = [.3 .3 .3];

do_pos = true;
do_neg = true;
%w = zeros(n_region, n_region);
do_highlight = false;

cortex_alpha = 0.03;
cerebellum_alpha = .1;
normfactor_input = [];

pos_edge_color = [215,25,28]./255;
neg_edge_color = [43,131,186]./255;
edge_alpha = 1;
% pos_edge_color = 'r';
% neg_edge_color = 'b';
cortex = which('surf_BrainMesh_ICBM152.mat');
do_cerebellum = true;
do_ellipsoid = false;
node_alpha = 1;

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            % functional commands
            case {'sphere'}
                dosphere = true;
            case {'radius'}
                group_radius = varargin{i+1};
            case {'colors'}
                group_cols = varargin{i+1};
            case {'edge_weights'}
                w_cell = varargin{i+1};
            case {'hl_node_edge'}
                do_highlight = true;
                radius_hl = varargin{i+1};
                radius_unhl = varargin{i+2};
                color_unhl = varargin{i+3};
            case {'positive_only'}
                do_pos = true;
                do_neg = false;
            case {'negative_only'}
                do_pos = false;
                do_neg = true;
            case {'cortex_alpha'}
                cortex_alpha = varargin{i+1};
            case {'cerebellum_alpha'}
                cerebellum_alpha = varargin{i+1};
            case {'pos_edge_color'}
                pos_edge_color = varargin{i+1};
            case {'neg_edge_color'}
                neg_edge_color = varargin{i+1};
            case {'edge_alpha'}
                edge_alpha = varargin{i+1};
            case {'node_alpha'}
                node_alpha = varargin{i+1};
            case {'norm_factor'}
                normfactor_input = varargin{i+1};
            case {'center'}
                centers = varargin{i+1};
                n_region = size(centers,1);
            case {'inflated'}
                cortex = which('surf_BrainMesh_ICBM152_smoothed.mat');
            case {'left'}
                cortex = which('surf_BrainMesh_ICBM152Left.mat');
            case {'right'}
                cortex = which('surf_BrainMesh_ICBM152Right.mat');
            case {'left_inflated'}
                cortex = which('surf_BrainMesh_ICBM152Left_smoothed.mat');
            case {'right_inflated'}
                cortex = which('surf_BrainMesh_ICBM152Right_smoothed.mat');
            case {'very_inflated'}
                cortex = which('surf_workbench_very_inflated_32k.mat');
            case {'left_very_inflated'}
                cortex = which('surf_workbench_very_inflated_32k_Left.mat');
            case {'right_very_inflated'}
                cortex = which('surf_workbench_very_inflated_32k_Right.mat');
            case {'nocerebellum'}
                do_cerebellum = false;
            case {'ellipsoid'}
                do_ellipsoid = true;
                
        end
    end
end

if ~exist('n_region', 'var')
    n_region = numel(r);
end
group = ones(n_region,1);

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            case {'group'}
                group = varargin{i+1};
        end
    end
end

% colors for each node
cols = group_cols(group,:);
try
    radius = group_radius(group,:);
catch
    if size(group_radius,1) == 1
        radius = repmat(group_radius, numel(group),1);
    else
        error('Check whether the number of radius is same with the number of groups.');
    end
end

if do_highlight
    if iscell(w_cell)
        w  = zeros(size(w_cell{1}));
        for i = 1:numel(w_cell)
            w = w + w_cell{i};
        end
    else
        w = w_cell;
    end
    if issymmetric(w) && do_pos && ~do_neg
        w_temp = w .* double(w > 0);
    elseif issymmetric(w) && ~do_pos && do_neg
        w_temp = w .* double(w < 0);
    elseif issymmetric(w) && do_pos && do_neg
        w_temp = w;
    else
        error('Please check your w. It''s not symmetric.');
    end
    
    regions_hl = find(any(w_temp))';
    regions_unhl = find(~any(w_temp))';
    radius(regions_hl) = radius_hl;
    radius(regions_unhl) = radius_unhl;
    cols(regions_unhl,:) = repmat(color_unhl, numel(regions_unhl), 1);
end

% cerebellum
% if do_cerebellum, cluster_surf(region, which('surf_BrainMesh_Cerebellum_by_SLF.mat'), 2); end
draw_surf(which('surf_BrainMesh_Cerebellum_by_SLF.mat'));
% whole-brain 
% cluster_surf(region, cortex, 2);
draw_surf(cortex);

% make glass brains
h = get(gca, 'children');
for i = 2:3
    if i == 2
        set(h(i), 'FaceAlpha', cortex_alpha);
    else
        if do_cerebellum, set(h(i), 'FaceAlpha', cerebellum_alpha); end
    end
end

out.h1 = h;

axis vis3d;

if ~exist('centers', 'var')
    centers = cat(1,r.mm_center);
end

% draw all nodes

[x,y,z] = sphere;

hold on;

for i = 1:size(centers,1)
    
    if do_ellipsoid
        h = surf(x *radius(i,1)+centers(i,1),y*radius(i,2)+centers(i,2),z*radius(i,3)+centers(i,3));
    else
        h = surf(x *radius(i)+centers(i,1),y*radius(i)+centers(i,2),z*radius(i)+centers(i,3));
    end
    h.FaceColor = cols(i,:);
    h.EdgeAlpha = 0;
    h.FaceAlpha = node_alpha;
    
    if dosphere
        h.FaceLighting = 'gouraud';
    else
        h.FaceLighting = 'none';
    end
    
    h.AmbientStrength = .4;
    h.DiffuseStrength = 1;
    h.SpecularStrength = .1;
    
end

out.h2 = h;

%% draw edges

if exist('w_cell', 'var') 
    if iscell(w_cell)
        for i = 1:numel(w_cell)
            w = w_cell{i};
            
            if iscell(pos_edge_color), pos_col = pos_edge_color{i}; else, pos_col = pos_edge_color; end
            if iscell(neg_edge_color), neg_col = neg_edge_color{i}; else, neg_col = neg_edge_color; end
            if iscell(edge_alpha), edge_alpha2 = edge_alpha{i}; else, edge_alpha2 = edge_alpha; end
            
            normfactor{i} = draw_edges(w, pos_col, neg_col, edge_alpha2, normfactor_input, do_pos, do_neg, centers);
            hold on;
        end
    else
        w = w_cell;
        w(logical(eye(size(w,1)))) = 0; % remove diagonal
        normfactor = draw_edges(w, pos_edge_color, neg_edge_color, edge_alpha, normfactor_input, do_pos, do_neg, centers);
    end
    
    out.normfactor = normfactor;
end

end

% --------------- SUBFUNCTIONS --------------- 

function normfactor = draw_edges(w, pos_edge_color, neg_edge_color, edge_alpha, normfactor_input, do_pos, do_neg, centers)


[a,b] = find(w~=0);

pos_ab = [];
neg_ab = [];

for i = 1:numel(a)
    if w(a(i), b(i)) > 0
        pos_ab = [pos_ab; [a(i) b(i) w(a(i), b(i))]];
    else
        neg_ab = [neg_ab; [a(i) b(i) w(a(i), b(i))]];
    end
end

% l2norm

l2normfun = @(x) sum(x .^ 2) .^ .5;
all_ab = [pos_ab; neg_ab];

if do_pos && ~isempty(pos_ab)
    
    normfactor = l2normfun(all_ab(:,3))./sqrt(size(all_ab,1));
    if ~isempty(normfactor_input), normfactor = normfactor_input; end
    
    % norm_lw = ((pos_ab(:,3)-min(pos_ab(:,3)))/max(pos_ab(:,3))+1).^3;
    norm_lw = (pos_ab(:,3)./normfactor).*3;
    
    for i = 1:size(pos_ab,1)
        
        hold on;
        xyz = [centers(pos_ab(i,1),:); centers(pos_ab(i,2),:)];
        line(xyz(:,1),xyz(:,2),xyz(:,3), 'color', [pos_edge_color edge_alpha], 'linewidth', norm_lw(i));
        
    end
    
end

if do_neg && ~isempty(neg_ab)
    
    normfactor = l2normfun(all_ab(:,3))./sqrt(size(all_ab,1));
    if ~isempty(normfactor_input), normfactor = normfactor_input; end
    
    % norm_lw = ((-neg_ab(:,3)-min(-neg_ab(:,3)))/max(-neg_ab(:,3))+1).^3;
    norm_lw = (neg_ab(:,3)./normfactor).*3;
    
    for i = 1:size(neg_ab,1)
        
        hold on;
        xyz = [centers(neg_ab(i,1),:); centers(neg_ab(i,2),:)];
        %         line(xyz(:,1),xyz(:,2),xyz(:,3), 'color', 'b', 'linewidth', -neg_ab(i,3)*30000);
        line(xyz(:,1),xyz(:,2),xyz(:,3), 'color', [neg_edge_color edge_alpha], 'linewidth', -norm_lw(i));
        
    end
end
end

function draw_surf(P)


actcolor = [1 1 1];
basecolor = [.5 .5 .5]
mind = 2;
ovlcolor = [0 1 1];

load(P);

%%figure
p = patch('Faces',faces,'Vertices',vertices,'FaceColor',[.5 .5 .5], ...
    'EdgeColor','none','SpecularStrength',.2,'FaceAlpha',1,'SpecularExponent',200);
lighting gouraud; camlight right
axis image;
lightRestoreSingle(gca);
%myLight = camlight(0,0);set(myLight,'Tag','myLight');
%set(gcf, 'WindowButtonUpFcn', 'lightFollowView');lightFollowView

view(135, 30);
drawnow

[~, ~] = getVertexColors([], p, actcolor, basecolor, mind, 'ovlcolor', ovlcolor);


% -------------------------------------------------------------------------
% * run color change
% -------------------------------------------------------------------------
fprintf(' Running color change.\n');
axis off;
set(gcf, 'color', 'w');

end