function out = vis_corr(r, varargin)

% Draw a nice correlation figure
%
% Usage:
% -------------------------------------------------------------------------
% function h = vis_corr(r, varargin)
%
% Inputs:
% -------------------------------------------------------------------------
% r     square matrix (doesn't need to be a correlation)
%
% Optional inputs: Enter keyword followed by variable with values
% 'colormap'      colormap (default uses a colormap from colorbrewer2)
% 'clim'          two numbers [CLOW CHIGH] for color scaling of imagesc 
%                 (see >> help imagesc)
% 'colorbar'      show colorbar
%
% Outputs:
% -------------------------------------------------------------------------
% h              graphic handles 
%
% Examples: 
% -------------------------------------------------------------------------
% % data
% r = rand(30,30);
% h = vis_corr(r, 'clim', [0 1])
%
% savename = 'example_corr.png';
% 
% try
%     pagesetup(gcf);
%     saveas(gcf, savename);
% catch
%     pagesetup(gcf);
%     saveas(gcf, savename);    
% end
%
% -------------------------------------------------------------------------
% Copyright (C) 2015  Wani Woo

% Programmers' notes:

col_map = [0.3686    0.3098    0.6353
    0.1961    0.5333    0.7412
    0.4000    0.7608    0.6471
    0.6706    0.8667    0.6431
    0.9020    0.9608    0.5961
    1.0000    1.0000    1.0000
    0.9961    0.8784    0.5451
    0.9922    0.6824    0.3804
    0.9569    0.4275    0.2627
    0.8353    0.2431    0.3098
    0.6196    0.0039    0.2588];

col_map = interp1(1:size(col_map,1), col_map, linspace(1, size(col_map,1), 64));

do_display = 1;
docolorbar = 0;
dolines_1st = 1;
dolines_2nd = 1;
dolines_out = 1;
do_triangle = 0;

lcolor = [.7 .7 .7; .3 .3 .3; 0 0 0];
lwidth = [0.5 1 1.5]; % one five outline

dogroupsort = 0;
dogroupline = 0;

glcolor = [1 0 0];
glwidth = 2;
display_group_mean = 0;
display_group_sum = 0;
display_group_color = 0;
dosmooth = 0;
do_sig = false;
sig_lwidth = 3;
sig_col = [0.7529 0 0;0 0.4392 0.7529];

r_descript = 'input r';
triangle_col = 'k';
triangle_width = 5;

dogrouptick = 0;
tickstyle = 'edge';
tickcolor = 'k';
tickwidth = 2;
ticklength = 5;
tickoffset = 1;
same_fig = 0;
no_triangle_line = 0;

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            % functional commands
            case {'colormap', 'colors'}
                col_map = varargin{i+1};
            case {'colorbar'}
                docolorbar = 1;
            case {'nolines'}
                dolines_1st = 0;
                dolines_2nd = 0;
                dolines_out = 0;
            case {'nolines_1st'}
                dolines_1st = 0;
            case {'nolines_2nd'}
                dolines_2nd = 0;
            case {'nolines_out'}
                dolines_out = 0;
            case {'line_color'}
                lcolor = varargin{i+1};
            case {'line_width'}
                lwidth = varargin{i+1};
            case {'group'}
                dogroupsort = 1;
                dogroupline = 1;
                group = varargin{i+1};
            case {'group_linecolor'}
                glcolor = varargin{i+1};
            case {'group_linewidth'}
                glwidth = varargin{i+1};
            case {'group_mean'}
                display_group_mean = 1; 
                dogroupline = 0;
            case {'group_sum'}
                display_group_sum = 1;
                dogroupline = 0;
            case {'group_color'}
                display_group_color = 1; 
                group_color = varargin{i+1};
            case {'group_tick'}
                dogrouptick = 1;
            case {'group_tickstyle'}
                tickstyle = varargin{i+1};
                if strcmp(tickstyle, 'center')
                    tickcentering = -0.5;
                elseif strcmp(tickstyle, 'edge')
                    tickcentering = 0;
                end
            case {'group_tickwidth'}
                tickwidth = varargin{i+1};
            case {'group_ticklength'}
                ticklength = varargin{i+1};
            case {'group_tickoffset'}
                tickoffset = varargin{i+1};
            case {'smooth'}
                dosmooth = 1; % 2d gaussian smoothing (imgaussfilt) with sigma = 1
            case {'triangle'}
                do_triangle = 1;
            case {'triangle_color'}
                triangle_col = varargin{i+1};
            case {'triangle_width'}
                triangle_width = varargin{i+1};
            case {'sig'}
                do_sig = true;
                sig_matrix = varargin{i+1};
            case {'sig_linewidth'}
                sig_lwidth = varargin{i+1};
            case {'sig_color'}
                sig_col = varargin{i+1};
            case {'noplot', 'nodisplay'}
                do_display = 0;
            case {'same_fig'}
                same_fig = 1;
            case {'no_triangle_line'}
                no_triangle_line = 1;
        end
    end
end

% if there is a group variable
if dogroupsort
    [sorted_group, idx] = sort(group);
else
    idx = 1:size(r,1);
end

% display group mean instead of the raw r
if display_group_mean || display_group_sum
    sorted_r = r(idx,idx);
    u_group = unique(sorted_group);
    
    n_c = histc(sorted_group, u_group);
    n_c = [0; n_c];
    
    for i = 1:numel(u_group)
        for j = 1:numel(u_group)
            
            conn_mat = sorted_r((sum(n_c(1:i))+1):sum(n_c(1:(i+1))), (sum(n_c(1:j))+1):sum(n_c(1:(j+1))));
            
            if i == j
                non_zero_n = size(conn_mat,1).*(size(conn_mat,1)-1)./2;
                if display_group_mean
                    if non_zero_n ~= 0
                        net_mat(i,j) = sum(sum(triu(conn_mat,1)))./non_zero_n;
                    else
                        net_mat(i,j) = 0;
                    end
                elseif display_group_sum
                    net_mat(i,j) = sum(sum(triu(conn_mat,1)));
                end
            else
                if display_group_mean
                    net_mat(i,j) = sum(sum(conn_mat))./numel(conn_mat);
                elseif display_group_sum
                    net_mat(i,j) = sum(sum(conn_mat));
                end
            end            
        end
    end
    
    r = net_mat;       % new group-level r to display
    r_descript = 'group mean';
    idx = 1:size(r,1); % not to use group sorted idx
end

if do_display
    % default: upper 97.5%, lower 2.5%
    m = max(prctile(r(:), 97.5), abs(prctile(r(:), 2.5)));
    if m == 0 
        tr = r; 
        tr(r==0) = [];
        m = max(prctile(tr(:), 97.5), abs(prctile(tr(:), 2.5)));
    end
    clim = [-m m];
    clim_descript = 'default: upper 97.5%, lower 2.5%';
    
    for i = 1:length(varargin)
        if ischar(varargin{i})
            switch varargin{i}
                % functional commands
                case {'clim'}
                    clim = varargin{i+1};
                    clim_descript = 'optional input';
            end
        end
    end
    
    % check if r is a square matrix
    if size(r,1) == size(r,2)
        size_r = size(r,1);
    else
        error('r should be a square matrix.');
    end
    size_max = size_r+.5;
    
    % close all;
    if ~same_fig
        h = figure;
    else
        h = gcf;
    end
    
    % imagesc
    if dosmooth
        imagesc(imgaussfilt(r(idx,idx),1), [clim(1) clim(2)]);
    else
        imagesc(r(idx,idx), [clim(1) clim(2)]);
    end
    
    if docolorbar
        colorbar;
        set(gcf, 'color', 'w');
    else
        set(gcf, 'color', 'w');
    end
    
    axis off;
    colormap(col_map);
    
    set(gca, 'xlim', [-.5 size_r+1.5], 'ylim', [-.5 size_r+1.5])
    hold on;
    
    if dolines_1st
        % one
        for i = 1.5:1:(size_r-.5)
            line([i i], [.5 size_max], 'color', lcolor(1,:), 'linewidth', lwidth(1));
            line([.5 size_max], [i i], 'color', lcolor(1,:), 'linewidth', lwidth(1));
        end
    end
    
    if do_sig
        % one
        [b, a, c] = find(sig_matrix~=0);
        for i = 1:numel(a)
            x = [a(i)-.5 a(i)-.5 a(i)+.5 a(i)+.5 a(i)-.5];
            y = [b(i)-.5 b(i)+.5 b(i)+.5 b(i)-.5 b(i)-.5];
            if c(i) == 1
                line(x, y, 'linewidth', sig_lwidth, 'color', sig_col(1,:));
            else
                line(x, y, 'linewidth', sig_lwidth, 'color', sig_col(2,:));
            end
        end
    end
    
    if dolines_2nd
        % five
        for i = .5:5:(size_r+.5)
            line([i i], [.5 size_max], 'color', lcolor(2,:), 'linewidth', lwidth(2));
            line([.5 size_max], [i i], 'color', lcolor(2,:), 'linewidth', lwidth(2));
        end
    end
    
    if dolines_out
        for i = [.5 size_r+.5]
            line([i i], [.5 size_max], 'color', lcolor(3,:), 'linewidth', lwidth(3));
            line([.5 size_max], [i i], 'color', lcolor(3,:), 'linewidth', lwidth(3));
        end
    end
    
    if dogrouptick
        if ~display_group_mean && ~display_group_sum
            n_c = histc(sorted_group, unique(sorted_group));
        else
            n_c = ones(size(n_c));
        end
        
        set(gca, 'Clipping', 'off');
        
        for i = 1:numel(n_c)
            if i == numel(n_c) && strcmp(tickstyle, 'edge')
                break;
            end
            tick_base = [0 ticklength] + tickoffset;
            ytick_x1 = -tick_base + 0.5;
            ytick_x2 = tick_base + sum(n_c) +  0.5;
            ytick_y = repmat(sum(n_c(1:i)), 1, 2) + tickcentering + 0.5;
            line(ytick_x1, ytick_y, 'color', tickcolor, 'linewidth', tickwidth);
            if ~do_triangle; line(ytick_x2, ytick_y, 'color', tickcolor, 'linewidth', tickwidth); end
            xtick_x = ytick_y;
            xtick_y1 = ytick_x2;
            xtick_y2 = ytick_x1;
            line(xtick_x, xtick_y1, 'color', tickcolor, 'linewidth', tickwidth);
            if ~do_triangle; line(xtick_x, xtick_y2, 'color', tickcolor, 'linewidth', tickwidth); end
            
        end
    end
    
    if dogroupline
        n_c = histc(sorted_group, unique(sorted_group));
        n_c = [0; n_c];
        for i = 1:(numel(n_c)-1)
            
            xy1 = [sum(n_c(1:i)) sum(n_c(1:i))]+0.5;
            xy2 = [sum(n_c(1:(i+1))) sum(n_c(1:(i+1)))]+0.5;
            
            line([xy1(1) xy1(1)], [xy1(2), xy2(2)], 'color', glcolor, 'linewidth', glwidth);
            line([xy1(1) xy2(1)], [xy1(2) xy1(2)], 'color', glcolor, 'linewidth', glwidth);
            line([xy2(1) xy2(1)], [xy1(2) xy2(2)], 'color', glcolor, 'linewidth', glwidth);
            line([xy1(1) xy2(1)], [xy2(2) xy2(2)], 'color', glcolor, 'linewidth', glwidth);
            
        end
    end
    
    if display_group_color
        if ~display_group_mean && ~display_group_sum
            n_c = histc(sorted_group, unique(sorted_group));
            n_c = [0; n_c];
            set(gca, 'Clipping', 'off');
            for i = 1:(numel(n_c)-1)
                line_base = (ticklength + tickoffset) * 2;
                yline_x1 = repmat(-line_base + 0.5, 1, 2);
                yline_x2 = repmat(line_base + sum(n_c) + 0.5, 1, 2);
                if i == 1
                    yline_y = [sum(n_c(1:i)) sum(n_c(1:i+1))-1] + 0.5;
                elseif i == numel(n_c)-1
                    yline_y = [sum(n_c(1:i))+1 sum(n_c(1:i+1))] + 0.5;
                else
                    yline_y = [sum(n_c(1:i))+1 sum(n_c(1:i+1))-1] + 0.5;
                end
                line(yline_x1, yline_y, 'color', group_color(i,:), 'linewidth', 4);
                if ~do_triangle; line(yline_x2, yline_y, 'color', group_color(i,:), 'linewidth', 4); end
                xline_x = yline_y;
                xline_y1 = yline_x2;
                xline_y2 = yline_x1;
                line(xline_x, xline_y1, 'color', group_color(i,:), 'linewidth', 4);
                if ~do_triangle; line(xline_x, xline_y2, 'color', group_color(i,:), 'linewidth', 4); end
            end
        else
            for i = 1:size(group_color,1)
                scatter(-.5, i, 500, group_color(i,:), 'filled');
                scatter(i, size(group_color,1)+1.5, 500, group_color(i,:), 'filled');
            end
        end
    end
    
    if do_triangle
        patch([-1, size(r,1)+2, size(r,1)+2], [-1, -1, size(r,1)+2], 'w', 'edgecolor', 'w');
        if ~no_triangle_line
            patch([.5, .5, size(r,1)+.5], [.5, size(r,1)+.5, size(r,1)+.5], 'w', 'edgecolor', triangle_col, 'facealpha', 0, 'LineWidth', triangle_width);
        end
    end
    
    out.h = h;
    out.color_lim = clim;
    out.color_lim_descript = clim_descript;
end
out.r = r(idx,idx);
out.r_descript = r_descript;

if docolorbar
    set(gcf, 'position', [680   558   475   410]);
else
    set(gcf, 'position', [680   558   431   410]);
end


end






