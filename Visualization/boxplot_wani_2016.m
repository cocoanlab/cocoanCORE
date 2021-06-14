function [handles, dot_locs] = boxplot_wani_2016(x, varargin)

% Draw a box plot with some additional useful features. (work in the
% matlab version since 2016 maybe 2015 as well, but not tested). 
%
% Usage:
% -------------------------------------------------------------------------
% boxplot_wani_2016(x, varargin)
%
% Inputs:
% -------------------------------------------------------------------------
% x    A data matrix or cell. If there are multiple groups, you can insert
%      a column of NaNs (matrix) or insert a empty matrix (cell) between
%      the two groups.
%      Also if you have multiple columns that have different number of data
%      points (e.g., one has 50, and the other has 38), you can just use
%      cell-type input or manually fill in NaNs (e.g., 12 NaNs).
% 
% -------------------------------------------------------------------------
% Optional inputs: Enter keyword followed by variable with values
%
% ['color', 'colors', cols]     cols: matrix of 3 x n (n columns)
% ['boxlinewidth', line_box]    scalar, linewidth for the box
% ['boxlinecolor', boxcols]        boxcols: matrix of 3 x n (n columns)
% ['linewidth', line_etc]       scalar, linewidth for all the lines other 
%                               than box
% ['axislinewidth', line_axis]  scalar, linewidth for the axis
% ['fontsize', font_size]       scalar, font size for the axis
% ['dorefline', ref]            scalar, draw a reference line at y = ref
% ['reflinewidth', line_ref]    scalar, linewidth for the refline
% ['reflinestyle', reflinestyle]   '-', '--', ':', etc.
% ['reflinecolor', reflinecol]  reflinecol: matrix of 3 x n (n columns)
% ['dotcolor', dotcols]         outlier dot colors: dotcols = matrix of 3 x n (n columns)
% ['dotsize', dotsize]          outlier dot size: dotsize = scalar, size of the outlier dots
% ['mediancolor', mdcols]       mdcols: matrix of 3 x n (n columns)
% ['medianlinewidth', line_md]  scalar, linewidth for the median line
% ['samefig']                   
% ['violin']                    draw violin plot (only line)
% ['dots']                      data dots: show data dots
% ['dot_alpha', alpha]          data dots: alpha for data dots, default = .4
% ['dot_size', dot_size]        data dot size: dot_size = scalar, default = 40
% ['bw', bw]                    bandwidth of violin plot, default = [];
% ['data_dotcolor', color]      datapoint dot colors
%
% example:
% 
% x = rand(100,5);
% col =  [0.3765    0.2902    0.4824
%     0.2157    0.3765    0.5725
%     0.4667    0.5765    0.2353
%     0.8941    0.4235    0.0392
%     0.5843    0.2157    0.2078];
%
% % example 1
% boxplot_wani_2016(x, 'color', col, 'refline', 0.5);
% % example 2 : with thiner lines
% boxplot_wani_2016(x, 'color', col, 'refline', 0.5, 'linewidth', 2, 'boxlinewidth', 3);
% % example 3 : with violin plot and data dots
% boxplot_wani_2016(x, 'color', col, 'refline', 0.5, 'linewidth', 2, 'boxlinewidth', 3, 'violin', 'dots');
%
% savename = 'example_box.pdf';
% 
% try
%     pagesetup(gcf);
%     saveas(gcf, savename);
% catch
%     pagesetup(gcf);
%     saveas(gcf, savename);
% end

if iscell(x)
    coln = numel(x);
    eachn = cellfun(@numel, x);
    maxn = max(eachn);
    x_new = NaN(maxn, coln);
    for i = 1:coln
        x_new(1:eachn(i),i) = x{i};
    end
    x = x_new;
end

coln = size(x,2);
colud = repmat([0 0 0], coln*2, 1); % default color = black
colud2 = repmat([.5 .5 .5], coln*2, 1); % default back color = black
colud3 = flipud(colud2); 

dorefline = 0;
line_etc = 2;
line_md = 2;
line_box = 3;
line_ref = 1.5;
line_axis = 1.5;
font_size = 25;
dotcolor = [0.7608    0.3020         0];
dotsize = 12;
reflinestyle = '--';
reflinecol = [.7 .7 .7];
boxlinestyle = {'-'};
mdcol = [0.7608 0.3020 0];
usesamefig = false;
facealpha = 1;
doviolin = 0;
violinalpha = 0;
dodots = 0;
dot_size = 40;
dot_alpha = .4;
bw = [];
data_dotcolor = [];
use_onedotcolor = false;
do_compact = false;
do_box = true;
do_boxtrans = false;

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            % functional commands
            case {'color', 'colors'}
                colud2 = flipud(varargin{i+1}); 
                colud = repmat(colud2, 100, 1);
                colud3 = flipud(colud2); 
            case {'box_trans'}
                do_boxtrans = 1;
                facealpha = varargin{i+1};
            case {'refline'}
                dorefline = 1;
                ref = varargin{i+1};
            case {'linewidth'}
                line_etc = varargin{i+1};
            case {'boxlinewidth'}
                line_box  = varargin{i+1};
            case {'boxlinecolor'}
               colud = repmat(flipud(varargin{i+1}), 100, 1);
            case {'reflinewidth'}
                line_ref  = varargin{i+1};
            case {'reflinestyle'}
                reflinestyle  = varargin{i+1};
            case {'reflinecolor'}
                reflinecol  = varargin{i+1};
            case {'axislinewidth'}
                line_axis  = varargin{i+1};
            case {'fontsize'}
                font_size = varargin{i+1};
            case {'dotcolor'}
                dotcolor = varargin{i+1};
            case {'dotsize'}
                dotsize = varargin{i+1};
            case {'mediancolor'}
                mdcol = varargin{i+1};
            case {'medianlinewidth'}
                line_md = varargin{i+1};
            case {'samefig'}
                usesamefig = true;
            case {'violin'}
                doviolin = 1;
            case {'dots'}
                dodots = 1;
            case {'dot_alpha'}
                dot_alpha = varargin{i+1};
            case {'dot_size'}
                dot_size = varargin{i+1};
            case {'bw'}
                bw = varargin{i+1};
            case {'data_dotcolor'}
                data_dotcolor = varargin{i+1};
                use_onedotcolor = true;
            case {'compact'}
                do_compact = true;
                mdcol = 'none';
            case {'nobox'}
                do_box = false;
            case {'violin_alpha'}
                violinalpha = varargin{i+1};
        end
    end
end

if ~usesamefig
    create_figure('box_plot');
end

if do_box
    if ~do_compact
        boxplot(x); % using boxplot default
    else
        boxplot(x, 'PlotStyle','compact'); % using boxplot default
    end
    
    % h = get(get(gca, 'children'), 'children');
    h = findobj(gca, 'Tag', 'Box');
    
    k=0;
    % if iscell(h), h = h{1}; end
    for i = 1:length(h)
        if isequal(get(h(i), 'color'), [0 0 1])
            k = k+1;
            patchdata.x{k} = get(h(i), 'xdata');
            patchdata.y{k} = get(h(i), 'ydata');
        end
    end
    
    handles.boxplot = h;
    
    clf;
    
    for j = 1:2 % just twice
        for i = 1:numel(patchdata.x)
            hh = patch(patchdata.x{i}, patchdata.y{i}, colud2(i,:), 'EdgeColor', colud2(i,:), 'FaceAlpha', facealpha);
        end
        
        hold on;
        if ~do_compact
            boxplot(x); % using boxplot default
        else
            boxplot(x, 'PlotStyle','compact'); % using boxplot default
        end
        
        set(gca, 'fontSize', font_size, 'lineWidth', line_axis, 'xlim', ...
            [0.2 coln+.8], 'XTickLabelMode', 'auto', 'XTickMode', 'auto');
        set(gcf, 'position', [50   159   105*coln   291]);
        h = findobj(gca, 'Type', 'Line');
        %h = h{1};
        k=0;
        for i = 1:length(h)
            set(h(i), 'lineWidth', line_etc)
            if strcmp(h(i).Tag, 'Box')
                k = k+1;
                if numel(boxlinestyle) == 1
                    % set(h(i), 'color', colud(k,:), 'linewidth', line_box, 'linestyle', boxlinestyle{1});
                    if ~do_boxtrans
                        set(h(i), 'color', colud(k,:), 'linewidth', line_box, 'linestyle', boxlinestyle{1});
                    else
                        set(h(i), 'color', 'w', 'linewidth', 0.0001, 'linestyle', boxlinestyle{1});
                    end
                else
                    % set(h(i), 'color', colud(k,:), 'linewidth', line_box, 'linestyle', boxlinestyle{k});
                    if ~do_boxtrans
                        set(h(i), 'color', colud(k,:), 'linewidth', line_box, 'linestyle', boxlinestyle{1});
                    else
                        set(h(i), 'color', 'w', 'linewidth', 0.0001, 'linestyle', boxlinestyle{1});
                    end
                end
            end
        end
        
        for i = 1:length(h) %(3*coln)
            % set(h(i), 'lineWidth', line_etc)
            
            if strcmp(h(i).Tag, 'Outliers')
                set(h(i), 'marker', '.', 'markerSize', dotsize, 'MarkerEdgeColor', dotcolor)
            end
        end
        
        if j == 1
            if dorefline
                for ii = 1:numel(ref)
                    l = refline([0 ref(ii)]);
                    set(l, 'color', reflinecol, 'linestyle', reflinestyle, 'linewidth', line_ref);
                end
            end
        end
    end
    
    handles.boxplot_others = h;
end

if doviolin
    x_cell = enforce_cell_array(x);
    hold on;
    if ~isempty(bw)
%         violinplot(x_cell, 'facecolor', colud3, 'edgecolor', colud3, ...
%             'x', 1:numel(x_cell), 'mc', [0.3686    0.3098    0.6353], 'medc', mdcol, 'nopoints', 'facealpha', 0, 'linewidth', 1.5, 'bw', bw);
        h = violinplot(x_cell, 'facecolor', colud3, 'edgecolor', colud3, ...
            'x', 1:numel(x_cell), 'mc', 'none', 'medc', mdcol, 'nopoints', 'facealpha', violinalpha, 'linewidth', 1.5, 'bw', bw);
    else
%         violinplot(x_cell, 'facecolor', colud3, 'edgecolor', colud3, ...
%             'x', 1:numel(x_cell), 'mc', [0.3686    0.3098    0.6353], 'medc', mdcol, 'nopoints', 'facealpha', 0, 'linewidth', 1.5);
        h = violinplot(x_cell, 'facecolor', colud3, 'edgecolor', colud3, ...
            'x', 1:numel(x_cell), 'mc', 'none', 'medc', mdcol, 'nopoints', 'facealpha', violinalpha, 'linewidth', 1.5);
    end
    legend off
    
    handles.violin = h;
end

x_cell = enforce_cell_array(x);
xvalues = get_violin_points(1:size(x,2), x);

dot_locs.x = xvalues;
dot_locs.y = x_cell;

if dodots
    for i = 1:numel(xvalues)
        if ~use_onedotcolor
            data_dotcolor = colud3(i,:);
        end
        h = scatter(xvalues{i}, x_cell{i}, dot_size, data_dotcolor, 'filled', 'MarkerFaceAlpha', dot_alpha);
    end
    handles.dots = h;
end

hh = findobj('Tag', 'Median');

for i = 1:numel(hh)
    set(hh(i), 'color', mdcol, 'linewidth', line_md, 'linestyle', '-');
end

handles.median = hh;

set(gca, 'xtick', find(sum(isnan(x))~=size(x,1)), 'xticklabel', ' ',...
    'box', 'off',  'TickLength', [.015 .015], 'TickDir', 'out');
set(gca, 'fontSize', font_size, 'lineWidth', line_axis, 'xlim', ...
    [0.2 coln+.8], 'XTickMode', 'auto');
set(gcf, 'position', [50   159   105*coln   291]);

end


function Y = enforce_cell_array(Y)

k = size(Y, 2);

if ~iscell(Y)
    for i = 1:k
        Ytmp{i} = Y(:, i);
        Ytmp{i}(isnan(Ytmp{i})) = [];
    end
    Y = Ytmp;
end
end % function


function xvalues = get_violin_points(x, Y)
% x = vector of x positions for each "column"
% Y = cell array of input data, one cell per "column"
%
% from violinplot.m

nbins = 10;

k = size(Y, 2);

xvalues = cell(1, k);

% Enforce cell, no NaNs
% ------------------------------------------------
Y = enforce_cell_array(Y);

% calculate  density values
% ------------------------------------------------
for i = 1:k
    
    [f, u, bb]=ksdensity(Y{i});
    
    f=f/max(f)*0.3; %normalize
    F(:,i) = f;
    U(:,i) = u;
    
end

% get x positions
% ------------------------------------------------

for i = 1:k
    
    myx = x(i);     % x-value for this bar in plot
    
    myU = U(:, i);  % x-values of ksdensity output
    myF = F(:, i);  % y-values (density) of ksdensity output
    
    myY = Y{i};     % data points
    mybins = linspace(min(myY), max(myY), nbins);
    
    % starting and ending values
    st = [-Inf mybins(1:end-1)];
    en = [mybins];
    
    for j = 1:nbins
        % define points within a bin or 'slab'
        
        whpoints = myY > st(j) & myY <= en(j);
        
        if sum(whpoints) == 0, continue, end
        
        whu = myU > st(j) & myU <= en(j);
        
        mylimit(j) = mean(myF(whu),'omitnan')/1.5;  % average density for this 'slab' of points
        
        my_xvals = linspace(myx - mylimit(j), myx + mylimit(j), sum(whpoints))';
        
        if length(my_xvals) == 1, my_xvals = myx;  end
        
        ylocs = myY(whpoints);
        
        xlocs = my_xvals(1:length(ylocs));
        
        % save in original point list
        xvalues{i}(whpoints, 1) = xlocs;
        
    end % slab
    
end % column

end % function
