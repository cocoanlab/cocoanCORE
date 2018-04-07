function boxplot_wani(x, varargin)

% Draw a box plot with some additional useful features. (work in the
% matlab version up to 2014a). 
%
% Usage:
% -------------------------------------------------------------------------
% boxplot_wani(x, varargin)
%
% Inputs:
% -------------------------------------------------------------------------
% x    A data matrix. If there are multiple groups, you can insert a column
%      of NaNs between two groups. Also if you have multiple columns that
%      have different number of data points (e.g., one has 50, and the 
%      other has 38), you can just fill in NaNs (e.g., 12 NaNs). 
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
% ['dotcolor', dotcols]         dotcols: matrix of 3 x n (n columns)
% ['dotsize', dotsize]          scalar, size of the outlier dots
% ['mediancolor', mdcols]       mdcols: matrix of 3 x n (n columns)
% ['medianlinewidth', line_md]  scalar, linewidth for the median line
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
% boxplot_wani(x, 'color', col, 'refline', 0.5);
% % example 2 : with thiner lines
% boxplot_wani(x, 'color', col, 'refline', 0.5, 'linewidth', 2, 'boxlinewidth', 3);
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

coln = size(x,2);
colud = repmat([0 0 0], coln, 1); % default color = black
colud2 = repmat([.9 .9 .9], coln, 1); % default back color = black

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

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            % functional commands
            case {'color', 'colors'}
                colud2 = flipud(varargin{i+1}); 
                colud = colud2;
            case {'refline'}
                dorefline = 1;
                ref = varargin{i+1};
            case {'linewidth'}
                line_etc = varargin{i+1};
            case {'boxlinewidth'}
                line_box  = varargin{i+1};
            case {'boxlinecolor'}
                colud  = flipud(varargin{i+1});
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
        end
    end
end


create_figure('box_plot');
boxplot(x); % using boxplot default

h = get(get(gca, 'children'), 'children');

k=0;
for i = (3*coln+1):length(h)
    if isequal(get(h(i), 'color'), [0 0 1])
        k = k+1;
        patchdata.x{k} = get(h(i), 'xdata');
        patchdata.y{k} = get(h(i), 'ydata');
    end
end

clf;

for j = 1:2 % just twice
    for i = 1:numel(patchdata.x)
        patch(patchdata.x{i}, patchdata.y{i}, colud2(i,:), 'EdgeColor', colud2(i,:));
    end
    
    hold on;
    boxplot(x); % because of the patch, do this again
    set(gca, 'fontSize', font_size, 'lineWidth', line_axis, 'xlim', ...
        [0.2 coln+.8], 'XTickLabelMode', 'auto', 'XTickMode', 'auto');
    set(gcf, 'position', [50   159   105*coln   291]);
    h = get(get(gca, 'children'), 'children');
    h = h{1};
    k=0;
    for i = (3*coln+1):length(h)
        set(h(i), 'lineWidth', line_etc)
        if isequal(get(h(i), 'color'), [0 0 1])
            k = k+1;
            if numel(boxlinestyle) == 1
                set(h(i), 'color', colud(k,:), 'linewidth', line_box, 'linestyle', boxlinestyle{1});
            else
                set(h(i), 'color', colud(k,:), 'linewidth', line_box, 'linestyle', boxlinestyle{k});
            end
        end
    end
    
    for i = 1:(3*coln)
        set(h(i), 'lineWidth', line_etc)
        if i > coln && i <= (coln*2)
            set(h(i), 'marker', '.', 'markerSize', dotsize, 'MarkerEdgeColor', dotcolor)
        end
        
        if isequal(get(h(i), 'color'), [1 0 0])
            set(h(i), 'color', mdcol, 'linewidth', line_md, 'linestyle', '-');
        end
    end
    
    if j == 1
        if dorefline
            l = refline([0 ref]);
            set(l, 'color', reflinecol, 'linestyle', reflinestyle, 'linewidth', line_ref);
        end
    end
end

set(gca, 'xtick', find(sum(isnan(x))~=size(x,1)), 'xticklabel', ' ',...
    'box', 'off',  'TickLength', [.015 .015], 'TickDir', 'out');

end

