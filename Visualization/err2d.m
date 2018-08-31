function h = err2d(x, y, xe, ye, varargin)

% Draw a scatter plot with 2d error bars
%
% Usage:
% -------------------------------------------------------------------------
% function h = err2d(x, y, xe, ye, varargin)
%
% Inputs:
% -------------------------------------------------------------------------
% x     mean x
% y     mean y
%
% xe    standard error of x
% ye    standard error of y
%
% Optional inputs: Enter keyword followed by variable with values
% 'ylim'     set up y-axis limit manually
% 'colors'   colors of the dots; default: [0.8353 0.2431 0.3098];
% 'colors_refline'    colors for the reference line 
%                     default: [0.9569 0.4275 0.2627];
% 'refline'  draw reference line (regression); default: no refline
% 'reflinestyle'  linestyle for refline; default: '--'
% 'reflinewidth'  linewidth for refline; default: 1
%
% Outputs:
% -------------------------------------------------------------------------
% h              graphic handles 
%   h{1}         entire figure 
%   h{2}         dots    
%   h{3}         reference line (empty if there is no refline)
%   h{4}         2d error bars
%
% Examples: 
% -------------------------------------------------------------------------
% % data
% X = rand(20,20); Y = rand(20,20);
% x = mean(X);
% y = mean(Y);
% xe = ste(X);
% ye = ste(Y);
%
% h = err2d(x, y, xe, ye, 'refline')
%
% savename = 'example_2derrplot.pdf';
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

doman_ylim = 0;
dorefline = 0;
reflinest = '--';
reflinew = 1;
eb_width = 1.5;
m_size = 60;

colors = repmat([0.8353    0.2431    0.3098], numel(x), 1);
colors_ref = [0.9569    0.4275    0.2627];

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            % functional commands
            case {'ylim'}
                doman_ylim = 1;
                ylim = varargin{i+1};
            case {'colors', 'color'}
                colors = varargin{i+1};
            case {'colors_refline', 'colors_ref'}
                colors_ref = varargin{i+1};
            case {'refline'}
                dorefline = 1;
            case {'reflinestyle'}
                reflinest = varargin{i+1};
            case {'reflinewidth'}
                reflinew = varargin{i+1};
            case {'errorbar_width'}
                eb_width = varargin{i+1};
            case {'marker_size'}
                m_size = varargin{i+1};
        end
    end
end


h{1} = create_figure('2d_plot');
    
h{2} = scatter(x,y, m_size, colors, 'filled');

if dorefline
    h{3} = refline;
    set(h{3}, 'Color', colors_ref, 'linewidth', reflinew, 'linestyle', reflinest);
end
% 
% if doindiv
%     for i = 1:subjn
%         scatter(Xbins(i,:), Ybins(i,:));
%     end
% end

xmin = min(x-xe) - range(x)*.05;
xmax = max(x+xe) + range(x)*.05;
% ymin = min(y-ye) - range(y)*.05;
% ymax = max(y+ye) + range(y)*.05;

xmin2 = min(x-xe) - range(x)*.1;
xmax2 = max(x+xe) + range(x)*.1;
ymin2 = min(y-ye) - range(y)*.1;
ymax2 = max(y+ye) + range(y)*.1;


for i = 1:numel(x)
    h{4}{i} = ploterr(x(i),y(i),xe(i),ye(i));
    set(h{4}{i}(1), 'marker', '.', 'color', colors(i,:), 'markersize', 1);
    set(h{4}{i}(2), 'color', colors(i,:), 'linewidth', eb_width);
    set(h{4}{i}(3), 'color', colors(i,:), 'linewidth', eb_width);
    xdata = get(h{4}{i}(2), 'xData');
    xdata(4:5) = xdata(1:2); xdata(7:8) = xdata(1:2);
    set(h{4}{i}(2), 'xdata', xdata);
    
    ydata = get(h{4}{i}(3), 'yData');
    ydata(4:5) = ydata(1:2); ydata(7:8) = ydata(1:2);
    set(h{4}{i}(3), 'ydata', ydata);
    hold on;
end

if dorefline
    xdata = get(h{3}, 'xdata');
    ydata = get(h{3}, 'ydata');
    slope = (ydata(2)-ydata(1))./(xdata(2) - xdata(1));
    intercept = ydata(2) - xdata(2).*slope;
    set(h{3}, 'xdata', [xmin xmax], 'ydata', [xmin*slope+intercept xmax*slope+intercept])
end

if doman_ylim
    set(gca, 'xlim', [xmin2 xmax2], 'ylim', ylim, 'linewidth', 1.2, 'fontsize', 18);
else
    set(gca, 'xlim', [xmin2 xmax2], 'ylim', [ymin2 ymax2], 'linewidth', 1.2, 'fontsize', 18);
end

