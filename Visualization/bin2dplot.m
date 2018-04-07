function [h, stats] = bin2dplot(X, Y, varargin)

% Draw a 2d plot with 2d error bars with binned data (multilevel)
%
% Usage:
% -------------------------------------------------------------------------
% [h, {stats}] = bin2dplot(X, Y, varargin)
%
% Inputs:
% -------------------------------------------------------------------------
% X     cell array: each cell for each subject contains any x-axis values
%       For example, it can be temperature (stimulus intensity)
% Y     cell array: each cell for each subject contains any y-axis values
%       that match with X. For example, it can be pain ratings. 
%
% Optional inputs: Enter keyword followed by variable with values
% 'ylim'     set up y-axis limit manually
% 'nbins'    the number of bins for x-axis (default: 4)
% 'colors'   colors of the dots; default: [0.8353 0.2431 0.3098];
% 'colors_refline'    colors for the reference line 
%                     default: [0.9569 0.4275 0.2627];
% 'refline'  draw reference line (regression); default: no refline
% 'reflinestyle'  linestyle for refline; default: '--'
% 'reflinewidth'  linewidth for refline; default: 1
%
% 'stats'    running glmfit_multilevel on X and Y
% 'covs'     covariates 
% 'resid'    draw plot with residualized X with covariates
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
% for i = 1:10, X{i} = rand(20,1); Y{i} = rand(20,1); end
% h = bin2dplot(X, Y, 'nbins', 10, 'refline');
%
% savename = 'example_bin2dplot.pdf';
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
nbins = 4;
reflinest = '--';
reflinew = 1;
do_resid = 0;
do_stat = 0;
dosameplot = 0;

colors = [0.8353    0.2431    0.3098];
colors_ref = [0.9569    0.4275    0.2627];

subjn = numel(X);

X_cov = cell(subjn,1);
covariates= cell(subjn,1);

Xbins = NaN(subjn, nbins);
Ybins = NaN(subjn, nbins);

for i = 1:subjn, covariates{i} = []; end

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
            case {'covs'}
                covariates = varargin{i+1};
            case {'resid'}
                do_resid = 1;
            case {'stats', 'stat'}
                do_stat = 1;
            case {'nbins'}
                nbins = varargin{i+1};
            case {'sameplot'}
                dosameplot = 1;

        end
    end
end

for i = 1:subjn
    X_cov{i} = [X{i} covariates{i}]; 
    
    if do_resid
        X{i} = resid(covariates{i}, X{i});
    end
end

if do_stat
    stats = glmfit_multilevel(Y, X_cov, [], 'verbose', 'weighted', 'boot', 'nresample', 10000);
end

clear Xbins Ybins;

if do_resid
    for i = 1:subjn
        Y{i} = resid(covariates{i}, Y{i});
    end
end

for i = 1:subjn
    
    % get binidx
    [~, sort_idx] = sort(X{i});
    binidx = zeros(numel(X{i}),1);
    if numel(unique(X{i})) == nbins
        u = unique(X{i});
        for j = 1:numel(u)
            binidx(X{i} == u(j)) = j;
        end
        
        for j = 1:nbins
            Xbins(i,j) = mean(X{i}(binidx==j));
            Ybins(i,j) = mean(Y{i}(binidx==j));
        end
        
    else
        algo = {'ceil', 'floor'};
        for j = 1:(nbins-1)
            algon = double(rand>.5)+1;
            eval(['tri_n = ' algo{algon} '(numel(X{i})./nbins);']);
            binidx(find(binidx==0, 1, 'first'):(find(binidx==0, 1, 'first')+tri_n-1)) = ...
                repmat(j, tri_n, 1);
        end
        binidx(binidx==0) = nbins;
        for j = 1:nbins
            Xbins(i,j) = mean(X{i}(sort_idx(binidx==j)));
            Ybins(i,j) = mean(Y{i}(sort_idx(binidx==j)));
        end
    end
    
end

x = nanmean(Xbins);
xe = ste(Xbins);

y = nanmean(Ybins);
ye = ste(Ybins);

if ~dosameplot
    h{1} = create_figure('2d_plot');
end
    
h{2} = scatter(x,y, 100, colors, 'filled');

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
    set(h{4}{i}(1), 'marker', '.', 'color', colors, 'markersize', 1);
    set(h{4}{i}(2), 'color', colors, 'linewidth', 2);
    set(h{4}{i}(3), 'color', colors, 'linewidth', 2);
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

end
