function [h, stats] = ind_xylines(X, Y, varargin)

% Draw individual and group lines (from mean-2std to mean+2std)
%
% Usage:
% -------------------------------------------------------------------------
% [h, stats] = ind_xylines(X, Y, varargin)
%
% Inputs:
% -------------------------------------------------------------------------
% X     cell array: each cell for each subject contains any x-axis values
%       For example, it can be temperature (stimulus intensity)
% Y     cell array: each cell for each subject contains any y-axis values
%       that match with X. For example, it can be pain ratings. 
%
% Optional inputs: Enter keyword followed by variable with values
% 'ind_colors'     colors of individual lines   
% 'ind_linewidth'  linewidth of individual lines
%
% 'grp_colors'     color of a group line
% 'grp_linewidth'  linewidth of a group line
% 'nogrp'      no group line
% 'covs'           covariates 
% 'samefig'        draw the plot in the same fig
%
% Outputs:
% -------------------------------------------------------------------------
% h              graphic handles 
% stats          glmfit_multilevel results
%
% Examples: 
% -------------------------------------------------------------------------
% % data
% for i = 1:10, X{i} = rand(20,1); Y{i} = rand(20,1); end
% [h, stats] = ind_xylines(X, Y, 'ind_colors', [.5 .5 .5], 'ind_linewidth', 1, ...
%          'grp_colors', 'k', 'grp_linewidth', 2.5);
%
% savename = 'example_ind_xylines.pdf';
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
% Copyright (C) 2016  Wani Woo

% Programmers' notes:

ind_colors = [.5 .5 .5];
ind_linewidth = 1;
ind_linestyle = '-';

dogrp = true;
doindline = true;
grp_colors = [0 0 0];
grp_linewidth = 2;
grp_linestyle = '-';

subjn = numel(X);

X_cov = cell(subjn,1);
covariates= cell(subjn,1);
dosamefig = 0;

dogrpseband = false;
dogrpseline = false;

grp_selinestyle = '-';

for i = 1:subjn, covariates{i} = []; end

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            % functional commands
            case {'ind_colors'}
                ind_colors = varargin{i+1};
            case {'grp_colors'}
                grp_colors = varargin{i+1};
            case {'ind_linewidth'}
                ind_linewidth = varargin{i+1};
            case {'grp_linewidth'}
                grp_linewidth = varargin{i+1};
            case {'ind_linestyle'}
                ind_linestyle = varargin{i+1};
            case {'grp_linestyle'}
                grp_linestyle = varargin{i+1};
            case {'covs'}
                covariates = varargin{i+1};
            case {'samefig'}
                dosamefig = 1;
            case {'nogrp'}
                dogrp = false;
            case {'grpseband'}
                dogrpseband = true;
            case {'grpseline'}
                dogrpseline = true;
            case {'grp_selinestyle'}
                grp_selinestyle = varargin{i+1};
            case {'noind'}
                doindline = false;
        end
    end
end

for i = 1:subjn
    X_cov{i} = [X{i} covariates{i}]; 
end

if ~dosamefig
    create_figure('ind_xylines');
else
    hold on;
end

for i = 1:subjn
    meanX = nanmean(X_cov{i});
    stdX = nanstd(X_cov{i});
    
    if ~isempty(covariates{i})
        stdX(2:end) = 0;
    end
    
    stats.B{i} = glmfit(X_cov{i}, Y{i});
    
    newX = [meanX'-2*stdX' meanX'+2*stdX'];
    stats.newY{i} = glmval(stats.B{i}, newX, 'identity')';
    
    stats.newX{i} = newX(1,:);
    
    if ~isempty(covariates{i})
        covmean{i} = newX(2,:);
    end
    
    if doindline
        h.ind(i,1) = line(stats.newX{i}, stats.newY{i}, 'color', ind_colors, 'linewidth', ind_linewidth, 'linestyle', ind_linestyle);
    end
end

stats.glmfit_multilevel_stats = glmfit_multilevel(Y, X_cov, [], 'weighted', 'verbose');

stats.grpX = mean(cat(1,stats.newX{:}));

if isempty(covariates{1})
    stats.grpY = stats.glmfit_multilevel_stats.mean*[ones(1,2); stats.grpX];
else
    stats.grpY = stats.glmfit_multilevel_stats.mean*[ones(1,2); stats.grpX; mean(cat(1,covmean{:}))];
end

if dogrpseband && dogrpseline
    h.grp = plot_grperrorband(stats.glmfit_multilevel_stats.first_level.beta, ...
        stats.glmfit_multilevel_stats.W, stats.grpX, grp_colors, grp_selinestyle, 'fill', 'line');
elseif ~dogrpseband && dogrpseline
    h.grp = plot_grperrorband(stats.glmfit_multilevel_stats.first_level.beta, ...
        stats.glmfit_multilevel_stats.W, stats.grpX, grp_colors, grp_selinestyle, 'line');
elseif dogrpseband && ~dogrpseline
    h.grp = plot_grperrorband(stats.glmfit_multilevel_stats.first_level.beta, ...
        stats.glmfit_multilevel_stats.W, stats.grpX, grp_colors, grp_selinestyle, 'fill');
end

if dogrp
    h.grp = line(stats.grpX, stats.grpY, 'color', grp_colors, 'linewidth', grp_linewidth, 'linestyle', grp_linestyle);
end

end


function h = plot_grperrorband(betas, w, x, color, linestyle, varargin)

doline = false;
dofill = false;

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            % functional commands
            case {'line'}
                doline = true;
            case {'fill'}
                dofill = true;
        end
    end
end

% mean line
wmean = @(paths, w) ((w ./ sum(w))'*paths)';
% means = wmean(betas', w); % should be slope then intercept
% y = means(1) + means(2) * x; % mean(2) is slope, and mean(1) is intercept

% SE lines (bootstrapped; multilevel)
xx = linspace(x(1), x(2), 50);

boots = 5000;
means =  bootstrp(boots, wmean, betas', w);
yy = zeros(boots, 50);
for i = 1:boots
    yy(i,:) = means(i, 2) * xx + means(i, 1);
end

ub = prctile(yy, 95);
lb = prctile(yy, 5);

if doline    
    h.ub = plot(xx, ub, 'color', color, 'LineWidth', 1.5, 'linestyle', linestyle);
    h.lb = plot(xx, lb, 'color', color, 'LineWidth', 1.5, 'linestyle', linestyle);
end

if dofill
    h.fhan = fill([xx xx(end:-1:1)], [lb ub(end:-1:1)], color);
    h.fhan.LineStyle = 'none';
    set(h.fhan, 'FaceAlpha', .3)
end

% plot(x, y, '-', 'Color', 'k', 'LineWidth', 3);

% xb = mean(x); s = resid. std. stat.se(2) is se of slope
%se_mean = (stat.se(2).^2 * (xx-xb).^2 + s.^2/length(x)) .^ .5;

% drawnow

end

