function results = principal_gradients_plot(thresh_img, varargin)
% 
% =========================================================================
%
%                   STILL WORKING (Be careful).                             
%                          
% =========================================================================
% Draw line plot using principal_gradients_axis map. The map was made using
% fifty-six resting-state connectivity (3mm whole-brain voxel-level
% connectivity).
%
% :Usage:
% ::
%
%    results = principal_gradients_plot(thresh_img, varargin)
% 
% :Output: 
%  
%   **results:**
%       results.bins_Percentage: 
%       results.bins_nVox:
%:
% :Inputs:
%
%   **thresh_img:**
%        thresholded image files (e.g., thresh_img = 'thresh_p_05.nii';)
%
%
% :Optional Inputs: Enter keyword followed by variable with values
% 
%   **'other_comp':**
%        use other gradient component (default: 1 ; first gradients) 
%
%   **'noplots':**
%        do not show sepplot (defualt: true). If ture, the variable
%        'results.bins_results_Percent' will use for drawing plot.
%
%   **'numbins:**
%        change nubmer of bins (default: 20).
%
% :Examples: 
%
%   % ** Use first gradient map
%   results = principal_gradients_plot(thresh_img);
%   pagesetup(gcf);
%   savenames = 'temp_map.pdf'
%   saveas(gcf, savenames);
%
%   % **Use other component with different number of bins
%   numComp = 2; % 1:10
%   results = principal_gradients_plot(thresh_img,'other_comp',numComp,'numBins',10); 
%

%% Parse various argumetns 
doplot = true;
pgmapdir = which('ten_prinicpal_gradients_volumn_DE.nii'); % In cocoanCORE/Canonical_brains
comp_num = 1; % first gradient axis (transmoal - unimodal)
numBins  = 20;
for i = 1:length(varargin)
    if ischar(varargin{i})
        switch lower(varargin{i})
            % functional commands
%             case {'gradientmap'}
%                 pgmapdir = varargin{i+1};
            case {'noplots'}
                doplot = false;
            case {'other_comp'}
                comp_num = varargin{i+1};
            case {'numbins'}
                numBins = varargin{i+1};
            otherwise
                error('Unknown arguments')                 
        end
    end
end
%% Load gradients maps 
temp_pgmap = fmri_data(pgmapdir,which('gray_matter_mask.nii')); 
pgmap = temp_pgmap.get_wh_image(comp_num); 
pgmap.dat = pgmap.dat - min(pgmap.dat); % make data positive
%% compare voxel space and size 
isdiff = compare_space(pgmap, thresh_img);

if isdiff == 1 || isdiff == 2 % diff space, not just diff voxels
    % == 3 is ok, diff non-empty voxels
    
    % Both work, but resample_space does not require going back to original
    % images on disk.    
    mask = resample_space(pgmap, thresh_img,'nearest');
    
    % tor added may 1 - removed voxels was not legal otherwise
    %mask.removed_voxels = mask.removed_voxels(mask.volInfo.wh_inmask);
    % resample_space is not *always* returning legal sizes for removed
    % vox? maybe this was updated to be legal
    
    if length(mask.removed_voxels) == mask.volInfo.nvox
        disp('Warning: resample_space returned illegal length for removed voxels. Fixing...');
        mask.removed_voxels = mask.removed_voxels(mask.volInfo.wh_inmask);
    end
    pgmap = mask;
end
%% Binning
t = [];
voxel_percent = pgmap.prctile(0:(100./numBins):100);
bins_mask = []; 
%bins_gmmask = fmri_data(which('gray_matter_mask.nii'));
bins_gmmask = pgmap;
bins_gmmask.dat = 1:size(pgmap.dat,1);
nBinsVox = []; 
for i = 1:numBins
     temp_bins = [];      
     temp_bins =  ((voxel_percent(i) < pgmap.dat) &  ( pgmap.dat) <= voxel_percent(i+1));
     t(:,i)=double(temp_bins);
     bins_gmmask.dat = double(temp_bins);
     bins_mask{i} = bins_gmmask;       
end
%% Average -activation 
bins_results_nVox = []; 
bins_results_Percent = []; 
lengthmasked = []; 
for bin_i = 1:numBins
    temp_resutls = []; 
    %lengthmasked{bin_i} = length(find(thresh_img.dat ~=0));
    temp_results = apply_mask(thresh_img, bins_mask{bin_i});
    bins_results_nVox(:,bin_i) = sum(temp_results.dat ~=0);    
end
bins_results_Percent = bins_results_nVox./sum(bins_results_nVox);
%% results
results.bins_nVox = bins_results_nVox;
results.bins_Percentage = bins_results_Percent;

%% Ploting
if doplot 
    create_figure('map-wise percentil');
    set(gcf,'position',[185   726   449   165]);
    mSize = 50; 
    sepplot2(1:numBins, bins_results_Percent, 0.7,'markercolor',[245 23 123]./255,'linewidth',2,'markersize',mSize);
    box off;
    set(gca, 'linewidth', 2,'xlim',[0.7 numBins+0.3], 'XTick','','XTickLabel', '', 'YTick', [0:0.1:1],'YTickLabel',[0:0.1:1], ...
    'tickdir', 'out', 'ticklength', [.02 0.02], 'fontsize', 18)%,'ylim',[-0.15 1.5]);% ,'XTickLabelRotation',360-45);
end


end


%% SUBFUNCTION 
function [mkr,h] = sepplot2(x, y, prop, varargin)
% Draw points and a shorter line plots between points. This function almost
% similar with sepplot in CanlabCore. 
%
% :Usage:
% ::
%
%    [mkr,h] = sepplot(x, y, prop, varargin)
%
% :Inputs:
%
%   **x, y:**
%        The function plots vector Y against vector X
%
%   **prop:**
%        The proportion of the lines: prop can be between 0 and 1
%
% :Optional Inputs: Enter keyword followed by variable with values
%
%   **'linecolor':**
%        followed by color (e.g., 'color', [.5 .5 .5]) (default = black)
%
%   **'linewidth':**
%        followed by a number for linewidth (e.g., 'linewidth', 2) (default = .5)
%
%   **'linestyle':**
%        linestyle, e.g., followed by '-', '--', ':' (default = '-')
%
% :Output:
%
%   **mkr:**
%       graphic handles for markers
%  
%   **h:**
%        graphic handles for lines
%
% ::
%
%    x = 1:5; % x values
%    y = [32 40 55 84 130]; % mean
%    e = [6 6 6 6 6]; % standard error of the mean
%
%    create_figure(y_axis);
%    set(gcf, 'Position', [1   512   268   194]);
%    col = [0.3333    0.6588    1.0000];
%    markercol = col-.2;
%
%    h = errorbar(x, y, e, 'o', 'color', 'k', 'linewidth', 1.5, 'markersize', 7, 'markerfacecolor', col);
%    hold on;
%    sepplot2(x, y, .75, 'color', col, 'linewidth', 2);
%    errorbar_width(h, x, [0 0]);
%
%    set(gca, 'xlim', [.5 5.5], 'linewidth', 1.5);
%
%    try
%        pagesetup(gcf);
%        saveas(gcf, 'example.pdf');
%    catch
%        pagesetup(gcf);
%        saveas(gcf, 'example.pdf');
%    end
%
% ..
%    Copyright (C) 2014  Wani Woo (modified by Suhwan Gim, 2020)
% ..
dcol = [1 1 1]; % Marker color default (white)
lcol = [0 0 0]+0.1; % Line color default (black)
linew = .5; % linewidth default
lines = '-'; % linestyle default
edgeCol = 'none';
mSize = 36;
for i = 1:numel(varargin)
    if ischar(varargin{i})
        switch lower(varargin{i})
            case 'markercolor'
                dcol = varargin{i+1};
            case 'linecolor'
                lcol = varargin{i+1};
            case 'linewidth'
                linew = varargin{i+1};
            case 'linestyle'
                lines = varargin{i+1};
            case 'edgecolor'
                edgeCol = varargin{i+1};
            case 'markersize'
                mSize = varargin{i+1};
            otherwise
                error('Unknown inputs');
        end
    end
end
hold on;
%mkr=plot(x,y,'Marker','o','MarkerFaceColor',dcol,'MarkerEdgeColor', [edgeCol 0.8],'LineStyle','none');%,'MarkerSize',mSize);
mkr = scatter(x,y,mSize,dcol ,'o','filled','MarkerEdgeColor',[edgeCol 0.8]);
for i = 1:(numel(x)-1)
    xstep = (x(i+1)-x(i)).*((1-prop)/2);
    newx = [x(i)+xstep x(i+1)-xstep];
    
    ystep = (y(i+1)-y(i)).*((1-prop)/2);
    newy = [y(i)+ystep y(i+1)-ystep];
    
    h(i) = plot(newx, newy, 'color', lcol, 'linewidth', linew, 'linestyle', lines);
end

end
