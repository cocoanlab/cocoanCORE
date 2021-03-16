function results = principal_gradients_plot(thresh_img, varargin)
% 
% =========================================================================
%
%                   STILL WORKING (Be careful).                             
%                          
% =========================================================================
% Draw line plot using principal_gradients_axis map. The principal gradient
% map was calulated using HCP datasets (n = about 1000, 2mm whole-brain
% voxel-level connectivity). 
%  * You can use the previous one (3mm whole-brain voxel-level
%  connectivity) to compare. See below optional inputs. 
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
%   **plot (component 1) :** 
%       x axis: transmodal area <      ---------    > unimodal area 
%       y axis: percentage 
% :Inputs:
%
%   **thresh_img:**
%        thresholded image files 
%           (e.g., thresh_img = 'thresh_p_05.nii';
%                  thresh_img = {'thresh_p_05.nii';'thresh_p_01.nii'})
%               or 
%        fmri_data object 
%           (e.g., thresh_img = fmri_data(which('thresh_p_05.nii');
%        
%
%
% :Optional Inputs: Enter keyword followed by variable with values
% 
%   **'n56':**
%        use previous gradient map (3mm, n=56) 
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
% :See also:
%   example_principal_gradient_plot.m
%
% :Examples: 
%
%   % ** Use first gradient map
%   thresh_img = {'thresh_p_05.nii', 'thresh_p_01.nii'};
%   results = principal_gradients_plot(thresh_img);
%   pagesetup(gcf);
%   savenames = 'temp_map.pdf'
%   saveas(gcf, savenames);
%
%   % **Use other component with different number of bins
%   numComp = 2; % 1:10
%   results = principal_gradients_plot(thresh_img,'other_comp',numComp,'numBins',10); 
%

% History (yyyymmdd):
% 20210316 change prinipal_gradient_map 
%          (n=56, 3mm connectivity, one dataset -> n=1000, 2mm connectivity, HCP)
%
% The principal 


%% Default option 
doplot = true;
pgmapdir = which('Volumetric_hcp_gradients_GSP_90_DE_wholebrain_2mm.nii');  % made using 2mm HCP whole-brain voxel-level connectivity (n=1000)
mapn56 = false;
comp_num = 1; % first gradient axis (transmoal - unimodal)
numBins  = 20;
nImgs = numel(thresh_img);
mCol = hsv(nImgs);
mCol = mCol(randperm(nImgs),:);
%mCol = [245 23 123]./255;
use_exist = false;
%% Parse argumetns 
for i = 1:length(varargin)
    if ischar(varargin{i})
        switch lower(varargin{i})
            % functional commands
            case {'n56'}
                mapn56 = true;
                pgmapdir = which('ten_prinicpal_gradients_volumn_DE.nii');  % made using 3mm whole-brain voxel-level connectivity (n=56)
            case {'markercolor'}
                mCol = varargin{i+1};
            case {'noplots'}
                doplot = false;
            case {'other_comp'}
                comp_num = varargin{i+1};
            case {'numbins'}
                numBins = varargin{i+1};
            case {'samefigure'}
                use_exist = true;
            otherwise
                error('Unknown arguments')                 
        end
    end
end
%% Load gradients maps 
temp_pgmap = fmri_data(pgmapdir,which('gray_matter_mask.nii')); 
% aligning  
pgmap = temp_pgmap.get_wh_image(comp_num); 
if ~mapn56    
    pgmap.dat = pgmap.dat.*-1;
end
zeroidx = (pgmap.dat == 0);
pgmap = remove_empty(pgmap, zeroidx);
pgmap.dat = pgmap.dat - min(pgmap.dat); % make data positive


%% Load data 
if isobject(thresh_img) % 
    nImgs = 1;
    pattern_data = apply_mask(thresh_img,pgmap); 
else         
    thresh_img = thresh_img(:);
    pattern_data = fmri_data(thresh_img, which('gray_matter_mask.nii'));    
end
%% compare voxel space and size 
isdiff = compare_space(pgmap, pattern_data);

if isdiff == 1 || isdiff == 2 % diff space, not just diff voxels
    % == 3 is ok, diff non-empty voxels
    
    % Both work, but resample_space does not require going back to original
    % images on disk.    
    mask = resample_space(pgmap, pattern_data,'nearest');
    
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
%pattern_data = apply_mask(pattern_data, pgmap);
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
%% Average
bins_results_nVox = [];
bins_results_Percent = [];
lengthmasked = [];
for img_i = 1:nImgs
    for bin_i = 1:numBins
        temp_results = [];
        %lengthmasked{bin_i} = length(find(thresh_img.dat ~=0));
        temp_results = apply_mask(pattern_data.get_wh_image(img_i), bins_mask{bin_i});
        bins_results_nVox{img_i}(:,bin_i) = sum(temp_results.dat ~=0);
    end
    bins_results_Percent{img_i} = bins_results_nVox{img_i}./sum(bins_results_nVox{img_i});
end

%% results
results.bins_nVox = bins_results_nVox; % number of voxels 
results.bins_Percentage = bins_results_Percent; % percentage
results.binned_map = bins_mask; % binned mask 
%% Ploting
if doplot 
    if ~use_exist
        create_figure('map-wise percentil');
    end
    set(gcf,'position',[185   726   449   165]);
    mSize = 50; 
    hold on;
    for img_i = 1:nImgs       
        sepplot2(1:numBins, bins_results_Percent{img_i}, 0.7,'markercolor',mCol(img_i,:),'linewidth',2,'markersize',mSize);
    end
    hold off;
    box off;
    set(gca, 'linewidth', 2,'xlim',[0.7 numBins+0.3], 'XTick','','XTickLabel', '', 'YTick', [0:0.1:1],'YTickLabel',[0:0.1:1], ...
    'tickdir', 'out', 'ticklength', [.02 0.02], 'fontsize', 18,'ylim', [0 max(cat(2,bins_results_Percent{:}))]);% ,'XTickLabelRotation',360-45);
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
