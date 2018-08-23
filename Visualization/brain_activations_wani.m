function [out, o2] = brain_activations_wani(r, varargin)

% This function diplay brain activations on a inflated brain and few 
% saggital, axial slices. Cocoan style activation visualization.
%
% :Usage:
% ::
%
%    [out, o2] = brain_activations_wani(cl, varargin)
%
% :Inputs:
%
%   **r:**
%        region object/activation map
%
% :Optional Inputs:
%
%   **inflated:**
%        not recommended
%        use inflated brain. We use the 32k inflated brain surface from HCP
%        connectome workbench. (Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii and 
%        Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii)
%
%   **very_inflated (default):**
%        recommended
%        use freesurfer inflated brain with Thomas Yeo group's RF_ANTs mapping
%        from MNI to Freesurfer. (https://doi.org/10.1002/hbm.24213)
%
%   **very_inflated_workbench:**
%        use very inflated brain. We also use the 32k inflated brain surface 
%        from HCP connectome workbench. 
%        (Q1-Q6_R440.R.very_inflated.32k_fs_LR.surf.gii and 
%         Q1-Q6_R440.L.very_inflated.32k_fs_LR.surf.gii)
%   
%   **depth:**
%        depth for surface map, (e.g., 'depth', 3)
%        default is 2 mm
%
%   **color:**
%        if you want to use one color for blobs, you can specify color
%        using this option.
%
%   **axial_slice_range:**
%        followed by axial slice range in a cell
%            e.g., 'axial_slice_range', {[-10 30]}
%        You can also define spacing in the same cell.
%            e.g., 'axial_slice_range', {[-10 30], 6}
%        The default range is [-20 25] with the spacing of 10.
%
%   **outline:**
%        draw outline, default linewidth: 2
%
%   **surface_only:**
%        you can use this option to draw only surface
%
%   **montage_only:**
%        you can use this option to draw only montage
%
%   **x1:**
%        you can specify the sagittal slice numbers using this option
%        e.g., 'x1', [-5 4] (default: [-5 5])
%
%   **x2:**
%        you can specify the second set of sagittal slice numbers using
%        this option. If you don't want to draw this, just leave it blank.
%        e.g., 'x2', [-41 41] or 'x2', []  (default: [-37 37]);
%
%   **z:**
%        you can specify the axial slice numbers
%        e.g., 'z', [-25 -15 -6 13 22]    
%        (default: slice range between z = -20 and 25 with spacing 10 mm)
%
%   **squeeze_x1:**
%        you can specify squeeze percentage using this. default: 40
%        e.g., 'squeeze_x1', 0 or 'squeeze_x1', 30
%
%   **squeeze_x2:**
%        you can specify squeeze percentage using this. default: 50
%
%   **squeeze_z:**
%        you can specify squeeze percentage using this. default: 20
%
%   **pruned:**
%        if you have pruned version of map, you can use this option. 
%        currently only works with (e.g., -3, -2, -1, 1, 2, 3)
%
%   **cmaprange:**
%        you can use this option to specify cmaprange. (see help of
%        addblob.m to see more details about cmaprange)
%
%   **o2:**
%        if you want to reuse montages underlay that are already exist, you
%        can simply provide o2 (fmridisplay object) as an input. It
%        automatically check whether there is an input that is fmridisplay
%        object, and reuse those montage. 

use_inflated = false;
use_veryinflated = true;
use_veryinflated_wb = false;
do_color = false;
depth = 2;
do_montage = true;
do_surface = true;

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            % functional commands
            case {'inflated'}
                use_inflated = true;
                use_veryinflated = false;
            case {'very_inflated'}
                % this is a default
                % use_veryinflated = true;
            case {'very_inflated_workbench'}
                use_veryinflated_wb = true;
                use_veryinflated = false;
            case {'color'}
                do_color = true;
                color = varargin{i+1};
            case {'depth'}
                depth = varargin{i+1};
            case {'surface_only'}
                do_montage = false;
            case {'montage_only'}
                do_surface = false;
        end
    end
end

%% RIGHT

out = [];

if do_surface
    poscm = colormap_tor([0.96 0.41 0], [1 1 0]);  % warm
    negcm = colormap_tor([.23 1 1], [0.11 0.46 1]);  % cools
    
    if ~do_color
        if use_inflated
            out.h_surf_R = cluster_surf(r ,which('surf_workbench_inflated_32k_Right.mat'), depth, 'heatmap', 'colormaps', poscm, negcm);
            out.h_surf_L = cluster_surf(r ,which('surf_workbench_inflated_32k_Left.mat'), depth, 'heatmap', 'colormaps', poscm, negcm);
        elseif use_veryinflated
            out.h_surf_R = cluster_surf(r, 'fsavg_right', depth, 'heatmap', 'colormaps', poscm, negcm);
            out.h_surf_L = cluster_surf(r, 'fsavg_left', depth, 'heatmap', 'colormaps', poscm, negcm);
        elseif use_veryinflated_wb
            out.h_surf_R = cluster_surf(r ,which('surf_workbench_very_inflated_32k_Right.mat'), depth, 'heatmap', 'colormaps', poscm, negcm);
            out.h_surf_L = cluster_surf(r ,which('surf_workbench_very_inflated_32k_Left.mat'), depth, 'heatmap', 'colormaps', poscm, negcm);
        end
    else
        if use_inflated
            out.h_surf_R = cluster_surf(r ,which('surf_workbench_inflated_32k_Right.mat'), depth, {color});
            out.h_surf_L = cluster_surf(r ,which('surf_workbench_inflated_32k_Left.mat'), depth, {color});
        elseif use_veryinflated
            out.h_surf_R = cluster_surf(r, 'fsavg_right', depth, {color});
            out.h_surf_L = cluster_surf(r, 'fsavg_left', depth, {color});
        elseif use_veryinflated_wb
            out.h_surf_R = cluster_surf(r ,which('surf_workbench_very_inflated_32k_Right.mat'), depth, {color});
            out.h_surf_L = cluster_surf(r ,which('surf_workbench_very_inflated_32k_Left.mat'), depth, {color});
        end

    end
    
    out.h = get(gca, 'children');
    set(out.h(2), 'BackFaceLighting', 'lit')
    %
    camlight(-90,-20);
    
    axis vis3d;
    
    set(gcf, 'Position', [532    25   931   930]);
    
    view(90, 0);
end

%% Montage: canlab visualization

% disply overlay
if do_montage
    o2 = brain_montage(r, varargin);
end

end

function o2 = brain_montage(r, vars)

% default
axial_slice_range = [-20 25];
dooutline = false;
spacing = 10;
do_slice_range = true;
do_color = false;
x1 = [-5 5]';
x2 = [-37 37]';

do_pruned = false;
reuse_o2 = false;
do_cmaprange = false;

squeeze_x1 = 40;
squeeze_x2 = 50;
squeeze_z = 20;

% parsing varargin

for i = 1:length(vars)
    if ischar(vars{i})
        switch vars{i}
            % functional commands
            case {'axial_slice_range'}
                axial_slice_range = vars{i+1}{1};
                if numel(vars{i+1}) == 2
                    spacing = vars{i+2}{2};
                end
            case {'outline'}
                dooutline = true;
            case {'color'}
                do_color = true;
                color = vars{i+1};
            case {'x1'}
                x1 = vars{i+1};
                if size(x1,1) == 1, x1 = x1'; end
            case {'x2'}
                x2 = vars{i+1};
                if size(x2,1) == 1, x2 = x2'; end
            case {'z'}
                do_slice_range = false;
                z = vars{i+1};
                if size(z,1) == 1, z = z'; end
            case {'squeeze_x1'}
                squeeze_x1 = varargin{i+1};
            case {'squeeze_x2'}
                squeeze_x2 = varargin{i+1};
            case {'squeeze_z'}
                squeeze_z = varargin{i+1};
            case {'pruned'}
                do_pruned = true;
            case {'cmaprange'}
                do_cmaprange = true;
                cmaprange = varargin{i+1};
        end
    else
        if isa(vars{i}, 'fmridisplay')
            reuse_o2 = true;
            o2 = vars{i};
        end 
    end
end

if ~reuse_o2
    
    o2 = fmridisplay;
    xyz1 = x1;
    xyz2 = x2;
    
    xyz1(:, 2:3) = 0;
    xyz2(:, 2:3) = 0;
    
    o2 = montage(o2, 'saggital', 'wh_slice', xyz1, 'onerow', 'brighten', .5);
    
    if ~isempty(x2)
        o2 = montage(o2, 'saggital', 'wh_slice', xyz2, 'onerow', 'brighten', .5);
    end
    
    if do_slice_range
        o2 = montage(o2, 'axial', 'slice_range', axial_slice_range, 'onerow', 'spacing', spacing, 'brighten', .5);
    else
        xyz3 = z;
        xyz3(:, 2:3) = 0;
        xyz3 = xyz3(:,[2 3 1]);
        
        o2 = montage(o2, 'axial', 'wh_slice', xyz3, 'onerow', 'brighten', .5);
    end
    
    squeeze_axes_percent(o2.montage{1}.axis_handles, squeeze_x1); 
    squeeze_axes_percent(o2.montage{2}.axis_handles, squeeze_x2); 
    squeeze_axes_percent(o2.montage{3}.axis_handles, squeeze_z); 
    
end

%%
o2 = removeblobs(o2);

if ~do_color
    if ~do_pruned && ~do_cmaprange
        o2 = addblobs(o2, r, 'splitcolor', {[.23 1 1], [0.17 0.61 1], [0.99 0.46 0], [1 1 0]}); % A&B
    elseif do_pruned
        o2 = addblobs(o2, r, 'splitcolor', {[.23 1 1], [0.17 0.61 1], [0.99 0.46 0], [1 1 0]}, 'cmaprange', [-2.8 -1.2 1.2 2.8]);
    elseif do_cmaprange
        o2 = addblobs(o2, r, 'splitcolor', {[.23 1 1], [0.17 0.61 1], [0.99 0.46 0], [1 1 0]}, 'cmaprange', cmaprange);
    end
else
    o2 = addblobs(o2, r, 'color', color); % A&B
end

if dooutline, o2 = addblobs(o2, r, 'outline', 'linewidth', 2, 'outline_color', [0 0 0]); end

end
