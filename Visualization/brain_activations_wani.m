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
%   **axial_slice_range:**
%        followed by axial slice range in a cell
%            e.g., 'axial_slice_range', {[-10 30]}
%        You can also define spacing in the same cell.
%            e.g., 'axial_slice_range', {[-10 30], 6}
%        The default range is [-20 25] with the spacing of 10.
%
%   **outline:**
%        draw outline, linewidth: 2


use_inflated = false;
use_veryinflated = true;
use_veryinflated_wb = false;
axial_slice_range = [-20 25];
dooutline = false;
spacing = 10;
do_color = false;

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
            case {'axial_slice_range'}
                axial_slice_range = varargin{i+1}{1};
                if numel(varargin{i+1}) == 2
                    spacing = varargin{i+2}{2};
                end
            case {'outline'}
                dooutline = true;
            case {'color'}
                do_color = true;
                color = varargin{i+1};
        end
    end
end

%% RIGHT

poscm = colormap_tor([0.96 0.41 0], [1 1 0]);  % warm
negcm = colormap_tor([.23 1 1], [0.11 0.46 1]);  % cools

if ~do_color
    if use_inflated
        out.h_surf_R = cluster_surf(r ,which('surf_workbench_inflated_32k_Right.mat'), 3, 'heatmap', 'colormaps', poscm, negcm);
        out.h_surf_L = cluster_surf(r ,which('surf_workbench_inflated_32k_Left.mat'),3, 'heatmap', 'colormaps', poscm, negcm);
    elseif use_veryinflated
        out.h_surf_R = cluster_surf(r, 'fsavg_right', 2, 'heatmap', 'colormaps', poscm, negcm);
        out.h_surf_R = cluster_surf(r, 'fsavg_left', 2, 'heatmap', 'colormaps', poscm, negcm);
    elseif use_veryinflated_wb
        out.h_surf_R = cluster_surf(r ,which('surf_workbench_very_inflated_32k_Right.mat'), 3, 'heatmap', 'colormaps', poscm, negcm);
        out.h_surf_L = cluster_surf(r ,which('surf_workbench_very_inflated_32k_Left.mat'), 3, 'heatmap', 'colormaps', poscm, negcm);
    end
else
    if use_inflated
        out.h_surf_R = cluster_surf(r ,which('surf_workbench_inflated_32k_Right.mat'), 3, {color});
        out.h_surf_L = cluster_surf(r ,which('surf_workbench_inflated_32k_Left.mat'), 3, {color});
    elseif use_veryinflated
        out.h_surf_R = cluster_surf(r, 'fsavg_right', 2, 'heatmap', {color});
        out.h_surf_R = cluster_surf(r, 'fsavg_left', 2, 'heatmap', {color});
    elseif use_veryinflated_wb
        out.h_surf_R = cluster_surf(r ,which('surf_workbench_very_inflated_32k_Right.mat'), 3, {color});
        out.h_surf_L = cluster_surf(r ,which('surf_workbench_very_inflated_32k_Left.mat'), 3, {color});
    end
end

out.h = get(gca, 'children');
set(out.h(2), 'BackFaceLighting', 'lit')
% 
camlight(-90,-20);

axis vis3d;

set(gcf, 'Position', [532    25   931   930]);

view(90, 0);

%% Montage: canlab visualization

% disply overlay
o2 = fmridisplay;

xyz1 = [-5;5];
xyz2 = [-37;37];

xyz1(:, 2:3) = 0;
xyz2(:, 2:3) = 0;

o2 = montage(o2, 'saggital', 'wh_slice', xyz1, 'onerow', 'brighten', .5);
o2 = montage(o2, 'saggital', 'wh_slice', xyz2, 'onerow', 'brighten', .5);
o2 = montage(o2, 'axial', 'slice_range', axial_slice_range, 'onerow', 'spacing', spacing, 'brighten', .5);

squeeze_axes_percent(o2.montage{1}.axis_handles, 40);
squeeze_axes_percent(o2.montage{2}.axis_handles, 50);
squeeze_axes_percent(o2.montage{3}.axis_handles, 20);

%%
o2 = removeblobs(o2);
if ~do_color
    o2 = addblobs(o2, r, 'splitcolor', {[.23 1 1], [0.17 0.61 1], [0.99 0.46 0], [1 1 0]}); % A&B
else
    o2 = addblobs(o2, r, 'color', color); % A&B
end
if dooutline
    o2 = addblobs(o2, r, 'outline', 'linewidth', 2, 'outline_color', [0 0 0]); % A&B
end

end
