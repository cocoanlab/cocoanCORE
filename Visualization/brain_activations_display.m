function [out, o2] = brain_activations_display(r, varargin)

% This function diplay brain activations on a inflated brain and few 
% saggital, axial slices. Cocoan style activation visualization.
%
% :Usage:
% ::
%
%    [out, o2] = brain_activations_display(cl, varargin)
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
%        depth for surface map, (e.g., 'depth', 4)
%        default is 3 mm
%
%   **color:**
%        if you want to use one color for blobs, you can specify color
%        using this option.
%
%   **region_color (or region_colors):**
%        if you want to use one color for each region, you can specify
%        region colors using this option. The input should have the same 
%        number of rows with the region, i.e., # region x 3
%
%   **custom_color (or custom_colors):**
%        if you want to define colors for all voxels on your own, 
%        you can specify voxel colors using this option. The input should 
%        have the same number of rows with the voxels, i.e., # voxel x 3
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
%        you can use this option to draw only surface. As a default, it
%        will show only lateral view. 
%
%   **surface_all:**
%        You can use this option to draw lateral and also medial view of
%        surface maps. 
% 
%   **all:**
%        This will put 2 surface and 12 montage maps (4 sagittal and 8 axial) 
%        in one row. This will be good when you want to explore the 
%        activation maps. 
%
%   **all2:**
%        This will show 2 surface and 10 montage maps (4 sagittal and 6 axial) 
%        in one row. This will be good when you want to use the figure for 
%        publicaition. 
%
%   **all_xyz**
%        You can also specify x and z for sagittal and axial slices for "all" 
%        option using this options. 
%        E.g., 'all_xyz', [-5 2 -35 35 -30:12:60]
%              first four will be used as x''s and the eight numbers after 
%              that will be used as z. More than 12 numbers will be ignored.
%
%   **all2_xyz**
%        You can also specify x and z for sagittal and axial slices for "all2" 
%        option using this options. 
%        E.g., 'all_xyz', [-5 2 -35 35 -30:12:60]
%              first four will be used as x''s and the six numbers after 
%              that will be used as z. More than 10 numbers will be ignored.
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
%   **y:**
%        you can specify the coronal slice numbers using this option. 
%        e.g., 'y', [-10 10] (default: []);
%
%   **z:**
%        you can specify the axial slice numbers
%        e.g., 'z', [-25 -15 -6 13 22]    
%        (default: slice range between z = -20 and 25 with spacing 10 mm)
%
%   **squeeze_x1:**
%        you can specify squeeze percentage for x1 using this. default: 40
%        e.g., 'squeeze_x1', 0 or 'squeeze_x1', 30
%
%   **squeeze_x2:**
%        you can specify squeeze percentage for x2 using this. default: 50
%
%   **squeeze_y:**
%        you can specify squeeze percentage for y using this. default: 30
%
%   **squeeze_z:**
%        you can specify squeeze percentage for z using this. default: 20
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
%
%   **colorbar:**
%        Show colorbar. default: false
%
%   **colorbar_fontsize:**
%        Font size of colorbar. default: 14
%
%   **prioritize_last:**
%        For determining colors of each vertex, prioritize the colors of
%        the voxels that are drawn last. Without specifying this, colors
%        are determined based on the colors of nearest voxels.
%
%
%  Examples:
%
%  % % Yeo 10 network, surface only
%  % gray_mask = fmri_data(which('Yeo_10networks_4mm.nii'));
%  % load(which('Schaefer_Net_Labels_r265.mat'));
%  % brain_activations_display(region(gray_mask, 'unique_mask_values'), 'surface_only', 'region_color', Schaefer_Net_Labels.ten_network_col);
%
%  % % Yeo 10 network, surface and montage
%  % gray_mask = fmri_data(which('Yeo_10networks_4mm.nii'));
%  % load(which('Schaefer_Net_Labels_r265.mat'));
%  % brain_activations_display(region(gray_mask, 'unique_mask_values'), 'all2', 'region_color', Schaefer_Net_Labels.ten_network_col);
%
%  % % SIIPS1 mask, surface only
%  % SIIPS1_mask = fmri_data(which('nonnoc_v11_4_137subjmap_weighted_mean.nii'));
%  % brain_activations_display(region(SIIPS1_mask, 'contiguous_regions'), 'surface_only', 'colorbar');
%
%  % % SIIPS1 mask, surface and montage
%  % SIIPS1_mask = fmri_data(which('nonnoc_v11_4_137subjmap_weighted_mean.nii'));
%  % brain_activations_display(region(SIIPS1_mask, 'contiguous_regions'), 'all2', 'colorbar');
%

global surface_style color depth poscm negcm do_color do_all all_style do_custom_color do_region_color prioritize_last

surface_style = 'veryinflated';
do_color = false;
do_custom_color = false;
do_region_color = false;
depth = 3;
do_montage = true;
do_surface = true;
do_medial_surface = false;
% do_all = false;
do_all = true;
all_style = 'v1';
do_colorbar = false;
colorbar_fontsize = 14;
prioritize_last = false;

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            % functional commands
            case {'inflated'}
                surface_style = 'inflated';
            case {'very_inflated'}
                % this is a default
                % surface_style = 'veryinflated';
            case {'very_inflated_workbench'}
                surface_style = 'veryinflated_wb';
            case {'color'}
                do_color = true;
                color = varargin{i+1};
                prioritize_last = false;
            case {'depth'}
                depth = varargin{i+1};
            case {'surface_only'}
                do_montage = false;
                do_all = false; 
            case {'surface_all'}
                do_medial_surface = true;
            case {'montage_only'}
                do_surface = false;
                do_all = false; 
            case {'all'}
                do_all = true;
                all_style = 'v1';
                disp('***************************************************************************************************************');
                disp('You selected to ''all'' option. It will draw four sagittal slices and eight axial slices with two surface maps.');
                disp('If you want to specify x and z for sagittal and axial slices, you can use ''all_xyz'' option.');
                disp('E.g., ''all_xyz'', [-5 2 -35 35 -30:12:60], first four will be used as x''s and the eight numbers after that');
                disp('      will be used as z. More than 12 numbers will be ignored.');
                disp('***************************************************************************************************************');
            case {'all2'}
                do_all = true;
                all_style = 'v2';     
                disp('***************************************************************************************************************');
                disp('You selected to ''all2'' option. It will draw four sagittal slices and six axial slices with two surface maps.');
                disp('If you want to specify x and z for sagittal and axial slices, you can use ''all2_xyz'' option.');
                disp('E.g., ''all2_xyz'', [-5 2 -35 35 -30:12:60], first four will be used as x''s and the six numbers after that');
                disp('      will be used as z. More than 10 numbers will be ignored.');
                disp('***************************************************************************************************************');
            case {'custom_color', 'custom_colors'} % surface only
                do_custom_color = true;
                color = varargin{i+1}; % this input should have the same 
                                        % number of rows with the voxel
                                        % i.e., # voxel x 3
                prioritize_last = false;
            case {'region_color', 'region_colors'}
                do_region_color = true;
                color = varargin{i+1}; % this input should have the same 
                                        % number of rows with the region
                                        % i.e., # region x 3
                prioritize_last = false;
                                        
            case {'colorbar'}
                do_colorbar = true;
                
            case {'colorbar_fontsize'}
                colorbar_fontsize = varargin{i+1};
                
            case {'prioritize_last'}
                prioritize_last = true;

        end
    end
end

out = [];
o2 = [];
s = get(0,'ScreenSize');

%% RIGHT

poscm = colormap_tor([0.96 0.41 0], [1 1 0]);  % warm
negcm = colormap_tor([.23 1 1], [0.11 0.46 1]);  % cools

if do_custom_color 
    poscm = color;
    negcm = [];
elseif do_region_color
    if numel(r) ~= size(color,1)
        error('The number of region and the rows of colors are different');
    end
    numVox = cat(1,r(:).numVox);
    poscm = [repelem(color(:,1), numVox) repelem(color(:,2), numVox) repelem(color(:,3), numVox)];
    negcm = [];
end

if do_surface && ~do_all    
    
    set(gcf, 'position', [1 s(4)/1.5 s(3)/3 s(4)/2]); % figure size
    
    if do_medial_surface        
        axes_positions = {[0.02 0.5 .46 .5], [0.52 0.5 .46 .5], [0.02 0.1 .46 .5], [0.52 0.1 .46 .5]};
    else
        axes_positions = {[0.02 0 .46 1], [0.52 0 .46 1]};
    end
    
    axes('Position', axes_positions{1});
    out = draw_surface(r, out, 'left');
    surface_light(gca);
    view(-90, 0);
    
    axes('Position', axes_positions{2});
    out = draw_surface(r, out, 'right');
    if strcmp(surface_style,'veryinflated')
        surface_light(gca);
    else
        camlight(-90,-20); axis vis3d;
    end
    view(90, 0);
    
    if do_medial_surface
        
        axes_h = get(gcf, 'Children');
        axes_new_h(1) = copyobj(axes_h(2), gcf);
        axes_new_h(2) = copyobj(axes_h(1), gcf);
        
        axes(axes_new_h(1));
        set(axes_new_h(1), 'Position', axes_positions{3});
        set(axes_new_h(1).Children(3), 'BackFaceLighting', 'reverselit');
        view(90, 0);
        
        axes(axes_new_h(2));    
        set(axes_new_h(2), 'Position', axes_positions{4});
        set(axes_new_h(2).Children(3), 'BackFaceLighting', 'reverselit');
        view(-90, 0);
        
    end
    
elseif do_all
    
    set(gcf, 'position', [0 s(4)/3*2 s(3) s(4)/3]);
    
    switch all_style
        
        case 'v1'
            axes_positions = {[0.87 0 .1 1], [0.01 0 .1 1]};
        case 'v2'
            axes_positions = {[0.85 0 .13 1], [0.02 0 .13 1]};
    end
    
    axes('Position', axes_positions{1});
    out = draw_surface(r, out, 'right');
    surface_light(gca);
    view(90, 0);
    
    axes('Position', axes_positions{2});
    out = draw_surface(r, out, 'left');
    surface_light(gca);
    view(-90, 0);
    
end

%% Montage: canlab visualization

% disply overlay
if do_montage
    o2 = brain_montage(r, varargin);
end

%% Colorbar
if do_colorbar
    
    if do_surface && ~do_all 
        p = 0.05 * 3; % unit for determining position and size of colorbar
        q = 0.7; % scaling factor of height: figure size of do_all is 3/2 of do_surface.
    elseif do_all
        p = 0.05; % unit for determining position and size of colorbar
        q = 1; % scaling factor of height: figure size of do_all is 3/2 of do_surface.
    end
    
    if isfield(out.colorbar, 'pos') && isfield(out.colorbar, 'neg')
        sf = 1 / (1 + 2*p); % scaling factor of overall size; here two colorbars
        cb_ax_pos = {[(1+p*2/5)*sf, (1-q*sf)/2, (p*2/5)*sf, q*sf], ...
            [(1+p*7/5)*sf, (1-q*sf)/2, (p*2/5)*sf, q*sf]};
        cb_pos = {[(1+p*3/5)*sf, (1-q*0.6*sf)/2, (p*1/5)*sf, q*0.6*sf], ...
            [(1+p*8/5)*sf, (1-q*0.6*sf)/2, (p*1/5)*sf, q*0.6*sf]};
        cb_lim = {out.colorbar.pos([1 end],1), out.colorbar.neg([1 end],1)};
        cb_map = {out.colorbar.pos(:,2:4), out.colorbar.neg(:,2:4)};
    elseif isfield(out.colorbar, 'pos') && ~isfield(out.colorbar, 'neg')
        sf = 1 / (1 + p); % scaling factor of overall size; here one colorbar
        cb_ax_pos = {[(1+p*2/5)*sf, (1-q*sf)/2, (p*2/5)*sf, q*sf]};
        cb_pos = {[(1+p*3/5)*sf, (1-q*0.6*sf)/2, (p*1/5)*sf, q*0.6*sf]};
        cb_lim = {out.colorbar.pos([1 end],1)};
        cb_map = {out.colorbar.pos(:,2:4)};
    elseif ~isfield(out.colorbar, 'pos') && isfield(out.colorbar, 'neg')
        sf = 1 / (1 + p); % scaling factor of overall size; here one colorbar
        cb_ax_pos = {[(1+p*2/5)*sf, (1-q*sf)/2, (p*2/5)*sf, q*sf]};
        cb_pos = {[(1+p*3/5)*sf, (1-q*0.6*sf)/2, (p*1/5)*sf, q*0.6*sf]};
        cb_lim = {out.colorbar.neg([1 end],1)};
        cb_map = {out.colorbar.neg(:,2:4)};
    else
        disp('No information for colorbar.');
        return;
    end
    
    fig_pos = get(gcf, 'position');
    set(gcf, 'position', [fig_pos(1:2), fig_pos(3:4) * sf]);
    axes_h = findobj('type', 'axes');
    for i = 1:numel(axes_h)
        set(axes_h, 'units', 'points');
    end
    set(gcf, 'position', [fig_pos(1:3), fig_pos(4) * sf]);
    for i = 1:numel(axes_h)
        set(axes_h, 'units', 'normalized');
    end
    for i = 1:numel(cb_ax_pos)
        cb_ax = axes('Position', cb_ax_pos{i});
        axis off;
        colormap(cb_ax, cb_map{i});
        if cb_lim{i}(1) == cb_lim{i}(2)
            caxis([cb_lim{i}(1)-eps*100 cb_lim{i}(1)+eps*100]);
            colorbar('Position', cb_pos{i}, 'AxisLocation', 'in', 'Tickdirection', 'out', 'TickLength', 0.015, 'FontSize', colorbar_fontsize, 'Ticks', cb_lim{i}(1));
        else
            caxis(cb_lim{i});
            colorbar('Position', cb_pos{i}, 'AxisLocation', 'in', 'Tickdirection', 'out', 'TickLength', 0.015, 'FontSize', colorbar_fontsize);
        end
    end

end

end

function o2 = brain_montage(r, vars)

global color do_color do_region_color

% default

dooutline = false;
do_pruned = false;
reuse_o2 = false;
do_cmaprange = false;

% parsing varargin

for i = 1:length(vars)
    if ischar(vars{i})
        switch vars{i}
            % functional commands
            case {'outline'}
                dooutline = true;
            case {'pruned'}
                do_pruned = true;
            case {'cmaprange'}
                do_cmaprange = true;
                cmaprange = vars{i+1};
        end
    else
        if isa(vars{i}, 'fmridisplay')
            reuse_o2 = true;
            o2 = vars{i};
        end 
    end
end

if ~reuse_o2 
    o2 = draw_montage(vars);
end

%%
o2 = removeblobs(o2);

if ~do_color
    if ~do_region_color
        if ~do_pruned && ~do_cmaprange
            o2 = addblobs(o2, r, 'splitcolor', {[.23 1 1], [0.17 0.61 1], [0.99 0.46 0], [1 1 0]}); % A&B
        elseif do_pruned
            o2 = addblobs(o2, r, 'splitcolor', {[.23 1 1], [0.17 0.61 1], [0.99 0.46 0], [1 1 0]}, 'cmaprange', [-2.8 -1.2 1.2 2.8]);
        elseif do_cmaprange
            o2 = addblobs(o2, r, 'splitcolor', {[.23 1 1], [0.17 0.61 1], [0.99 0.46 0], [1 1 0]}, 'cmaprange', cmaprange);
        end
    elseif do_region_color
        for ii = 1:numel(r)
            o2 = addblobs(o2, r(ii), 'color', color(ii,:)); 
        end
    end
else
    o2 = addblobs(o2, r, 'color', color); % A&B
end

if dooutline, o2 = addblobs(o2, r, 'outline', 'linewidth', 2, 'outline_color', [0 0 0]); end

end

function surface_light(gca)

out.h = get(gca, 'children');
set(out.h(2), 'BackFaceLighting', 'lit')
camlight(-90,-20);
axis vis3d;

end

function out = draw_surface(r, out, hemisphere)

global surface_style color depth poscm negcm do_color prioritize_last

switch hemisphere
    
    case 'left'
        
        if ~do_color
            switch surface_style
                case 'inflated'
                    [out.h_surf_L, out.colorbar] = cluster_surf_cocoan(r, 'underlay', which('surf_workbench_inflated_32k_Left.mat'), 'depth', depth, 'colormaps', poscm, negcm, 'prioritize_last', prioritize_last);
                case 'veryinflated'
                    [out.h_surf_L, out.colorbar] = cluster_surf_cocoan(r, 'underlay', 'fsavg_left', 'depth', depth, 'colormaps', poscm, negcm, 'prioritize_last', prioritize_last);
                case 'veryinflated_wb'
                    [out.h_surf_L, out.colorbar] = cluster_surf_cocoan(r, 'underlay', which('surf_workbench_very_inflated_32k_Left.mat'), 'depth', depth, 'colormaps', poscm, negcm, 'prioritize_last', prioritize_last);
            end
        else
            switch surface_style
                case 'inflated'
                    [out.h_surf_L, out.colorbar] = cluster_surf_cocoan(r, 'underlay', which('surf_workbench_inflated_32k_Left.mat'), 'depth', depth, 'color', color, 'prioritize_last', prioritize_last);
                case 'veryinflated'
                    [out.h_surf_L, out.colorbar] = cluster_surf_cocoan(r, 'underlay', 'fsavg_left', 'depth', depth, 'color', color, 'prioritize_last', prioritize_last);
                case 'veryinflated_wb'
                    [out.h_surf_L, out.colorbar] = cluster_surf_cocoan(r, 'underlay', which('surf_workbench_very_inflated_32k_Left.mat'), 'depth', depth, 'color', color, 'prioritize_last', prioritize_last);
            end
        end
    case 'right'
        
        if ~do_color
            switch surface_style
                case 'inflated'
                    [out.h_surf_R, out.colorbar] = cluster_surf_cocoan(r, 'underlay', which('surf_workbench_inflated_32k_Right.mat'), 'depth', depth,  'colormaps', poscm, negcm, 'prioritize_last', prioritize_last);
                case 'veryinflated'
                    [out.h_surf_R, out.colorbar] = cluster_surf_cocoan(r, 'underlay', 'fsavg_right', 'depth', depth, 'colormaps', poscm, negcm, 'prioritize_last', prioritize_last);
                case 'veryinflated_wb'
                    [out.h_surf_R, out.colorbar] = cluster_surf_cocoan(r, 'underlay', which('surf_workbench_very_inflated_32k_Right.mat'), 'depth', depth, 'colormaps', poscm, negcm, 'prioritize_last', prioritize_last);
            end
        else
            switch surface_style
                case 'inflated'
                    [out.h_surf_R, out.colorbar] = cluster_surf_cocoan(r, 'underlay', which('surf_workbench_inflated_32k_Right.mat'), 'depth', depth, 'color', color, 'prioritize_last', prioritize_last);
                case 'veryinflated'
                    [out.h_surf_R, out.colorbar] = cluster_surf_cocoan(r, 'underlay', 'fsavg_right', 'depth', depth, 'color', color, 'prioritize_last', prioritize_last);
                case 'veryinflated_wb'
                    [out.h_surf_R, out.colorbar] = cluster_surf_cocoan(r, 'underlay', which('surf_workbench_very_inflated_32k_Right.mat'), 'depth', depth, 'color', color, 'prioritize_last', prioritize_last);
            end
        end
end
        
end

function o2 = draw_montage(vars)

global do_all all_style

% default
axial_slice_range = [-20 25];
spacing = 10;
do_slice_range = true;
x1 = [-5 5]';
x2 = [-37 37]';
y = [];
do_label = false;

squeeze_x1 = 40;
squeeze_x2 = 50;
squeeze_y = 30;
squeeze_z = 20;
fontsize = 15;

for i = 1:length(vars)
    if ischar(vars{i})
        switch vars{i}
            % functional commands
            case {'axial_slice_range'}
                axial_slice_range = vars{i+1}{1};
                if numel(vars{i+1}) == 2
                    spacing = vars{i+2}{2};
                end
            case {'x1'}
                x1 = vars{i+1};
                if size(x1,1) == 1, x1 = x1'; end
            case {'x2'}
                x2 = vars{i+1};
                if size(x2,1) == 1, x2 = x2'; end
            case {'y'}
                y = vars{i+1};
                if size(y,1) == 1, y = y'; end
            case {'z'}
                do_slice_range = false;
                z = vars{i+1};
                if size(z,1) == 1, z = z'; end
            case {'squeeze_x1'}
                squeeze_x1 = vars{i+1};
            case {'squeeze_x2'}
                squeeze_x2 = vars{i+1};
            case {'squeeze_y'}
                squeeze_y = vars{i+1};
            case {'squeeze_z'} 
                squeeze_z = vars{i+1};
            case {'label', 'labels'}
                do_label = true;
            case {'fontsize'}
                fontsize = vars{i+1};
        end
    end
end

o2 = fmridisplay('overlay',which('keuken_2014_enhanced_for_underlay.img'));

if ~do_all
    
    xyz1 = x1;
    xyz2 = x2;
    xyz3 = y;
    
    xyz1(:, 2:3) = 0;
    xyz2(:, 2:3) = 0;
    
    xyz3(:, 2:3) = 0;
    xyz3 = xyz3(:,[2 1 3]);
    
    o2 = montage(o2, 'saggital', 'wh_slice', xyz1, 'onerow', 'brighten', .5);
    
    if ~isempty(x2)
        o2 = montage(o2, 'saggital', 'wh_slice', xyz2, 'onerow', 'brighten', .5);
    end
    
    if ~isempty(y)
        o2 = montage(o2, 'coronal', 'wh_slice', xyz3, 'onerow', 'brighten', .5);
    end
    
    if do_slice_range
        o2 = montage(o2, 'axial', 'slice_range', axial_slice_range, 'onerow', 'spacing', spacing, 'brighten', .5);
    else
        xyz4 = z;
        xyz4(:, 2:3) = 0;
        xyz4 = xyz4(:,[2 3 1]);
        
        o2 = montage(o2, 'axial', 'wh_slice', xyz4, 'onerow', 'brighten', .5);
    end
    
    k = 1;
    squeeze_axes_percent(o2.montage{k}.axis_handles, squeeze_x1);
    
    if ~isempty(x2)
        k = k + 1;
        squeeze_axes_percent(o2.montage{k}.axis_handles, squeeze_x2);
    end
    
    if ~isempty(y)
        k = k + 1;
        squeeze_axes_percent(o2.montage{k}.axis_handles, squeeze_y);
    end
    
    k = k + 1;
    squeeze_axes_percent(o2.montage{k}.axis_handles, squeeze_z);
    
else
    % predefined styles
    switch all_style
        case 'v1'

            s_first = [0.12 0.26];
            s_interval = [0.05 0.05];
            
            a_initial = 0.4;
            a_interval = 0.055;
            
            axes_positions = {[s_first(1) 0 .1 1], [s_first(1)+s_interval(1) 0 .1 1], ... % sagittal 1
                [s_first(2) 0 .1 1], [s_first(2)+s_interval(2) 0 .1 1], ... % sagittal 2
                [a_initial 0 .075 1], [a_initial+a_interval 0 .075 1], [a_initial+a_interval*2 0 .075 1], [a_initial+a_interval*3 0 .075 1], ... % axial
                [a_initial+a_interval*4 0 .075 1], [a_initial+a_interval*5 0 .075 1],  [a_initial+a_interval*6 0 .075 1], [a_initial+a_interval*7 0 .075 1]};
            
            xyz = [-2 2 -37 37 ... % sagittal
                -20 -10 0 10 20 30 40 50]; % axial
            if any(strcmp(vars, 'all_xyz')), xyz = vars{find(strcmp(vars, 'all_xyz'))+1}; end
            
            texts{1} = [38 -50;10 -50;35, -60;15, -60;-30, -115;-15, -115;-5, -115;-10, -115;-10, -115;-10, -115;-10, -115;-10, -115];
            for i = 1:numel(xyz)
                if i == 1 || i == 3
                    texts{2}{i} = ['x = ' num2str(xyz(i))];
                elseif i == 5
                    texts{2}{i} = ['z = ' num2str(xyz(i))];
                else
                    texts{2}{i} = num2str(xyz(i));
                end
            end
            
            for i = 1:numel(axes_positions)
                axh = axes('Position', axes_positions{i});
                if i < 5
                    o2 = montage(o2, 'saggital', 'wh_slice', [xyz(i) 0 0], 'onerow', 'brighten', .5, 'existing_axes', axh);
                else
                    o2 = montage(o2, 'axial', 'wh_slice', [0 0 xyz(i)], 'onerow', 'brighten', .5, 'existing_axes', axh);
                end
                if do_label, text(texts{1}(i,1), texts{1}(i,2), texts{2}{i}, 'fontsize', fontsize); end
            end
            
        case 'v2'
            
            s_first = [0.15 0.32];
            s_interval = [0.07 0.05];
            
            a_initial = 0.46;
            a_interval = 0.06;
            
            axes_positions = {[s_first(1) 0 .12 1], [s_first(1)+s_interval(1) 0 .12 1], ... % sagittal 1
                [s_first(2) 0 .1 1], [s_first(2)+s_interval(2) 0 .1 1], ... % sagittal 2
                [a_initial 0 .08 1], [a_initial+a_interval 0 .08 1], [a_initial+a_interval*2 0 .08 1], ... % axial
                [a_initial+a_interval*3 0 .08 1], [a_initial+a_interval*4 0 .08 1], [a_initial+a_interval*5 0 .08 1]};
            
            
            xyz = [-2 2 -37 37 ... % sagittal
                -20 -10 0 10 20 30]; % axial
            
            if any(strcmp(vars, 'all2_xyz')), xyz = vars{find(strcmp(vars, 'all2_xyz'))+1}; end
            
            texts{1} = [35 -50;25, -50;35, -60;15, -60;-25, -115;-15, -115;-2, -115;-7, -115;-7, -115;-7, -115];
            for i = 1:numel(xyz)
                if i == 1 || i == 3
                    texts{2}{i} = ['x = ' num2str(xyz(i))];
                elseif i == 5
                    texts{2}{i} = ['z = ' num2str(xyz(i))];
                else
                    texts{2}{i} = num2str(xyz(i));
                end
            end
            
            
            for i = 1:numel(axes_positions)
                axh = axes('Position', axes_positions{i});
                if i < 5
                    o2 = montage(o2, 'saggital', 'wh_slice', [xyz(i) 0 0], 'onerow', 'brighten', .5, 'existing_axes', axh);
                else
                    o2 = montage(o2, 'axial', 'wh_slice', [0 0 xyz(i)], 'onerow', 'brighten', .5, 'existing_axes', axh);
                end
                if do_label, text(texts{1}(i,1), texts{1}(i,2), texts{2}{i}, 'fontsize', fontsize); end
            end
            
    end
end
end