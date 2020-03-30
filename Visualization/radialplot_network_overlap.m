function h = radialplot_network_overlap(thresh_img, varargin)

% Draw radialplot using network overlaps
%
% :Usage:
% ::
%
%    h = radialplot_network_overlap(thresh_img, varargin)
%
% :Inputs:
%
%   **thresh_img:**
%        thresholded image files (e.g., thresh_img = 'thresh_p_05.nii';)
%        You can visualize one image for activation and deactivation maps.
%        You can also visualize multiple images with distinct colors. In this case,
%           this function calculates the overlaps with non-zero voxels. 
%
%
% :Optional Inputs: Enter keyword followed by variable with values
%
%   **'nolabel' or 'noname':**
%        do not show the network labels
%
%   **'omit':**
%        omit a specific network, 
%        when you want to exclude visual: e.g., 'omit', 'visual'
%        when you want to exclude visual and brainstem: e.g., 'omit', {'visual', 'brainstem'}
%
%   **'long' or 'longlabel':**
%        show the long network labels. The default is using the short
%        labels. 
%
%   **'color' or 'colors':**
%        you can specify the colors (the first cell should be the color for
%        positive, and the second cell should be the color for negative)
%        e.g., 'colors', {[1 0 0], [0 0 1]}; % pos neg
%
%  **'normalize':**
%        make the sum of overlap probability into 1 (for each map)
%

dolabel = true;
omit = 'asdlkjfldkj';
do_long = false;
colors = {[255 195 2]./255, [0 198 255]./255}; % pos neg
do_posneg = true;
do_normalize= false;

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            % functional commands
            case {'noname', 'nolabel'}
                dolabel = false;
            case {'omit'}
                omit = varargin{i+1};
            case {'long', 'longlabel'}
                do_long = true;
            case {'color', 'colors'}
                colors = varargin{i+1};
            case {'normalize'}
                do_normalize= true;
        end
    end
end

if size(thresh_img,1)>1
    do_posneg = false;
end

img = which('rBucknerlab_7clusters_SPMAnat_Other_combined.img');
% this image is in CanlabCore/canlab_canonical_brains/Combined_multiatlas_ROI_masks/rBucknerlab_7clusters_SPMAnat_Other_combined.img

if isempty(img)
    warning('No file: rBucknerlab_7clusters_SPMAnat_Other_combined.img');
    error('Please check whether CanlabCore is in your path!');
end

mask = fmri_data(img, which('gray_matter_mask.img'));

network_names = {'Visual', 'Somatomotor', 'dAttention', 'vAttention', ...
    'Limbic', 'Frontoparietal', 'Default', 'Thalamus', 'Hipp/Amy', 'Brainstem'};

% Overlaps with large-scale networks with a polar plot
pattern_data = fmri_data(thresh_img, which('gray_matter_mask.img'));

isdiff = compare_space(pattern_data, mask);

if isdiff == 1 || isdiff == 2 % diff space, not just diff voxels
    % == 3 is ok, diff non-empty voxels
    
    % Both work, but resample_space does not require going back to original
    % images on disk.
    %mask = resample_to_image_space(mask, dat);
    mask = resample_space(mask, pattern_data);
    
    % tor added may 1 - removed voxels was not legal otherwise
    %mask.removed_voxels = mask.removed_voxels(mask.volInfo.wh_inmask);
    % resample_space is not *always* returning legal sizes for removed
    % vox? maybe this was updated to be legal
    
    if length(mask.removed_voxels) == mask.volInfo.nvox
        disp('Warning: resample_space returned illegal length for removed voxels. Fixing...');
        mask.removed_voxels = mask.removed_voxels(mask.volInfo.wh_inmask);
    end
    
end

pattern_thresh = pattern_data.dat;

dat = [mask.dat==1 | mask.dat==8 | mask.dat==15 ...
    mask.dat==2 | mask.dat==9 | mask.dat==16 ...
    mask.dat==3 | mask.dat==10 | mask.dat==17 ...
    mask.dat==4 | mask.dat==11 | mask.dat==18 ...
    mask.dat==5 | mask.dat==12 | mask.dat==19 ...
    mask.dat==6 | mask.dat==13 | mask.dat==20 ...
    mask.dat==7 | mask.dat==14 | mask.dat==21 ...
    mask.dat>=22 & mask.dat<=35 ...
    mask.dat>=36 & mask.dat<=47 ...
    mask.dat==49];

% calculate posterior probability of observing thresholded regions given each network

if do_posneg
    any_posneg = [any(pattern_thresh>0) any(pattern_thresh<0)];
    
    overlap_posneg = zeros(numel(network_names), 2);
    for i = find(any_posneg)
        if i == 1
            overlap_posneg(:,i) = canlab_pattern_similarity(dat, pattern_thresh>0, 'posterior_overlap', 'ignore_missing');
        else
            overlap_posneg(:,i) = canlab_pattern_similarity(dat, pattern_thresh<0, 'posterior_overlap', 'ignore_missing');
        end
    end
else
    for i = 1:numel(thresh_img)
        overlap_posneg(:,i) = canlab_pattern_similarity(dat, pattern_thresh(:,i)~=0, 'posterior_overlap', 'ignore_missing');
    end
end

% plotting (polar plot)
close all;
if do_long
    network_label = network_names;
else
    network_label = {'Vis', 'Som', 'dAtt', 'vAtt', 'Lim', 'FP', 'Def', 'Thal', 'Hipp/Amy', 'BS'};
end

overlap_posneg(contains(lower(network_names), lower(omit)),:) = []; 
network_label(contains(lower(network_names), lower(omit))) = []; 

if ~dolabel
    network_label = repmat({''}, 1, numel(network_label)); 
end

create_figure('radialplot_network_overlap');%  set(gcf, 'color', 'w');
if do_posneg
    [h.line, h.fill, h.ang] = tor_polar_plot({overlap_posneg(:,any_posneg).*100}, colors(any_posneg), {network_label}, 'nonumbers');
else
    if ismatrix(colors)
        colors_org = colors; clear colors;
        for i = 1:size(colors_org,1), colors{i} = colors_org(i,:); end
    end
    
    if do_normalize
        overlap_posneg = overlap_posneg./repmat(sum(overlap_posneg), size(overlap_posneg,1), 1);
    end
    [h.line, h.fill, h.ang] = tor_polar_plot({overlap_posneg.*100}, colors, {network_label}, 'nonumbers');
end

end