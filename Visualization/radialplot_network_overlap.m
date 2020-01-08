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
%   **'long' or 'longlabel:**
%        show the long network labels. The default is using the short
%        labels. 
%
%   **'color' or 'colors:**
%        you can specify the colors (the first cell should be the color for
%        positive, and the second cell should be the color for negative)
%        e.g., 'colors', {[1 0 0], [0 0 1]}; % pos neg
%

dolabel = true;
omit = 'asdlkjfldkj';
do_long = false;
colors = {[255 195 2]./255, [0 198 255]./255}; % pos neg

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
        end
    end
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

% Overlaps with large-scale networks with a polar plot
pattern_data = fmri_data(thresh_img, which('gray_matter_mask.img'));

pattern_thresh = pattern_data.dat;

% calculate posterior probability of observing thresholded regions given each network

do_posneg = [any(pattern_thresh>0) any(pattern_thresh<0)];

overlap_posneg = zeros(numel(network_names), 2);
for i = find(do_posneg)
    if i == 1
        overlap_posneg(:,i) = canlab_pattern_similarity(dat, pattern_thresh>0, 'posterior_overlap', 'ignore_missing');
    else
        overlap_posneg(:,i) = canlab_pattern_similarity(dat, pattern_thresh<0, 'posterior_overlap', 'ignore_missing');
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

[h.line, h.fill, h.ang] = tor_polar_plot({overlap_posneg(:,do_posneg).*100}, colors(do_posneg), {network_label}, 'nonumbers');

end