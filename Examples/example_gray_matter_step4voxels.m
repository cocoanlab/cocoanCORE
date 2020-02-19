% Example of how to use the 'gray_matter_step4voxels' for searchlight 
% analysis and mapping the results on whole brain again.
%
% 'gray_matter_step4voxels.mat' can be used for simplified searchlight
% analysis because it has only 3297 searchlights and instead of using whole
% voxels as a center of each searchlihgt.
%
% Copyright Byeol Kim, 2020
%% load brain data you want to use
% modeldir = fullfile(fastprojectdir, '/analysis/imaging/first_level/model07_spm_sing_tri_remv_vw');
load(which('gray_matter_step4voxels.mat'))
mask = which('gray_matter_mask.img');

% Brain Data loading example
% dat = fmri_data(filenames(fullfile(modeldir, 'sub-fast001', 'view_*.nii')), mask);
% dat = remove_empty(dat);

% Brain Data loading example
dat = fmri_data;
dat = apply_mask(dat, mask);
dat = remove_empty(dat);

%% Work with each searchlight

if size(dat.dat,1) == 211339
    for vox_i = 1:numel(vox_to_run)      % 3297
        seachlight = dat.dat(fillsphere_idx_211339(:,vox_i),:);
        
        % do analysis you want
    end
end

%% mapping the searchlight results on the whole brain
% If you get the results(3297xn), you can map the results on the whole
% brain with 'fill_step4voxels' function.
if size(result,1) == 3297
    new_dat = fill_step4voxels(result, 'values');
end

