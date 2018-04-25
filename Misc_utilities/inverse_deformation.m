function out = inverse_deformation(target_img, iy_nii_file, outputdir)

% Inverse warp. E.g., MNI to the native space.
%
% :Usage:
% ::
%
%    out = inverse_deformation(target_img, iy_nii_file, outputdir)
%
% :Inputs:
%
%   **target_img:**
%        The image in the standard (e.g., MNI) space that you want to
%        convert into the native space.
%
%   **iy_nii_file:**
%        the iy* nifti file name that is produced after the SPM12 
%        segmentation step. e.g., iy_sub-fast010_T1w.nii
%
%   **outputdir:**
%        the directory where you want to save the resulting image file
%        
%
% :Output:
%
%   **out:**
%        .matlabbatch         matlabbatch
%        .inv_deformation     
%        .target_image
%        .result_image
%
% :Examples:
% ::
%
%   iy_nii_file =  '/Volumes/cocoanlab/accumbens_sync/data/FAST/imaging/preprocessed/sub-fast010/anat/iy_sub-fast010_T1w.nii';
%   target_img = which('BN_Atlas_274_plus_brainstem_r280_wani.nii');
%   outputdir = '/Volumes/cocoanlab/accumbens_sync/data/FAST/imaging/preprocessed/sub-fast010/anat';
%   out = inverse_deformation(target_img, iy_nii_file, outputdir)
%
% ..
% Copyright (C) 2018  Choong-Wan Woo, Cocoan Lab
%
% ..

% Programmers' notes:
%   

matlabbatch{1}.spm.spatial.normalise.write.subj.def{1} = iy_nii_file;
matlabbatch{1}.spm.spatial.normalise.write.subj.resample{1} = target_img;
matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = NaN(2,3);
matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = NaN(1,3);
matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 4;
matlabbatch{1}.spm.spatial.normalise.write.woptions.prefix = 'inv_w';

out.matlabbatch = matlabbatch;

spm('defaults','fmri');
spm_jobman('initcfg');
spm_jobman('run', {matlabbatch});

out.inv_deformation = iy_nii_file;
out.target_image = target_img;

source = filenames(fullfile(fileparts(target_img), 'inv_w*'), 'char');

movefile(source, outputdir);

out.result_image = filenames(fullfile(outputdir, 'inv_w*'), 'char');
    
end