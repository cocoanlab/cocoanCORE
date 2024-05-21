%% Spatil smooth 3 mm with fslmaths  
% : https://kathleenhupfeld.com/how-to-smooth-images-in-fsl-its-different-from-spm/


%% SET DIRECTORY
basedir = '/Volumes/cocoanlab01/data/7T_HCP_emotion'; 
subdir = fullfile(basedir,'/imaging/preprocessed/sub-011/ses-12/func');
check_dir = fullfile(basedir,'/imaging/sanity_check');
nii_fname = filenames(fullfile(check_dir, 'sub-011_*HEAT*nii'));

%% How to calculate sigma for spatial smoothig
% FWHM = sigma * 2.354
FWHM = 3;
%  e.g.,   3 = sigma * 2.354
sigma = FWHM/2.354 ;
%% RUN fslmaths
setenv('FSLOUTPUTTYPE', 'NIFTI');
for i=1:length(nii_fname)
    tic;
    inputdir = nii_fname{i};
    [a,b,c] = fileparts(inputdir);
    %outputdir = sprintf('%s%s_%dsmooth%s',a,b,FWHM,c);
    outputdir = fullfile(a,[b '_smooth' num2str(FWHM) c]);
    
	system(sprintf('fslmaths %s -s %d %s',inputdir, sigma, outputdir));    
    toc;
end
%% load smoothed files
nii_fname = filenames(fullfile(check_dir, 'sub-011_*HEAT*smooth3.nii'));
dat = fmri_data(nii_fname{1});