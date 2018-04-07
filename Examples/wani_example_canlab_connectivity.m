% Basic setting
basedir = '/Volumes/engram/labdata/data/bmrk5/Imaging/tr460';
outputdir = '/Volumes/engram/labdata/projects/BMRK5/bmrk5_rest/rsfc/results_v2';

% Get images
imgs = filenames(fullfile(basedir, '*', 'Functional', 'Preprocessed', '*', 'sw*rest*nii'), 'char'); % 48 images so far 082114

% ROI masks
subjn = size(imgs,1);
roi_masks = which('Buckner17subcl_SPMAnat_Morel_PAG_combined.nii');

%% AN EXAMPLE TO DO PREPROCESSING + EXRACTING ROI VALUES

for j = 1:subjn
    
    subject_dir = fileparts(fileparts(fileparts(fileparts(imgs(j,:))))); % there might be a better way
    
    [~, subj] = fileparts(subject_dir);
    
    fprintf('\n\n******************************************************\n');
    fprintf('Working on subject: %s\n', subj);
    fprintf('******************************************************\n');
    
    mask = which('gray_matter_mask.img');
    
    % reading data
    dat = fmri_data(imgs(j,:), mask);
    % plot(dat); % checking
    % % Looks good.

    % get nuisance
    nuisance_file = fullfile(subject_dir, 'Functional/Preprocessed/Nuisance_covariates_R.mat');
    
    % setting for canlab_connectivity_preproc
    load(nuisance_file);
    run_n = 1; % resting is the first run
    dat.covariates = [R{1}{run_n} R{2}{run_n}];
    
    TR = 0.46;
    
    % main part
    [preprocessed_dat, roi_val] = canlab_connectivity_preproc(dat, 'vw', 'datdir', subject_dir, 'windsorize', 5, 'linear_trend', 'hpf', .008, TR, 'extract_roi', roi_masks);
    
    % save
    savename = fullfile(outputdir, ['dat_obj_' subj '_preproc_Bu17cl_' date '.mat']);
    save(savename, 'preprocessed_dat', 'roi_val');
    

end

%% AN EXAMPLE TO DO ONLY EXRACTING ROI VALUES (NO PREPROCESSING)

for j = 1:subjn
    
    subject_dir = fileparts(fileparts(fileparts(fileparts(imgs(j,:))))); % there might be a better way
    
    [~, subj] = fileparts(subject_dir);

    fprintf('\n\n******************************************************\n');
    fprintf('Working on subject: %s\n', subj);
    fprintf('******************************************************\n');
    
    mask = which('gray_matter_mask.img');
    
    % reading data
    dat = fmri_data(imgs(j,:), mask);
    
    % main part
    [dat, roi_val] = canlab_connectivity_preproc(dat, 'no_preproc', 'extract_roi', roi_masks);
    
    % save
    savename = fullfile(outputdir, ['dat_obj_' subj '_preproc_Bu17cl_' date '.mat']);
    save(savename, 'roi_val');
end


