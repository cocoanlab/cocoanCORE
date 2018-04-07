
basedir = '/Volumes/RAID1/labdata/current/BMRK3';
roidir = '/Volumes/RAID1/labdata/current/BMRK3/Imaging/ROI_masks';

xyz = [-30 46 26; %dlpfc
    -6 -2 60; %SMA
    -14 8 -8; %NAc
    -2 -24 -2; %PAG
    14 30 -22; %vmPFC
    -16 -18 2; %Thalamus
    2 -32 64; %SMC
    46 -32 20]; %S2

[cl, all_data] = sphere_roi_tool_2008(which('brainmask.nii'), 4, xyz, 'useexisting');

outputdir = roidir;
filename = {'reappraisal_dlpfc.img', 'reappraisal_SMA.img', 'aversive_neg_NAc.img', 'aversive_neg_PAG.img',...
    'aversive_neg_vmpfc.img', 'aversive_neg_thal.img', 'aversive_pos_SMC.img', 'aversive_pos_S2.img'};

for i = 1:length(cl), wani_make_mask(outputdir, filename{i}, cl(i)); end

cd(basedir);
load wani_33_variables
load EXPT_33

subjn = length(EXPT_33.subjects);

savename = fullfile(basedir, 'Imaging', 'roi_new_variables.mat');
load(savename);

for i = 1:length(cl)
    mask{i} = filenames(fullfile(roidir, filename{i}), 'absolute', 'char');
end

% extract_roi_averages can take one mask image that has multiple ROIs.
% Then, it extracted ROI values all together. 

for i = 1:subjn
    subj = EXPT_33.subjects{i};
    subjectdir = fullfile(basedir, 'Imaging', subj, 'Models/Model12_Trial_vw');
    cd(subjectdir);
    imgs = filenames('beta_00*.img','char','absolute');
    for j = 1:length(mask)
        cll = extract_image_data(imgs(1:97,:), mask{j});
        [~,mask_name,~] = fileparts(mask{j});
        eval(['roi_new.' mask_name '{i} = nanmean(cll'')'';']);
    end
end

load(fullfile(basedir,'high_vif_trials.mat'));

for k = 1:subjn
    for j = 1:length(mask)
        [~,mask_name,~] = fileparts(mask{j});
        eval(['roi_new.' mask_name '{k}(high_medium_vif_trials{k}) = [];']);
    end
end

cd(fullfile(basedir, 'Imaging'));
save(savename, 'roi_new');