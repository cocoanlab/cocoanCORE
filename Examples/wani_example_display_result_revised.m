%% This script is a generic script to display mediation or robust results
% You can use a part of this script to display whatever in the way you
% want.

% This script contains
% 1. display Mediation and robust results for three levels of threshold (standard way)
% 2. Only positive (activations) or negative predictors (deactivations)
% 3. extract data only in the mask image
% 4. display only one level thresholded image
% 5. display conjunctions

% Modules
% A) 1-1. get data (Module A) 
% B) 1-2. display on overlay (Module B)
% C) 2-1-3. choose only positive or negative (Module C)
% D) 3-1. get data from mask (Module D)
% E) 4-1. get data at only one threshold level (Module E)
% F) 5-1. get conjunction data (Module F)
% G) save png file (Module G)
% H) write an image with data (Module H) e.g.) to display with caret



%% =======================================================================================
%  1. display Mediation and robust results for three levels of threshold (standard way)
% ========================================================================================

clear all;

% ---------------------------
%  1-1. get data (Module A)
% ---------------------------

%% 1-1-1. input modeldir or go to the model dir directly.
modeldir = '/Volumes/RAID1/labdata/current/BMRK3/Imaging/Imagination_mediation/mediation_Model_a7_imgup_vs_std_boot_vw';
cd(modeldir);

%% 1-1-2. thresholding and assign different values

% Mediation results:
% choose among 1) X-M_effect.img (X-M_pvals.img), 2) M-Y_effect, 3) X-M-Y_effect

% Robust results:
% choose among 1) rob_beta_0001.img (rob_p_0001.img), 2) rob_beta_0002.img (rob_p_0002.img)

% "final_data" should have 1 for voxels at p < 0.01, 2 for voxels at p <
% .005, and 3 for voxels at p < .001 (k >5). 

effect_img = 'X-M-Y_effect.img'; % or 'rob_beta_0001.img'
pval_img = 'X-M-Y_pvals.img'; % or 'rob_p_0001.img'
mask_img = 'mask.img'; % or 'rob_beta_0001.img'
prim_p = [.001 .0005 .0001];k = [37 29 17];

% ================= threshold p < .01 unc  ================================== 
dat = statistic_image('image_names', effect_img, 'type', 'beta');
[datp, volInfo] = iimg_threshold(pval_img, 'mask', mask_img);

datp(datp == 0 & isnan(dat.dat)) = NaN;
dat.p = datp;
dat.dat(dat.dat > 0) = 1;
dat.dat(dat.dat < 0) = -1; 

dat = threshold(dat, prim_p(1), 'unc', 'k', k(1));
final_dat = dat.dat .* dat.sig;


% ================= threshold p < .005 unc  ================================== 
dat = replace_empty(dat);
dat = threshold(dat, prim_p(2), 'unc', 'k', k(2));
final_dat((dat.dat .* dat.sig) == 1) = 2;
final_dat((dat.dat .* dat.sig) == -1) = -2;


% ================= threshold p < .001, k > 5 ================================== 
dat = replace_empty(dat);
dat = threshold(dat, prim_p(3), 'unc', 'k', k(3));
final_dat((dat.dat .* dat.sig) == 1) = 3;
final_dat((dat.dat .* dat.sig) == -1) = -3;


% ================= apply an overlay mask ================================== 
data = fmri_data(which('mask_colin27.img'), effect_img);


% ================= prepare the frame ====================================== 
dat = replace_empty(dat);
dat = threshold(dat, prim_p(1), 'unc', 'k', k(1));
dat.dat = final_dat .* data.dat; 
dat.sig = dat.sig .* data.dat;


% ================= option to save an image with data or as a mask ========= 
% ================= (e.g., to display with caret) ========================== 
%
% % mask image
% cl = region(dat);
% outputdir = '/Volumes/RAID1/labdata/current/BMRK3/Results_PPTs/temp';
% filename = 'model_a5_pathb_001.img';
% wani_make_mask(outputdir, filename, cl);
%
%
% % image with data (Module H)
%
% dat.fullpath = fullfile(outputdir, filename);
% write(dat);
%

%% second option when the blobs are really small

% need to do overlay first
% this doesn't work for caret...

%% new

prim_p = [];
k = [];

effect_img = 'rob_beta_0001.img'; % or 'rob_beta_0001.img'
pval_img = 'rob_p_0001.img'; % or 'rob_p_0001.img'
mask_img = 'rob_beta_0001.img'; % or 'rob_beta_0001.img'

dat = statistic_image('image_names', effect_img, 'type', 'beta');
[datp, volInfo] = iimg_threshold(pval_img, 'mask', mask_img);

datp(datp == 0 & isnan(dat.dat)) = NaN;
dat.p = datp;


dat.dat(dat.dat > 0) = 1;
dat.dat(dat.dat < 0) = -1;

dat = threshold(dat, prim_p(1), 'unc', 'k', k(1));
final_dat = dat.dat .* dat.sig;

if length(k) > 1
    for i = 2:length(k)
        % ================= threshold p < .005 unc  ==================================
        dat = replace_empty(dat);
        dat = threshold(dat, prim_p(i), 'unc', 'k', k(i));
        final_dat((dat.dat .* dat.sig) == 1) = i;
        final_dat((dat.dat .* dat.sig) == -1) = -i;
    end
end

% ================= apply an overlay mask ================================== 
data = fmri_data(which('mask_colin27.img'), effect_img);

% ================= prepare the frame ====================================== 
dat = replace_empty(dat);
dat = threshold(dat, prim_p(1), 'unc', 'k', k(1));
dat.dat = final_dat .* data.dat .* dat.sig; 


%% threshold and draw

% resdir = '/Volumes/RAID1/labdata/current/BMRK3/Imaging/Imagination_mediation/mediation_Model_a1_boot_imgupdown_vw';
resdir = '/Volumes/RAID1/labdata/current/BMRK3/Imaging/Imagination_mediation/mediation_Model_a10_boot_imgupdown_vw_wopexp';
cd(resdir);

effect_img = 'X-M-Y_effect.img'; % or 'rob_beta_0001.img'
pval_img = 'X-M-Y_pvals.img'; % or 'rob_p_0001.img'
mask_img = 'mask.img'; % or 'rob_beta_0001.img'
data = fmri_data(which('mask_colin27.img'), effect_img);

% read data
o2 = removeblobs(o2);

for i = 1:2 % pos and neg
    dat = statistic_image('image_names', effect_img, 'type', 'beta');
    [datp, volInfo] = iimg_threshold(pval_img, 'mask', mask_img);
    datp(datp == 0 & isnan(dat.dat)) = NaN;
    dat.p = datp;
    dat.dat(dat.dat > 0) = 1;
    dat.dat(dat.dat < 0) = -1;
    dat = threshold(dat, 0.01, 'unc', 'k', 1);
    
    % mask
    dat.dat = dat.dat .* data.dat;
    dat.sig = dat.sig .* data.dat;
    dat.p = dat.p .* data.dat;
    
    % ================= only positive or negative... ===========================
    if i == 1
        k = dat.dat < 0; dat.dat(k) = 0; dat.p(k) = 1;
    elseif i == 2
        k = dat.dat > 0; dat.dat(k) = 0; dat.p(k) = 1;
    end
    
    % ================= threshold p < .01 unc  ==================================
    dat = replace_empty(dat);
    dat = threshold(dat, 0.01, 'unc', 'k', 1);
    cl = region(dat);
    if i == 1, o2 = addblobs(o2, cl, 'color', [1 0.38 0.27]);
    elseif i == 2, o2 = addblobs(o2, cl, 'color', [0.05 0.45 1]); end
    final_dat = dat.sig;
    
    % ================= threshold p < .005 unc  ==================================
    dat = replace_empty(dat);
    dat = threshold(dat, 0.005, 'unc', 'k', 1);
    cl = region(dat);
    if i == 1, o2 = addblobs(o2, cl, 'color', [1 0.7 0.17]);
    elseif i == 2, o2 = addblobs(o2, cl, 'color', [0 0.75 1]); end
    final_dat = final_dat | dat.sig;
    
    % ================= threshold p < .001, k > 5 ==================================
    dat = replace_empty(dat);
    dat = threshold(dat, 0.001, 'unc', 'k', 3);
    cl = region(dat);
    if i == 1, o2 = addblobs(o2, cl, 'color', [1 1 0]);
    elseif i == 2, o2 = addblobs(o2, cl, 'color', [0 1 1]); end
    final_dat = final_dat | dat.sig;
    
    % ================= outline ======================
    dat = replace_empty(dat);
    dat = threshold(dat, 0.01, 'unc', 'k', 1); % whatever
    dat.dat = double(final_dat);
    dat.sig = final_dat;
    cl = region(dat);
    o2 = addblobs(o2, cl, 'color', [0 0 0], 'outline', 'linewidth', 1);
    
end


%% not separate pos and neg

% resdir = '/Volumes/RAID1/labdata/current/BMRK3/Imaging/Imagination_mediation/mediation_Model_a1_boot_imgupdown_vw';
resdir = '/Volumes/RAID1/labdata/current/BMRK3/Imaging/Imagination_mediation/mediation_Model_a10_boot_imgupdown_vw_wopexp';
cd(resdir);

effect_img = 'X-M-Y_effect.img'; % or 'rob_beta_0001.img'
pval_img = 'X-M-Y_pvals.img'; % or 'rob_p_0001.img'
mask_img = 'mask.img'; % or 'rob_beta_0001.img'
data = fmri_data(which('mask_colin27.img'), effect_img);

% read data
o2 = removeblobs(o2);

dat = statistic_image('image_names', effect_img, 'type', 'beta');
[datp, volInfo] = iimg_threshold(pval_img, 'mask', mask_img);
datp(datp == 0 & isnan(dat.dat)) = NaN;
dat.p = datp;
dat.dat(dat.dat > 0) = 1;
dat.dat(dat.dat < 0) = -1;
dat = threshold(dat, 0.01, 'unc', 'k', 1);

% mask
dat.dat(~logical(data.dat)) = 0;
dat.sig(~logical(data.dat)) = 0;
dat.p(~logical(data.dat)) = 1;


% ================= threshold p < .01 unc  ==================================
dat = replace_empty(dat);
dat = threshold(dat, 0.01, 'unc', 'k', 1);
cl = region(dat);
o2 = addblobs(o2, cl, 'color', [1 0.38 0.27]);
final_dat = dat.sig;

% ================= threshold p < .005 unc  ==================================
dat = replace_empty(dat);
dat = threshold(dat, 0.005, 'unc', 'k', 1);
cl = region(dat);
o2 = addblobs(o2, cl, 'color', [1 0.7 0.17]);
final_dat = final_dat | dat.sig;

% ================= threshold p < .001, k > 5 ==================================
dat = replace_empty(dat);
dat = threshold(dat, 0.001, 'unc', 'k', 3);
cl = region(dat);
o2 = addblobs(o2, cl, 'color', [1 1 0]);
final_dat = final_dat | dat.sig;

% ================= outline ======================
dat = replace_empty(dat);
dat = threshold(dat, 0.01, 'unc', 'k', 1); % whatever
dat.dat = double(final_dat);
dat.sig = final_dat;
cl = region(dat);
o2 = addblobs(o2, cl, 'color', [0 0 0], 'outline', 'linewidth', 1);


%% ---------------------------
%  1-2. display on overlay (Module B)
% ---------------------------

% disply overlay
o2 = fmridisplay;

xyz = [-45 -30 -20 -10 -6 -2 0 2 6 10 20 30 45]';
xyz(:, 2:3) = 0;

o2 = montage(o2, 'saggital', 'wh_slice', xyz, 'onerow');
o2 = montage(o2, 'axial', 'slice_range', [-40 50], 'onerow', 'spacing', 6);
% o2 = montage(o2, 'coronal', 'slice_range', [-70 50], 'onerow', 'spacing', 10);


% % =========== coarser display =========== 
% xyz = [-20 -8 -2 0 2 8 20]';
% xyz(:, 2:3) = 0;
% 
% o2 = montage(o2, 'saggital', 'wh_slice', xyz, 'onerow');
% o2 = montage(o2, 'axial', 'slice_range', [-40 50], 'onerow', 'spacing', 10);


% % =========== squeeze_axes!!! =========== 
% squeeze_axes(o2.montage{1}.axis_handles, 10);
% squeeze_axes(o2.montage{2}.axis_handles, 10);


% % =========== if you want to display a specific slice ===========
% for i = 1:length(mask) 
%     
%     r{i} = region(mask{i});
%     
%     xyz_x(i,:) = r{i}.mm_center(1);
%     xyz_z(i,:) = r{i}.mm_center(3);
% 
% end
% 
% xyz_x = unique(xyz_x); xyz_x(:,2:3) = 0;
% xyz_z = unique(xyz_z); xyz_z(:,2:3) = 0; xyz_z = xyz_z(:,[2 3 1]);
% 
% % xyz_x = xyz_x([8 10], :);
% 
% clear o2;
% o2 = fmridisplay;
% 
% o2 = montage(o2, 'saggital', 'wh_slice', xyz_x, 'onerow');
% o2 = montage(o2, 'axial', 'wh_slice', xyz_z, 'onerow');

%% ---------------------------
%  1-3. draw!! with cmaprange
% ---------------------------

cl = region(dat);
cl.Z = awgn(cl.Z, 60); cl.val = awgn(cl.val, 60); % add some noise 

o2 = removeblobs(o2);

o2 = addblobs(o2, cl, 'splitcolor', {[0 1 1] [0.05 0.45 1] [1 0.38 0.27] [1 1 0]},...
    'cmaprange', [-2.8 -1.2 1.2 2.8]);

% o2 = addblobs(o2, cl, 'color', [0 0 0], 'outline', 'linewidth', 1.2);

% o2 = addblobs(o2, cl, 'color', [0 0 0], 'outline', 'linewidth', 1);
% o2 = addblobs(o2, cl, 'color', rgb('Tomato') );
% o2 = addblobs(o2, cl, 'splitcolor', {[0 0 1] [.3 0 .8] [.8 .1 0] [1 1 0]}, 'cmaprange', [2.5 1.5 -1.5 -2.5]);
% o2 = addblobs(o2, cl, 'splitcolor', {[0 0 1] [0 0 1] [1 0 0] [1 0 0]}, 'cmaprange', [2.5 1.5 -1.5 -2.5]);

%% Display beta/t/p map all together using transparency
% color: beta, transparency: t, outline: p


% For robust regression results

resdir = '/Volumes/RAID1/labdata/current/BMRK3/Imaging/Contrasts_Model15b/robust/robust_img/robust0001';
cd(resdir);

effect_img = 'rob_beta_0001.img'; % or 'rob_beta_0001.img'
pval_img = 'rob_p_0001.img'; % or 'rob_p_0001.img'
mask_img = 'rob_beta_0001.img'; % or 'rob_beta_0001.img'
t_img = 'rob_tmap_0001.img';

% ================= color beta, trans tmap, outline p ======================
dat_trans = fmri_data(effect_img, which('mask_colin27.img'));
dat_t = fmri_data(t_img, which('mask_colin27.img'));
r = region(dat_trans);
r_tmap = region(dat_t);
prim_p = .000001;
df = [1 32];

o2 = removeblobs(o2);
% o2 = addblobs_wani(o2, r, 'splitcolor', {[0 0.8 1] [0 0 1] [1 0 0] [1 0.8 0]}, 'alphamap', r_tmap, 0.005, [1 32]);
o2 = addblobs_wani(o2, r, 'splitcolor', {[.3 0 1] [0 .8 1] [1 .2 0] [1 0.8 0]}, 'alphamap', r_tmap, prim_p, df);

disp('=======================================================');
fprintf('\nmax t = %.2f(pos), min_t = %.2f(pos)\n', spm_u(prim_p/2, df, 'T'),spm_u(.5/2, df, 'T'));
fprintf('\nmax beta = %.4f(pos), %.4f(neg), min beta = %.4f(pos), %.4f(neg)\n\n', prctile(dat_trans.dat(dat_trans.dat > 0), 95),...
    prctile(dat_trans.dat(dat_trans.dat < 0), 5), prctile(dat_trans.dat(dat_trans.dat > 0), 20), prctile(dat_trans.dat(dat_trans.dat < 0), 80));
disp('=======================================================');

data = fmri_data(which('mask_colin27.img'), effect_img);


% ================= outline cluster p < .05  ================================== 
dat = statistic_image('image_names', effect_img, 'type', 'beta');
[datp, volInfo] = iimg_threshold(pval_img, 'mask', mask_img);

datp(datp == 0 & isnan(dat.dat)) = NaN;
dat.p = datp;

dat = threshold(dat, 0.005, 'unc', 'k', 75);
final_dat = dat.dat .* dat.sig;
dat.dat = final_dat .* data.dat; 
dat.sig = dat.sig .* data.dat;

cl = region(dat);
o2 = addblobs(o2, cl, 'color', [0 0 0], 'outline', 'linewidth', 1.2);


%% Display beta/t/p map all together using transparency
% color: beta, transparency: t, outline: p

% For mediation results
modeldir = '/Volumes/RAID1/labdata/current/BMRK3/Imaging/Imagination_mediation/mediation_Model_d3_imgdown_boot_temp_rating';

cd(modeldir);

effect_img = 'M-Y_effect.img'; 
pval_img = 'M-Y_pvals.img'; 
mask_img = 'mask.img'; 
ste_img = 'M-Y_ste.img';

% ================= color beta, trans tmap, outline p ======================
dat_trans = fmri_data(effect_img, which('mask_colin27.img'));
dat_ste = fmri_data(ste_img, which('mask_colin27.img'));
dat_t = dat_trans; dat_t.dat = dat_trans.dat./dat_ste.dat;

r = region(dat_trans);
r_tmap = region(dat_t);

o2 = removeblobs(o2);
% o2 = addblobs_wani(o2, r, 'splitcolor', {[0 0.8 1] [0 0 1] [1 0 0] [1 0.8 0]}, 'alphamap', r_tmap, 0.005, [1 32]);
o2 = addblobs_wani(o2, r, 'splitcolor', {[.3 0 1] [0 .8 1] [1 .2 0] [1 0.8 0]}, 'alphamap', r_tmap, 0.005, [1 32]);

data = fmri_data(which('mask_colin27.img'), effect_img);


% ================= outline cluster p < .05  ================================== 
dat = statistic_image('image_names', effect_img, 'type', 'beta');
[datp, volInfo] = iimg_threshold(pval_img, 'mask', mask_img);

datp(datp == 0 & isnan(dat.dat)) = NaN;
dat.p = datp;

dat = threshold(dat, 0.005, 'unc', 'k', 75);
final_dat = dat.dat .* dat.sig;
dat.dat = final_dat .* data.dat; 
dat.sig = dat.sig .* data.dat;

cl = region(dat);
o2 = addblobs(o2, cl, 'color', [0 0 0], 'outline', 'linewidth', 1.2);



%% =======================================================================================
%  2. Only positive (activations) or negative predictors (deactivations)
% ========================================================================================

% ---------------------------
%  2-1. get data
% ---------------------------

%% 2-1-1. input modeldir or go to the model dir directly.

modeldir = '/Volumes/RAID1/labdata/current/BMRK3/Imaging/Imagination_mediation/mediation_Model_a1_boot_imgupdown_vw';
cd(modeldir);

%% 2-1-2. thresholding and assign different values

% Mediation results:
% choose among 1) X-M_effect.img (X-M_pvals.img), 2) M-Y_effect, 3) X-M-Y_effect

% Robust results:
% choose among 1) rob_beta_0001.img (rob_p_0001.img), 2) rob_beta_0002.img (rob_p_0002.img)

dat = statistic_image('image_names', 'M-Y_effect.img', 'type', 'beta');
[datp, volInfo] = iimg_threshold('M-Y_pvals.img', 'mask', 'mask.img');

datp(datp == 0 & isnan(dat.dat)) = NaN;
dat.p = datp;

dat.dat(dat.dat > 0) = 1;
dat.dat(dat.dat < 0) = -1; 

% threshold p < .01 unc
% "final_data" should have 1 for voxels at p < 0.01, 2 for voxels at p <
% .005, and 3 for voxels at p < .001 (k >5). 
dat = threshold(dat, 0.01, 'unc', 'k', 1);
final_dat = dat.dat .* dat.sig;

% threshold p < .005 unc
dat = replace_empty(dat);
dat = threshold(dat, 0.005, 'unc', 'k', 1);
final_dat = final_dat + dat.dat .* dat.sig;

% threshold p < .001, k > 5
dat = replace_empty(dat);
dat = threshold(dat, 0.001, 'unc', 'k', 5);
final_dat = final_dat + dat.dat .* dat.sig;

% prepare
dat = replace_empty(dat);
dat = threshold(dat, 0.01, 'unc', 'k', 1);
dat.dat = final_dat; 


%% 2-1-3. choose only positive or negative (Module C)
% ----------------------------------------------------
dat.sig = dat.dat > 0 & dat.sig; 
dat.sig = dat.dat < 0 & dat.sig; 
% ----------------------------------------------------

% % Option to save a mask image 
% 
% cl = region(dat);
% outputdir = '/Volumes/RAID1/labdata/current/BMRK3/Results_PPTs/temp';
% filename = 'model_a5_pathb_001.img';
% 
% dat.fullpath = fullfile(outputdir, filename);
% write(dat);








%% =======================================================================================
%  3. extract data only in the mask image
% ========================================================================================

% ---------------------------
%  3-1. get data from mask (Module D)
% ---------------------------

%% 3-1-1. input modeldir or go to the model dir directly.

modeldir = '/Volumes/RAID1/labdata/current/BMRK3/Imaging/Imagination_mediation/mediation_Model_a1_boot_imgupdown_vw';
cd(modeldir);

%% 3-1-2. get mask image
mask_img = '/Volumes/RAID1/labdata/current/BMRK3/Imaging/ROI_masks/model_a5_pathb_deact_newvol.img';
% mask_img = 'mask.img';

%% 3-1-3. get data only in mask image

% Mediation results:
% choose among 1) X-M_effect.img (X-M_pvals.img), 2) M-Y_effect, 3) X-M-Y_effect

% Robust results:
% choose among 1) rob_beta_0001.img (rob_p_0001.img), 2) rob_beta_0002.img (rob_p_0002.img)

dat = statistic_image('image_names', mask_img, 'type', 'beta');
dateffect = fmri_data('X-M-Y_effect.img', mask_img);
dat.dat = dateffect.dat;

datp = fmri_data('X-M-Y_pvals.img', mask_img);

datp(datp == 0 & isnan(dat.dat)) = NaN;
dat.p = datp.dat;

% thresholding
dat = threshold(dat, 0.005, 'unc', 'k', 10);


%% 3-1-3. choose only positive or negative (Module C)
% ----------------------------------------------------
dat.sig = dat.dat > 0 & dat.sig; % positive
dat.sig = dat.dat < 0 & dat.sig; % negative
% ----------------------------------------------------

cl = region(dat);

% % make as a mask
% outputdir = '/Volumes/RAID1/labdata/current/BMRK3/Imaging/ROI_masks';
% filename = 'deact_imgupdown.img';
% wani_make_mask(outputdir, filename, cl);








%% =======================================================================================
%  4. display only one level thresholded image
% ========================================================================================

% ---------------------------
%  4-1. get data at only one threshold level (Module E)
% ---------------------------

%% 4-1-1. input modeldir or go to the model dir directly.
modeldir = '/Volumes/RAID1/labdata/current/BMRK3/Imaging/Imagination_mediation/mediation_Model_a7_imgup_vs_std_boot_vw';
cd(modeldir);

%% 4-1-2. get data

% Mediation results:
% choose among 1) X-M_effect.img (X-M_pvals.img), 2) M-Y_effect, 3) X-M-Y_effect

% Robust results:
% choose among 1) rob_beta_0001.img (rob_p_0001.img), 2) rob_beta_0002.img (rob_p_0002.img)

dat = statistic_image('image_names', 'rob_beta_0001.img', 'type', 'beta');
[datp, volInfo] = iimg_threshold('rob_p_0001.img', 'mask', 'rob_beta_0001.img');

datp(datp == 0 & isnan(dat.dat)) = NaN;
dat.p = datp;

% threshold at p < .005, k > 5
dat = threshold(dat, 0.005, 'unc', 'k', 5);

cl = region(dat);

o2 = removeblobs(o2);

o2 = addblobs(o2, cl, 'splitcolor', {[0 1 1] [0.05 0.45 1] [1 0.38 0.27] [1 1 0]});
o2 = addblobs(o2, cl, 'color', [0 0 0], 'outline', 'linewidth', 1);






%% =======================================================================================
%  5. display conjunctions
% ========================================================================================

% ---------------------------
%  5-1. get conjunction data (Module F)
% ---------------------------

%% 5-1-1. input modeldirs
clear;
modeldir{1} = '/Volumes/RAID1/labdata/current/BMRK3/Imaging/Contrasts_Model13/robust/robust0001';
modeldir{2} = '/Volumes/RAID1/labdata/current/BMRK3/Imaging/Contrasts_Model13/robust/robust0002';
modeldir{3} = '/Volumes/RAID1/labdata/current/BMRK3/Imaging/Contrasts_Model13/robust/robust0003';

%% 5-1-2. get data

for i = 1:length(modeldir)
    cd(modeldir{i});
    dat{i} = statistic_image('image_names', 'rob_beta_0001.img', 'type', 'beta');
    [datp, volInfo] = iimg_threshold('rob_p_0001.img', 'mask', 'rob_beta_0001.img');
    
    datp(datp == 0 & isnan(dat.dat)) = NaN;
    dat{i}.p = datp;

    % thresholding
    dat{i} = threshold(dat{i}, 0.01, 'unc', 'k', 5);
end


%% 5-1-3. get areas (conjunction and other)

% image frame
dat_frame = statistic_image('image_names', 'rob_beta_0001.img', 'type', 'beta');
[datp, volInfo] = iimg_threshold('rob_p_0001.img', 'mask', 'rob_beta_0001.img');

datp(datp == 0 & isnan(dat.dat)) = NaN;
dat_frame.p = datp;
dat_frame = threshold(dat_frame, 0.01, 'unc', 'k', 1); %whatever

n = length(modeldir);

if n == 2
    
    % conjunction (A & B)
    dat_frame.sig = dat{1}.sig & dat{2}.sig;
    dat_frame.dat = dat_frame.sig;
    cl{1} = region(dat_frame);
    
    % A & ~B
    dat_frame.sig = dat{1}.sig & (~dat{2}.sig);
    dat_frame.dat = dat_frame.sig;
    cl{2} = region(dat_frame);
    
    % ~A & B
    dat_frame.sig = (~dat{1}.sig) & dat{2}.sig;
    dat_frame.dat = dat_frame.sig;
    cl{3} = region(dat_frame);
    
    if ~exist(o2, 'var') % Module B
        % disply overlay
        o2 = fmridisplay;
        
        xyz = [-45 -30 -20 -10 -6 -2 0 2 6 10 20 30 45]';
        xyz(:, 2:3) = 0;
        
        o2 = montage(o2, 'saggital', 'wh_slice', xyz, 'onerow');
        o2 = montage(o2, 'axial', 'slice_range', [-40 50], 'onerow', 'spacing', 6);
        
        % % coarse display
        % xyz = [-20 -8 -2 0 2 8 20]';
        % xyz(:, 2:3) = 0;
        %
        % o2 = montage(o2, 'saggital', 'wh_slice', xyz, 'onerow');
        % o2 = montage(o2, 'axial', 'slice_range', [-40 50], 'onerow', 'spacing', 10);
    end
    
    o2 = removeblobs(o2);
    o2 = addblobs(o2, cl{1}, 'color', rgb('Gold')); % A&B
    o2 = addblobs(o2, cl{2}, 'color', rgb('Tomato')); % A&~B
    o2 = addblobs(o2, cl{3}, 'color', rgb('DodgerBlue')); % ~A&B

    
elseif n == 3
    
    % conjunction (A & B & C)
    dat_frame.sig = dat{1}.sig & dat{2}.sig & dat{3}.sig;
    dat_frame.dat = dat_frame.sig;
    cl{1} = region(dat_frame);
    
    % A & B & ~C
    dat_frame.sig = dat{1}.sig & dat{2}.sig & ~dat{3}.sig;
    dat_frame.dat = dat_frame.sig;
    cl{2} = region(dat_frame);
    
    % ~A & B & C
    dat_frame.sig = ~dat{1}.sig & dat{2}.sig & dat{3}.sig;
    dat_frame.dat = dat_frame.sig;
    cl{3} = region(dat_frame);
    
    % A & ~B & C
    dat_frame.sig = dat{1}.sig & ~dat{2}.sig & dat{3}.sig;
    dat_frame.dat = dat_frame.sig;
    cl{4} = region(dat_frame);
        
    % A & ~B & ~C
    dat_frame.sig = dat{1}.sig & ~dat{2}.sig & ~dat{3}.sig;
    dat_frame.dat = dat_frame.sig;
    cl{5} = region(dat_frame);
        
    % ~A & B & ~C
    dat_frame.sig = ~dat{1}.sig & dat{2}.sig & ~dat{3}.sig;
    dat_frame.dat = dat_frame.sig;
    cl{6} = region(dat_frame);
    
    % ~A & ~B & C
    dat_frame.sig = ~dat{1}.sig & ~dat{2}.sig & dat{3}.sig;
    dat_frame.dat = dat_frame.sig;
    cl{7} = region(dat_frame);
    
    if ~exist(o2, 'var') % Module B
        % disply overlay
        o2 = fmridisplay;
        
        xyz = [-45 -30 -20 -10 -6 -2 0 2 6 10 20 30 45]';
        xyz(:, 2:3) = 0;
        
        o2 = montage(o2, 'saggital', 'wh_slice', xyz, 'onerow');
        o2 = montage(o2, 'axial', 'slice_range', [-40 50], 'onerow', 'spacing', 6);
        
        % % coarse display
        % xyz = [-20 -8 -2 0 2 8 20]';
        % xyz(:, 2:3) = 0;
        %
        % o2 = montage(o2, 'saggital', 'wh_slice', xyz, 'onerow');
        % o2 = montage(o2, 'axial', 'slice_range', [-40 50], 'onerow', 'spacing', 10);
    end
        
    o2 = removeblobs(o2);
    o2 = addblobs(o2, cl{1}, 'color', rgb('Gold')); % A&B&C
    o2 = addblobs(o2, cl{2}, 'color', rgb('Orchid')); % A&B&~C
    o2 = addblobs(o2, cl{3}, 'color', rgb('MediumSpringGreen')); % ~A&B&C
    o2 = addblobs(o2, cl{4}, 'color', rgb('DarkGoldenrod')); % A&~B&C
    o2 = addblobs(o2, cl{5}, 'color', rgb('Tomato')); % A&~B&~C
    o2 = addblobs(o2, cl{6}, 'color', rgb('DodgerBlue')); % ~A&B&~C    
    o2 = addblobs(o2, cl{7}, 'color', rgb('ForestGreen')); % ~A&~B&C    
end


%% save png file (Module G)
scn_export_papersetup(600);
saveas(gcf, ['Seed_' comb_names{i} 'axial.png']);


%% write a image (Module H)

outputdir = '/Volumes/RAID1/labdata/current/BMRK3/Results_PPTs/Figures/caret';
dat.fullpath = fullfile(outputdir, 'path_ab_temp.img');
write(dat);




%% use transparency

% Example

mask_outline = which('weights_NSF_grouppred_cvpcr_FDR05.img');
mask = which('weights_NSF_grouppred_cvpcr.img');
dat = fmri_data(mask);
dat_outline = fmri_data(mask_outline);

max_trans_pos = min(dat_outline.dat(dat_outline.dat>0));
max_pos = max(dat.dat(dat.dat>0));
min_pos = min(dat.dat(dat.dat>0));

max_trans_neg = max(dat_outline.dat(dat_outline.dat<0));
max_neg = min(dat.dat(dat.dat<0));
min_neg = max(dat.dat(dat.dat<0));


%%
clear o2;
o2 = fmridisplay;

xyz = [-15 0 15]';
xyz(:, 2:3) = 0;

o2 = montage(o2, 'saggital', 'wh_slice', xyz, 'onerow');

o2 = montage(o2, 'axial', 'slice_range', [-20 50], 'onerow', 'spacing', 10);

squeeze_axes_wani(o2.montage{1}.axis_handles, -40);

%%

cl = region(dat);

o2 = removeblobs(o2);

o2 = addblobs(o2, cl, 'splitcolor', {[0 .75 1] [0.05 0.45 1] [1 0.38 0.27] [1 .65 0]},...
    'cmaprange', [max_trans_neg min_neg min_pos max_trans_pos], 'trans');

% o2 = addblobs(o2, cl, 'splitcolor', {[0 1 1] [0.05 0.45 1] [1 0.38 0.27] [1 1 0]},...
%     'cmaprange', [d 0 0 a], 'trans');

cl = region(dat_outline);

o2 = addblobs(o2, cl, 'splitcolor', {[0 1 1] [0 .75 1] [1 .65 0] [1 1 0]},...
     'cmaprange', [max_trans_neg min_neg min_pos max_trans_pos]);

o2 = addblobs(o2, cl, 'color', [0 0 0], 'outline', 'linewidth', 1.2);


%% Peak coordinate Table

outputdir = '/Volumes/RAID1/labdata/current/BMRK3/Results_PPTs/Tables';
resdir = '/Volumes/RAID1/labdata/current/BMRK3/Imaging/Contrasts_Model15b/robust/robust_img/robust0001';
cd(resdir);

% from wani_display_result.m
% 1-1-2. thresholding and assign different values

% Mediation results:
% choose among 1) X-M_effect.img (X-M_pvals.img), 2) M-Y_effect, 3) X-M-Y_effect

% Robust results:
% choose among 1) rob_beta_0001.img (rob_p_0001.img), 2) rob_beta_0002.img (rob_p_0002.img)

effect_img = 'rob_beta_0001.img'; % or 'rob_beta_0001.img'
pval_img = 'rob_p_0001.img'; % or 'rob_p_0001.img'
mask_img = 'rob_beta_0001.img'; % or 'rob_beta_0001.img'

data = fmri_data(which('mask_colin27.img'), effect_img);

dat = statistic_image('image_names', effect_img, 'type', 'beta');
[datp, volInfo] = iimg_threshold(pval_img, 'mask', mask_img);
datp(datp == 0 & isnan(dat.dat)) = NaN;
dat.p = datp;
dat = threshold(dat, 0.01, 'unc', 'k', 1);

% mask
dat.dat(~logical(data.dat)) = 0;
dat.sig(~logical(data.dat)) = 0;
dat.p(~logical(data.dat)) = 1;
date_now = scn_get_datetime;

% ================= threshold p < .01 unc  ==================================
dat = replace_empty(dat);
dat = threshold(dat, 0.005, 'unc', 'k', 75);
dat1 = dat;
dat1.dat = -log10(dat1.p) .* dat1.sig .* (-(dat1.dat < 0) + (dat1.dat > 0));
dat1.fullpath = fullfile(resdir, 'temp_for_table.img');
write(dat1);
cl = region(dat1.fullpath);

savename = fullfile(outputdir, ['table2_reapp_005_' date_now '.txt']);
cluster_table_wani(cl, 1, 0, 'writefile', savename);


% ================= threshold p < .005 unc  ==================================
dat = replace_empty(dat);
dat = threshold(dat, 0.001, 'unc', 'k', 37);
dat1 = dat;
dat1.dat = -log10(dat1.p) .* dat1.sig .* (-(dat1.dat < 0) + (dat1.dat > 0));
dat1.fullpath = fullfile(resdir, 'temp_for_table.img');
write(dat1);
cl = region(dat1.fullpath);

savename = fullfile(outputdir, ['table2_reapp_001_' date_now '.txt']);
cluster_table_wani(cl, 1, 0, 'writefile', savename);


% ================= threshold p < .001, k > 5 ==================================
dat = replace_empty(dat);
dat = threshold(dat, 0.0005, 'unc', 'k', 29);
dat1 = dat;
dat1.dat = -log10(dat1.p) .* dat1.sig .* (-(dat1.dat < 0) + (dat1.dat > 0));
dat1.fullpath = fullfile(resdir, 'temp_for_table.img');
write(dat1);
cl = region(dat1.fullpath);

savename = fullfile(outputdir, ['table2_reapp_0005_' date_now '.txt']);
cluster_table_wani(cl, 1, 0, 'writefile', savename);


%% pick cluster

cl = region(which('atlas_labels_combined.img'), 'unique_mask_values'); %anat_lbpa_thal.img'); %% 'lpba40.spm5.avg152T1.label.nii');

cluster_orthviews(cl, 'unique');

clout = [];

% choose MPFC clusters
[clout,cl] = cluster_graphic_select(cl,clout);

% display
cluster_orthviews(clout, 'unique');

for i =5:length(clout)
    k = (clout(i).XYZmm(2,:) < -20);
    clout(i).XYZ(:,k) = [];
    clout(i).XYZmm(:,k) = [];
    clout(i).val(k) = [];
    clout(i).Z(k) = [];
    clout(i).numVox = length(clout(i).Z);
end


% save
outputdir = '/Volumes/RAID1/labdata/current/Metaanalysis_Anjali/Anjali_MPFC_subcortical_connectivity/data/ROI_masks';
filename = 'MPFC_mask.img';
wani_make_mask(outputdir, filename, clout);

%% on dream/show insula and montage

o2 = fmridisplay;

xyz_x1 = [-5; 3]; xyz_x1(:, 2:3) = 0;
xyz_x2 = [-40; 44]; xyz_x2(:, 2:3) = 0;

o2 = montage(o2, 'saggital', 'wh_slice', xyz_x1, 'onerow');
o2 = montage(o2, 'saggital', 'wh_slice', xyz_x2, 'onerow');
o2 = montage(o2, 'axial', 'slice_range', [-25 40], 'onerow', 'spacing', 8);

squeeze_axes_wani(o2.montage{1}.axis_handles, 20);
squeeze_axes_wani(o2.montage{2}.axis_handles, 40);

for i = 1:3
    h = get(o2.montage{i}.axis_handles(1), 'parent');
    pos = get(h, 'position');
    pos(1:2) = [0 100];
    set(h, 'position', pos);
end