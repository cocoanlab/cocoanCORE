%% setting 
basedir = '/dreamio3/wagerlab/labdata/current/Liz_political/fsl_mediation_images/images_for_political_mediation';
addpath(basedir);

outputdir = '/dreamio3/wagerlab/labdata/current/Liz_political/fsl_mediation_images/Analysis/brain_analyses/beta_analyses';


load(fullfile(basedir, 'political_mediation_variables.mat'));

% % make the mask (run once)
% 
% for i = 1:numel(med_vars.imitation_brain)
%     dat = fmri_data(med_vars.imitation_brain{i});
%     dat.dat = all(dat.dat~=0,2);
%     data(:,i) = dat.dat;
% end
% 
% dat.dat = all(data,2);
% dat.fullpath = fullfile(outputdir, 'LizPol_mask.nii');
% write(dat);
% 
% indx = find(sum(data) < 200000);
% 
% k1 = 0;
% k2 = 0;
% k3 = 0;
% for i = indx
%     for j = 1:size(med_vars.imitation_brain{i},1)
%         dat = fmri_data(med_vars.imitation_brain{i}(j,:));
%         if sum(dat.dat ~=0)<280000
%             k1 = k1+1;
%             err_img_names.v280{k1} = med_vars.imitation_brain{i}(j,:);
%         end
%         if sum(dat.dat ~=0)<200000
%             k2 = k2+1;
%             err_img_names.v200{k2} = med_vars.imitation_brain{i}(j,:);
%         end
%         
%         if sum(dat.dat ~=0)<100000
%             k3 = k3+1;
%             err_img_names.v100{k3} = med_vars.imitation_brain{i}(j,:);
%         end
%     end
%     
% end
% 
% % 
% dat.dat = all(data,2);
med_vars.imgs = med_vars.imitation_brain;

%% running!

mask = fullfile(outputdir, 'LizPol_mask.nii');
jobn = 60;

% models.fns{1} = 'mediation(med_vars.model_ideology, med_vars.accuracy_proportion_correct, M, ''boot'', ''bootsamples'', 10000);';
% models.fns{2} = 'mediation(med_vars.model_similarity, med_vars.accuracy_proportion_correct, M, ''boot'', ''bootsamples'', 10000);';
% 
% models.fns{3} = 'mediation(med_vars.model_ideology, med_vars.accuracy_proportion_correct, M, ''covs'', med_vars.model_race, ''boot'', ''bootsamples'', 10000);';
% models.fns{4} = 'mediation(med_vars.model_similarity, med_vars.accuracy_proportion_correct, M, ''covs'', med_vars.model_race,  ''boot'', ''bootsamples'', 10000);';
% 
% models.fns{5} = 'mediation(med_vars.model_race, med_vars.accuracy_proportion_correct, M, ''covs'', med_vars.model_ideology, ''boot'', ''bootsamples'', 10000);';
% models.fns{6} = 'mediation(med_vars.model_race, med_vars.accuracy_proportion_correct, M, ''covs'', med_vars.model_similarity,  ''boot'', ''bootsamples'', 10000);';
% 
models.fns{1} = 'mediation_threepaths(med_vars.model_ideology, med_vars.accuracy_proportion_correct, med_vars.model_similarity, M, ''boot'', ''bootsamples'', 10000);';
models.fns{2} = 'mediation_threepaths(med_vars.model_ideology, med_vars.accuracy_proportion_correct, med_vars.model_similarity, M, ''covs'', med_vars.model_race, ''boot'', ''bootsamples'', 10000);';

models.fns{3} = 'mediation_threepaths(med_vars.model_ideology, med_vars.accuracy_proportion_correct, M, med_vars.model_similarity, ''boot'', ''bootsamples'', 10000);';
models.fns{4} = 'mediation_threepaths(med_vars.model_ideology, med_vars.accuracy_proportion_correct, M, med_vars.model_similarity, ''covs'', med_vars.model_race, ''boot'', ''bootsamples'', 10000);';

models.savepaths{1} = [2,3,6];
models.savepaths{2} = [2,3,6];
models.savepaths{3} = [1,2,3,6];
models.savepaths{4} = [1,2,3,6];
% 
% models.name{1} = 'model13_LizPol_ideology_brain_accprop';
% models.name{2} = 'model14_LizPol_similarity_brain_accprop';
% models.name{3} = 'model15_LizPol_ideology_brain_accprop_cov_modelrace';
% models.name{4} = 'model16_LizPol_similarity_brain_accprop_cov_modelrace';
% models.name{5} = 'model17_LizPol_modelrace_brain_accprop_cov_ideology';
% models.name{6} = 'model18_LizPol_modelrace_brain_accprop_cov_similarity';

models.name{1} = 'model19_LizPol_ideology_simil_brain_accprop_again';
models.name{2} = 'model20_LizPol_ideology_simil_brain_accprop_cov_modelrace_again';

models.name{3} = 'model21_LizPol_ideology_brain_simil_accprop_again';
models.name{4} = 'model22_LizPol_ideology_brain_simil_accprop_cov_modelrace_again';

code_filename = fullfile(outputdir, 'LizPol_082314_v3.m');
study_scriptdir = outputdir;

%% run mediation

mediation_dream_wani(med_vars, models, jobn, mask, code_filename, study_scriptdir);

%% submit the job
% $for i in ls Liz*; do distmsub2011 $i; done

%% combine results
modeldir = '/dreamio3/wagerlab/labdata/current/Liz_political/fsl_mediation_images/Analysis/brain_analyses/beta_analyses';
setupfile = fullfile(modeldir, 'mediation_SETUP_LizPol_082314_v3.mat');
mediation_dream_combineresults_wani(modeldir, setupfile);

%% look at results
outputdir = '/dreamio3/wagerlab/labdata/current/Liz_political/fsl_mediation_images/Analysis/brain_analyses/beta_analyses';
modeldir = fullfile(outputdir, 'model19_LizPol_ideology_simil_brain_accprop_again');
cd(modeldir);

path = 'X-M1-M2-Y';
effect_img = [path '_effect.img']; % or 'rob_beta_0001.img'
pval_img = [path '_pvals.img']; % or 'rob_p_0001.img'
% mask = fullfile(outputdir, 'LizPol_ma/sk.nii');
mask = which('gray_matter_mask.img');
% mask = effect_img;

dat = fmri_data(effect_img, mask);
datp = fmri_data(pval_img, mask);
datp.dat(datp.dat == 0) = Inf;

dat_stat = statistic_image;
dat_stat.dat = dat.dat;
dat_stat.p = datp.dat;
dat_stat.volInfo = dat.volInfo;

% dat_stat = threshold(dat_stat, .001, 'unc', 'k',1);
dat_stat = threshold(dat_stat, .05, 'fdr', 'k',1);
final_dat = double(dat_stat.sig) .* sign(dat_stat.dat);

dat_stat = threshold(dat_stat, .001, 'unc', 'k',1);
final_dat = final_dat + double(dat_stat.sig) .* sign(dat_stat.dat);

dat_stat = threshold(dat_stat, .005, 'unc', 'k',1);
final_dat = final_dat + double(dat_stat.sig) .* sign(dat_stat.dat);

dat_stat.dat = final_dat;

r = region(dat_stat);

k = [];
for i = 1:numel(r)
    if any(r(i).val == 3) % || any(r(i).val == -3)
        k(end+1) = i;
    end
    
end

r = r(k);

dat1 = region2imagevec(r);


%%
dat_m1_m2 = fmri_data('M1-M2_effect.img', mask); dat_m1_m2 = region2imagevec(region(dat_m1_m2));
dat_m2_y = fmri_data('M2-Y_effect.img', mask); dat_m2_y = region2imagevec(region(dat_m2_y));

dat1 = replace_empty(dat1);
dat_m1_m2 = replace_empty(dat_m1_m2);
dat_m2_y = replace_empty(dat_m2_y);

idx = (dat1.dat ~=0);
path_prod = dat_m1_m2.dat(idx) .* dat_m2_y.dat(idx);
dat_m1_m2.dat = dat_m1_m2.dat .* idx;
dat_m2_y.dat = dat_m2_y.dat .* idx;

dat_comb = dat_m1_m2;
dat_comb.dat = [dat_m1_m2.dat dat_m2_y.dat];
orthviews(dat_comb);


dat1.fullpath = fullfile(modeldir, 'X-M1-M2-Y_thresholded_fdr05_001_005.nii');
% dat1.fullpath = fullfile(modeldir3, 'sim_X-M_thresholded_001_005_01.nii');
write(dat1);

% dat2 = region2imagevec(r);
% dat1 = region2imagevec(r);
% dat1 = replace_empty(dat1);
% dat2 = replace_empty(dat2);
% dat1.dat(:,2) = dat2.dat;
% orthviews(dat1);

%% ideo/simil related brain activity

outputdir = '/dreamio3/wagerlab/labdata/current/Liz_political/fsl_mediation_images/Analysis/brain_analyses/beta_analyses';
modeldir = fullfile(outputdir, 'model13_LizPol_ideology_brain_accprop'); % .001, .005. .01 
% modeldir = fullfile(outputdir, 'model14_LizPol_similarity_brain_accprop'); % .001, .005, .01

cd(modeldir);

path = 'X-M';
effect_img = [path '_effect.img']; % or 'rob_beta_0001.img'
pval_img = [path '_pvals.img']; % or 'rob_p_0001.img'
% mask = fullfile(outputdir, 'LizPol_ma/sk.nii');
mask = which('gray_matter_mask.img');
maskhdr = which('gray_matter_mask.hdr');

dat = fmri_data(effect_img, mask);
datp = fmri_data(pval_img, mask);
datp.dat(datp.dat == 0) = NaN;

% change to z
% datz = datp;
% datz.dat = datz.dat/2;
% datz.dat(dat.dat>0) = 1-datz.dat(dat.dat>0);
% datz.dat = icdf('normal', datz.dat, 0, 1);
% datz.fullpath = fullfile(modeldir, 'ide_X-M_zvals.nii');
% write(datz, 'mni');

% eval(['!smoothest -m ' datz.fullpath ' -z ' datz.fullpath ' -V']);
% FWHMx = 9.93299 mm, FWHMy = 10.0517 mm, FWHMz = 9.78092 mm
% using AlphaSim, k = 72
% relevant code: cl_ext = cl_ext_3dClustSim(.05, [.001 .0005, .0001], [], mask, [2 2 2], ClustSimDir, 'fwhm', [9.93 10.05 9.78]', 'twotail');

dat_stat = statistic_image;
dat_stat.dat = dat.dat;
dat_stat.p = datp.dat;
dat_stat.volInfo = dat.volInfo;

dat_stat = threshold(dat_stat, .001, 'unc', 'k', 72);
% orthviews(dat_stat);
dat_stat.dat = dat_stat.dat .* dat_stat.sig;

dat_stat.fullpath = fullfile(modeldir2, 'ide_X-M_thresholded_001_k72.nii');
write(dat_stat);

%% ideo X-M-Y

outputdir = '/dreamio3/wagerlab/labdata/current/Liz_political/fsl_mediation_images/Analysis/brain_analyses/beta_analyses';
modeldir = fullfile(outputdir, 'model13_LizPol_ideology_brain_accprop'); % .001, .005. .01 
% modeldir = fullfile(outputdir, 'model14_LizPol_similarity_brain_accprop'); % .001, .005, .01

cd(modeldir);

path = 'X-M-Y';
effect_img = [path '_effect.img']; % or 'rob_beta_0001.img'
pval_img = [path '_pvals.img']; % or 'rob_p_0001.img'
% mask = fullfile(outputdir, 'LizPol_ma/sk.nii');
mask = which('gray_matter_mask.img');
maskhdr = which('gray_matter_mask.hdr');

dat = fmri_data(effect_img, mask);
datp = fmri_data(pval_img, mask);
datp.dat(datp.dat == 0) = NaN;

% % change to z
% datz = datp;
% datz.dat = datz.dat/2;
% datz.dat(dat.dat>0) = 1-datz.dat(dat.dat>0);
% datz.dat = icdf('normal', datz.dat, 0, 1);
% datz.fullpath = fullfile(modeldir, 'ide_X-M-Y_zvals.nii');
% write(datz, 'mni');
% 
% eval(['!smoothest -m ' datz.fullpath ' -z ' datz.fullpath ' -V']);
% FWHMx = 7.72604 mm, FWHMy = nan mm, FWHMz = nan mm
% using AlphaSim, k = 45
% relevant code: cl_ext = cl_ext_3dClustSim(.05, [.001 .0005, .0001], [], mask, [2 2 2], ClustSimDir, 'fwhm', [7.73 7.73 7.73]', 'twotail');

dat_stat = statistic_image;
dat_stat.dat = dat.dat;
dat_stat.p = datp.dat;
dat_stat.volInfo = dat.volInfo;

dat_stat = threshold(dat_stat, .05, 'fdr', 'k',1);
final_dat = double(dat_stat.sig) .* sign(dat_stat.dat);

dat_stat = threshold(dat_stat, .001, 'unc', 'k',1);
final_dat = final_dat + double(dat_stat.sig) .* sign(dat_stat.dat);

dat_stat = threshold(dat_stat, .01, 'unc', 'k',1);
final_dat = final_dat + double(dat_stat.sig) .* sign(dat_stat.dat);

dat_stat.dat = final_dat;

r = region(dat_stat);

k = [];
for i = 1:numel(r)
    if any(r(i).val == 3)  %|| any(r(i).val == -3)
        k(end+1) = i;
    end
    
end

r = r(k);

dat1 = region2imagevec(r);
dat1.fullpath = fullfile(modeldir, 'ide_X-M_Y_thresholded_fdr05_001_01.nii');
write(dat1);

dat_xmy_p = fmri_data('X-M-Y_pvals.img', dat1.fullpath); 
dat_xmy_p = remove_empty(dat_xmy_p); 

dat_xmy_p.dat = -log10(dat_xmy_p.dat);
cl = region(dat_xmy_p);

savename = 'ide_xmy_table.txt';
cluster_table_wani(cl, 1, 0, [1 32], 'writefile', savename);

dat_x_m = fmri_data('X-M_effect.img', mask); dat_x_m = region2imagevec(region(dat_x_m));
dat_m_y = fmri_data('M-Y_effect.img', mask); dat_m_y = region2imagevec(region(dat_m_y));


dat1 = replace_empty(dat1);
dat_x_m = replace_empty(dat_x_m);

dat1.dat = dat1.dat .* sign(dat_x_m.dat);
dat1.fullpath = fullfile(modeldir, 'ide_X-M_Y_thresholded_fdr05_001_01_sign_xm.nii');
write(dat1);

dat1 = region2imagevec(r);
dat1 = replace_empty(dat1);
dat_m_y = replace_empty(dat_m_y);

dat1.dat = dat1.dat .* sign(dat_m_y.dat);
dat1.fullpath = fullfile(modeldir, 'ide_X-M_Y_thresholded_fdr05_001_01_sign_my.nii');
write(dat1);

idx = (dat1.dat ~=0);
path_prod = dat_x_m.dat(idx) .* dat_m_y.dat(idx);
dat_x_m.dat = dat_x_m.dat .* idx;
dat_m_y.dat = dat_m_y.dat .* idx;

dat_comb = dat_x_m;
dat_comb.dat = [dat_x_m.dat dat_m_y.dat];
orthviews(dat_comb);


%% similarity
outputdir = '/dreamio3/wagerlab/labdata/current/Liz_political/fsl_mediation_images/Analysis/brain_analyses/beta_analyses';
modeldir = fullfile(outputdir, 'model14_LizPol_similarity_brain_accprop'); % .001, .005. .01 

cd(modeldir);

path = 'X-M-Y';
effect_img = [path '_effect.img']; % or 'rob_beta_0001.img'
pval_img = [path '_pvals.img']; % or 'rob_p_0001.img'
% mask = fullfile(outputdir, 'LizPol_ma/sk.nii');
mask = which('gray_matter_mask.img');
maskhdr = which('gray_matter_mask.hdr');

dat = fmri_data(effect_img, mask);
datp = fmri_data(pval_img, mask);
datp.dat(datp.dat == 0) = NaN;

% % change to z
% datz = datp;
% datz.dat = datz.dat/2;
% datz.dat(dat.dat>0) = 1-datz.dat(dat.dat>0);
% datz.dat = icdf('normal', datz.dat, 0, 1);
% datz.fullpath = fullfile(modeldir, 'ide_X-M-Y_zvals.nii');
% write(datz, 'mni');
% 
% eval(['!smoothest -m ' datz.fullpath ' -z ' datz.fullpath ' -V']);
% FWHMx = 7.72604 mm, FWHMy = nan mm, FWHMz = nan mm
% using AlphaSim, k = 45
% relevant code: cl_ext = cl_ext_3dClustSim(.05, [.001 .0005, .0001], [], mask, [2 2 2], ClustSimDir, 'fwhm', [7.73 7.73 7.73]', 'twotail');

dat_stat = statistic_image;
dat_stat.dat = dat.dat;
dat_stat.p = datp.dat;
dat_stat.volInfo = dat.volInfo;

dat_stat = threshold(dat_stat, .05, 'fdr', 'k',1);
final_dat = double(dat_stat.sig) .* sign(dat_stat.dat);

dat_stat = threshold(dat_stat, .001, 'unc', 'k',1);
final_dat = final_dat + double(dat_stat.sig) .* sign(dat_stat.dat);

dat_stat = threshold(dat_stat, .01, 'unc', 'k',1);
final_dat = final_dat + double(dat_stat.sig) .* sign(dat_stat.dat);

dat_stat.dat = final_dat;

r = region(dat_stat);

k = [];
for i = 1:numel(r)
    if any(r(i).val == 3)  %|| any(r(i).val == -3)
        k(end+1) = i;
    end
    
end

r = r(k);

dat1 = region2imagevec(r);
dat1.fullpath = fullfile(modeldir, 'sim_X-M_Y_thresholded_fdr05_001_01.nii');
write(dat1);

dat_xmy_p = fmri_data('X-M-Y_pvals.img', dat1.fullpath); 
dat_xmy_p = remove_empty(dat_xmy_p); 

dat_xmy_p.dat = -log10(dat_xmy_p.dat);
cl = region(dat_xmy_p);

savename = 'sim_xmy_table.txt';
cluster_table_wani(cl, 1, 0, [1 32], 'writefile', savename);

dat_x_m = fmri_data('X-M_effect.img', mask); dat_x_m = region2imagevec(region(dat_x_m));
dat_m_y = fmri_data('M-Y_effect.img', mask); dat_m_y = region2imagevec(region(dat_m_y));

dat1 = replace_empty(dat1);
dat_x_m = replace_empty(dat_x_m);

dat1.dat = dat1.dat .* sign(dat_x_m.dat);
dat1.fullpath = fullfile(modeldir, 'sim_X-M_Y_thresholded_fdr05_001_01_sign_xm.nii');
write(dat1);

dat1 = region2imagevec(r);
dat1 = replace_empty(dat1);
dat_m_y = replace_empty(dat_m_y);

dat1.dat = dat1.dat .* sign(dat_m_y.dat);
dat1.fullpath = fullfile(modeldir, 'sim_X-M_Y_thresholded_fdr05_001_01_sign_my.nii');
write(dat1);


idx = (dat1.dat ~=0);
path_prod = dat_x_m.dat(idx) .* dat_m_y.dat(idx);
dat_x_m.dat = dat_x_m.dat .* idx;
dat_m_y.dat = dat_m_y.dat .* idx;

dat_comb = dat_x_m;
dat_comb.dat = [dat_x_m.dat dat_m_y.dat];
orthviews(dat_comb);


%% simil

outputdir = '/dreamio3/wagerlab/labdata/current/Liz_political/fsl_mediation_images/Analysis/brain_analyses/beta_analyses';
% modeldir = fullfile(outputdir, 'model13_LizPol_ideology_brain_accprop'); % .001, .005. .01 
modeldir = fullfile(outputdir, 'model14_LizPol_similarity_brain_accprop'); % .001, .005, .01

cd(modeldir);

path = 'X-M';
effect_img = [path '_effect.img']; % or 'rob_beta_0001.img'
pval_img = [path '_pvals.img']; % or 'rob_p_0001.img'
% mask = fullfile(outputdir, 'LizPol_ma/sk.nii');
mask = which('gray_matter_mask.img');
maskhdr = which('gray_matter_mask.hdr');

dat = fmri_data(effect_img, mask);
datp = fmri_data(pval_img, mask);
datp.dat(datp.dat == 0) = NaN;

% change to z
datz = datp;
datz.dat = datz.dat/2;
datz.dat(dat.dat>0) = 1-datz.dat(dat.dat>0);
datz.dat = icdf('normal', datz.dat, 0, 1);
datz.fullpath = fullfile(modeldir, 'sim_X-M_zvals.nii');
write(datz, 'mni');

eval(['!smoothest -m ' datz.fullpath ' -z ' datz.fullpath ' -V']);
% FWHMx = 10.0046 mm, FWHMy = nan mm, FWHMz = 9.84924 mm98;.dsvx
% using AlphaSim, k = 72
% relevant code: cl_ext = cl_ext_3dClustSim(.05, [.001 .0005, .0001], [], mask, [2 2 2], ClustSimDir, 'fwhm', [9.93 10.05 9.78]', 'twotail');
% .01, 236

dat_stat = statistic_image;
dat_stat.dat = dat.dat;
dat_stat.p = datp.dat;
dat_stat.volInfo = dat.volInfo;

dat_stat = threshold(dat_stat, .01, 'unc', 'k', 236);
dat_stat.dat = dat_stat.dat .* dat_stat.sig;

dat_stat.fullpath = fullfile(modeldir, 'sim_X-M_thresholded_01_k236.nii');
write(dat_stat);

%% ideo/simil related brain activity

outputdir = '/dreamio3/wagerlab/labdata/current/Liz_political/fsl_mediation_images/Analysis/brain_analyses/beta_analyses';
modeldir = fullfile(outputdir, 'model13_LizPol_ideology_brain_accprop'); % .001, .005. .01 
% modeldir = fullfile(outputdir, 'model14_LizPol_similarity_brain_accprop'); % .001, .005, .01

cd(modeldir);

path = 'X-M';
effect_img = [path '_effect.img']; % or 'rob_beta_0001.img'
pval_img = [path '_pvals.img']; % or 'rob_p_0001.img'
% mask = fullfile(outputdir, 'LizPol_ma/sk.nii');
mask = which('gray_matter_mask.img');
maskhdr = which('gray_matter_mask.hdr');

dat = fmri_data(effect_img, mask);
datp = fmri_data(pval_img, mask);
datp.dat(datp.dat == 0) = NaN;

% change to z
datz = datp;
datz.dat = datz.dat/2;
datz.dat(dat.dat>0) = 1-datz.dat(dat.dat>0);
datz.dat = icdf('normal', datz.dat, 0, 1);
datz.fullpath = fullfile(modeldir, 'X-M_zvals.nii');
write(datz, 'mni');

eval(['!smoothest -m ' datz.fullpath ' -z ' datz.fullpath ' -V']);
% FWHMx = 9.93299 mm, FWHMy = 10.0517 mm, FWHMz = 9.78092 mm
% using AlphaSim, k = 72

dat_stat = statistic_image;
dat_stat.dat = dat.dat;
dat_stat.p = datp.dat;
dat_stat.volInfo = dat.volInfo;

dat_stat = threshold(dat_stat, .001, 'unc', 'k', 72);
orthviews(dat_stat);

% dat_stat.fullpath = fullfile(modeldir2, 'ide_X-M_thresholded_001.nii');
% write(dat_stat);

%%
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


