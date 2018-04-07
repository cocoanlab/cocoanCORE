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

models.fns{1} = 'mediation(med_vars.model_ideology, med_vars.accuracy_mean, M, ''boot'', ''bootsamples'', 10000);';
models.fns{2} = 'mediation(med_vars.model_similarity, med_vars.accuracy_mean, M, ''boot'', ''bootsamples'', 10000);';
models.fns{3} = 'mediation(med_vars.model_ideology, med_vars.accuracy_sum, M, ''boot'', ''bootsamples'', 10000);';
models.fns{4} = 'mediation(med_vars.model_similarity, med_vars.accuracy_sum, M, ''boot'', ''bootsamples'', 10000);';

models.fns{5} = 'mediation(med_vars.model_ideology, med_vars.accuracy_mean, M, ''covs'', med_vars.model_race, ''boot'', ''bootsamples'', 10000);';
models.fns{6} = 'mediation(med_vars.model_similarity, med_vars.accuracy_mean, M, ''covs'', med_vars.model_race,  ''boot'', ''bootsamples'', 10000);';
models.fns{7} = 'mediation(med_vars.model_ideology, med_vars.accuracy_sum, M, ''covs'', med_vars.model_race, ''boot'', ''bootsamples'', 10000);';
models.fns{8} = 'mediation(med_vars.model_similarity, med_vars.accuracy_sum, M, ''covs'', med_vars.model_race, ''boot'', ''bootsamples'', 10000);';

models.fns{9} = 'mediation(med_vars.model_race, med_vars.accuracy_mean, M, ''covs'', med_vars.model_ideology, ''boot'', ''bootsamples'', 10000);';
models.fns{10} = 'mediation(med_vars.model_race, med_vars.accuracy_mean, M, ''covs'', med_vars.model_similarity,  ''boot'', ''bootsamples'', 10000);';
models.fns{11} = 'mediation(med_vars.model_race, med_vars.accuracy_sum, M, ''covs'', med_vars.model_ideology, ''boot'', ''bootsamples'', 10000);';
models.fns{12} = 'mediation(med_vars.model_race, med_vars.accuracy_sum, M, ''covs'', med_vars.model_similarity, ''boot'', ''bootsamples'', 10000);';

for i = 1:12
    models.savepaths{i} = [1,2,3,4,5];
end

models.name{1} = 'model1_LizPol_ideology_brain_accmean';
models.name{2} = 'model2_LizPol_similarity_brain_accmean';
models.name{3} = 'model3_LizPol_ideology_brain_accsum';
models.name{4} = 'model4_LizPol_similarity_brain_accsum';

models.name{5} = 'model5_LizPol_ideology_brain_accmean_cov_modelrace';
models.name{6} = 'model6_LizPol_similarity_brain_accmean_cov_modelrace';
models.name{7} = 'model7_LizPol_ideology_brain_accsum_cov_modelrace';
models.name{8} = 'model8_LizPol_similarity_brain_accsum_cov_modelrace';

models.name{9} = 'model9_LizPol_modelrace_brain_accmean_cov_ideology';
models.name{10} = 'model10_LizPol_modelrace_brain_accmean_cov_similarity';
models.name{11} = 'model11_LizPol_modelrace_brain_accsum_cov_ideology';
models.name{12} = 'model12_LizPol_modelrace_brain_accsum_cov_similarity';

code_filename = fullfile(outputdir, 'LizPol_081814_runmediation_wani.m');
study_scriptdir = outputdir;

%% run mediation

mediation_dream_wani(med_vars, models, jobn, mask, code_filename, study_scriptdir);

%% combine results
modeldir = '/dreamio3/wagerlab/labdata/current/Liz_political/fsl_mediation_images/Analysis/brain_analyses/beta_analyses';
mediation_dream_combineresults_wani(modeldir);

%% look at results
outputdir = '/dreamio3/wagerlab/labdata/current/Liz_political/fsl_mediation_images/Analysis/brain_analyses/beta_analyses';

% effect_img = 'X-M1-M2-Y_effect.img'; % or 'rob_beta_0001.img'
% pval_img = 'X-M1-M2-Y_pvals.img'; % or 'rob_p_0001.img'
% mask_img = 'X-M1-M2-Y_effect.img'; % or 'rob_beta_0001.img'
path = 'X-M-Y';
effect_img = [path '_effect.img']; % or 'rob_beta_0001.img'
pval_img = [path '_pvals.img']; % or 'rob_p_0001.img'
mask = fullfile(outputdir, 'LizPol_mask.nii');
% mask = effect_img;

dat = fmri_data(effect_img, mask);
datp = fmri_data(pval_img, mask);
datp.dat(datp.dat == 0) = Inf;

dat_stat = statistic_image;
dat_stat.dat = dat.dat;
dat_stat.p = datp.dat;
dat_stat.volInfo = dat.volInfo;

%dat_stat = threshold(dat_stat, .001, 'unc', 'k',5);
dat_stat = threshold(dat_stat, .05, 'fdr', 'k',1);
final_dat = double(dat_stat.sig) .* sign(dat_stat.dat);

dat_stat = threshold(dat_stat, .005, 'unc', 'k',1);
final_dat = final_dat + double(dat_stat.sig) .* sign(dat_stat.dat);

dat_stat = threshold(dat_stat, .01, 'unc', 'k',1);
final_dat = final_dat + double(dat_stat.sig) .* sign(dat_stat.dat);

dat_stat.dat = final_dat;

r = region(dat_stat);

k = [];
for i = 1:numel(r)
    if any(r(i).val == 3) || any(r(i).val == -3)
        k(end+1) = i;
    end
    
end

r = r(k);

orthviews(r);
%dat2 = region2imagevec(r);
%dat1 = region2imagevec(r);
%dat1 = replace_empty(dat1);
%dat2 = replace_empty(dat2);
%dat1.dat(:,2) = dat2.dat;
%orthviews(dat1);

%%


