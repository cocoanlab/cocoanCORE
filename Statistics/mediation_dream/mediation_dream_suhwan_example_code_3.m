%% ====================================================================== %
%
%                 Whole brain mediatoin with only "SELF ratings" 
% 
% ======================================================================= %
%% setting
addpath(genpath('/sas1/cocoanlab/data/SEMIC_for_HPC'));

wh_loc = 'HPC';
[dir,str] =  SEMIC_HPC_set_dir(wh_loc);
for i=1:numel(str), eval(str{i}); end
basedir = fullfile(dir.prj_dir,'/scripts/4_custom_function/Mediation_dream_test');
clear dir

%% running!
jobn = 13;
sec_i = 32;
mask = which('gm_mask_semic2.nii');


loadnames = fullfile('/sas1/cocoanlab/data','SEMIC','data','meta_data','HPC_SEMIC_191203_mediation_dream_2mmWholeBrain_only_self.mat');
load(loadnames,'med_vars')

c=1;
models.name{1} = 'model01_X:_Stim_Y:_angle_M:_brain_Cov:_cue';
models.name{2} = 'model02_X:_Cue_Y:_angle_M:_brain_Cov:_Stim';

models.fns{1} = ['mediation(med_vars.X1, med_vars.Y, M,''covs'', med_vars.X2, ''boot'', ''bootsamples'', 10000)'];
models.fns{2} = ['mediation(med_vars.X2, med_vars.Y, M,''covs'', med_vars.X1, ''boot'', ''bootsamples'', 10000)'];



% it should change directory or filename


outputdir = fullfile('/sas1/cocoanlab/data','SEMIC','191203_2mm_whole_brain_onlySELF','outputs');
if ~exist(outputdir, 'dir'), mkdir(outputdir); end

med_vars.imgs = med_vars.M;
%models.name{1} = sprintf('Model01_STIM_phase_%02d_of_08_whole_brain',sec_i);
%models.name{2} = sprintf('Model02_CUE_phase_%02d_of_08_whole_brain',sec_i);
for i = 1:length(models.fns)
    models.savepaths{i} = [1,2,3,4,5]; % including all path
end


code_filename = fullfile(outputdir, sprintf('SEMIC_mediation_brain_%02d_of_32_run_suhwan.m',sec_i));
study_scriptdir = fullfile('/sas1/cocoanlab/data/SEMIC','scripts');

%% make mediation distributed scripts
mediation_dream_suhwan(med_vars, models, jobn, mask, code_filename, study_scriptdir,'wh_loc',wh_loc);

%% combine results
%modeldir = '/dreamio3/wagerlab/labdata/current/Liz_political/fsl_mediation_images/Analysis/brain_analyses/beta_analyses';
%%study_scriptdir 

%modeldir = fullfile('/sas1/cocoanlab/data/SEMIC/','191023_2mm_whole_brain_onlySELF','outputss');
modeldir = fullfile('/sas1/cocoanlab/data/SEMIC/','191203_2mm_whole_brain_onlySELF','outputs');
sec_i = 32;
SETUP_dir = fullfile(modeldir,sprintf('mediation_SETUP_SEMIC_mediation_brain_%02d_of_32_run_suhwan.mat',sec_i));

mediation_dream_combine_results_suhwan(modeldir, SETUP_dir);
%% look at results
outputdir = pwd;% '/dreamio3/wagerlab/labdata/current/Liz_political/fsl_mediation_images/Analysis/brain_analyses/beta_analyses';

% effect_img = 'X-M1-M2-Y_effect.img'; % or 'rob_beta_0001.img'
% pval_img = 'X-M1-M2-Y_pvals.img'; % or 'rob_p_0001.img'
% mask_img = 'X-M1-M2-Y_effect.img'; % or 'rob_beta_0001.img'
path = 'X-M-Y';
effect_img = [path '_effect.nii']; % or 'rob_beta_0001.img'
pval_img = [path '_pvals.nii']; % or 'rob_p_0001.img'
%mask = which('gm_mask_semic2.nii');%fullfile(outputdir, 'LizPol_mask.nii');
mask = which('gray_matter_mask.nii');%fullfile(outputdir, 'LizPol_mask.nii');
% mask = effect_img;
dat = [];
dat_stat = [];

dat = fmri_data(effect_img, mask);
datp = fmri_data(pval_img, mask);
datp.dat(datp.dat == 0) = Inf;

dat_stat = statistic_image;
dat_stat.dat = dat.dat;
dat_stat.p = datp.dat;
dat_stat.volInfo = dat.volInfo;
%FDR(dat_stat.p,0.05)

r = pruning_img(dat_stat,[.05 .01 FDR(dat_stat.p,0.05)], [1 1 10]);
r = region(r);
r(cat(3,r(:).numVox)<20) = [];
%orthview(r);
brain_activations_wani(r,'all2','pruned');
%%
%dat_stat = threshold(dat_stat, .005, 'unc', 'k',5); % if sec_i =1
dat_stat = threshold(dat_stat, .05, 'fdr', 'k',10);
%dat_stat = threshold(dat_stat, .01, 'extent', dat,'k',10);
final_dat = double(dat_stat.sig) .* sign(dat_stat.dat);

dat_stat = threshold(dat_stat, .01, 'unc', 'k',1);
final_dat = final_dat + double(dat_stat.sig) .* sign(dat_stat.dat);

dat_stat = threshold(dat_stat, .05, 'unc', 'k',1);
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
r(cat(1,r(:).numVox) < 20) = [];
orthviews(r);
%dat2 = region2imagevec(r);
%dat1 = region2imagevec(r);
%dat1 = replace_empty(dat1);
%dat2 = replace_empty(dat2);
%dat1.dat(:,2) = dat2.dat;
%orthviews(dat1);

%%
brain_activations_wani(r,'all2','pruned');

