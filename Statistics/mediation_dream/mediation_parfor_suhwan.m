%%
basedir = '/Volumes/homeo/Temp_Suhwan/SEMIC';
%%
load('/Volumes/homeo/Temp_Suhwan/SEMIC/SEMIC_231020_other_mediation_med_vars_overall1secs.mat','med_vars')
%load(fullfile('/Users/suhwan/Dropbox/Projects/SEMIC/sync_NAS/sync/data/meditaion_2mm_OCT2023/SEMIC_231019_self_mediation_med_vars_overall1secs.mat'),'med_vars');
%% ADDpath
addpath(genpath('/Users/suhwan/Dropbox/github/CanlabTools'))
addpath(genpath('/Users/suhwan/Dropbox/github/CocoanTools/cocoanCORE'))
addpath(genpath('/Users/suhwan/Dropbox/github/external_toolbox/spm12'))
rmpath(genpath('/Users/suhwan/Dropbox/github/external_toolbox/spm12/external/'));
rmpath(genpath('/Users/suhwan/Dropbox/github/external_toolbox/spm12/external/fieldtrip/compat'));
addpath(genpath('/Users/suhwan/Dropbox/Projects/SEMIC2/scripts/mediation_dream_toolbox'))
rmpath(genpath('/Users/suhwan/Dropbox/github/CanlabTools/CanlabPrivate'))
%%
mask = fullfile(basedir, 'gm_mask_semic2.nii');
for sub_i = 1:58
    temp_imgs=strrep(med_vars.imgs{sub_i}, '/Volumes/sein/dropbox/projects/SEMIC/analysis/imaging/first_level',basedir);
    dat = fmri_data(temp_imgs,mask);    
    for vox_i = 1:length(dat.dat)
        med_vars.m1{vox_i}{sub_i} = dat.dat(vox_i,:)';% voxel by trials
    end
end

% %%
% tas = [];
% 
% 
% for vox_i = 1:length(dat.dat)
%     tas{vox_i}=cellfun(@transpose, med_vars.m1{vox_i},'UniformOutput',false);
%     disp(vox_i);
% end

%% runing two-path mediaion 
gend = []; 
[ gen_Apath,gen_Bpath,gen_Cppath,gen_Cpath,gen_ABpath, ...
    gen_Ap,gen_Bp,gen_Cpp,gen_Cp,gen_ABp] = deal(NaN(size(med_vars.m1,2),1));
tic;
parfor vox_i = 1:size(med_vars.m1,2)
    temp_gend = []; 
    [temp_gend2, temp_gend] = mediation(med_vars.X1, med_vars.Y, med_vars.m1{vox_i}', 'names', {'STIM' 'ratings' 'brain' },'boottop','bootsamples',10000,'covs', med_vars.X2);
    %[temp_gend2, temp_gend] = mediation(med_vars.X1, med_vars.Y, tas{vox_i}', 'names', {'STIM' 'ratings' 'brain' },'boottop','bootsamples',10000,'covs', med_vars.X2);
    
    gen_Ap(vox_i,:)=temp_gend.p(1);
    gen_Bp(vox_i,:)=temp_gend.p(2);
    gen_Cpp(vox_i,:)=temp_gend.p(3);
    gen_Cp(vox_i,:)=temp_gend.p(4);
    gen_ABp(vox_i,:)=temp_gend.p(5);
    
    gen_Apath(vox_i,:)=temp_gend.beta(1);
    gen_Bpath(vox_i,:)=temp_gend.beta(2);
    gen_Cppath(vox_i,:)=temp_gend.beta(3);
    gen_Cpath(vox_i,:)=temp_gend.beta(4);
    gen_ABpath(vox_i,:)=temp_gend.beta(5);
   disp(vox_i); 
end
toc;
med_res.gen_Ap = gen_Ap; %P value of path A
med_res.gen_Bp = gen_Bp; %P value of path A
med_res.gen_Cpp = gen_Cpp; %P value of path A
med_res.gen_Cp = gen_Cp; %P value of path A
med_res.gen_ABp = gen_ABp; %P value of path A

med_res.gen_Apath = gen_Apath; %P value of path A
med_res.gen_Bpath = gen_Bpath; %P value of path A
med_res.gen_Cppath = gen_Cppath; %P value of path A
med_res.gen_Cpath = gen_Cpath; %P value of path A
med_res.gen_ABpath = gen_ABpath; %P value of path A
%
med_res.maks = fmri_data(mask,mask);

% med_res.input_imagesc = conlists;
% med_res.desc = {'x = gender (female) , y= intercept of pain ratings, m = brain (contrastimages)'};

%%
%save(fullfile(basedir,'results','SEMIC_2mm_mediation_parfor_xSTIM.mat'),'med_res');
%save(fullfile(basedir,'results','SEMIC_2mm_mediation_parfor_OTHER_xSTIM.mat'),'med_res');

%%
fdrval = FDR([med_res.gen_ABp med_res.gen_Bp med_res.gen_Ap],0.05)
%%
% path ab
mask = med_res.maks;

dat = fmri_data(mask, mask);
datp = dat;

dat.dat = med_res.gen_ABpath;
datp.dat = med_res.gen_ABp;
% datp = apply_mask(datp, which('gray_matter_mask.nii'));
% dat = apply_mask(dat, which('gray_matter_mask.nii'));
dat.dat = dat.dat.*double(datp.dat < 0.05);
dat.dat = dat.dat.*double(datp.dat < FDR(datp.dat,0.05));
dat.dat = dat.dat.*double(datp.dat < fdrval);
orthviews(dat);
%brain_activations_display(region(dat));
%% PRUNING 
fdr_cor = fdrval;

clear stat_img 
stat_img = statistic_image;
stat_img.dat = med_res.gen_ABpath;
stat_img.p = med_res.gen_ABp;
stat_img.volInfo = mask.volInfo;
%cue_tstat = threshold(stat_img,0.05,'fdr','k',10);
stim_tstat = threshold(stat_img,fdrval3,'unc','k',5);
prun_stim = pruning_img(stat_img, [0.05 0.01 fdrval3],[1 1 10]);
brain_activations_display(region(stim_tstat))

%% path b1 and b2
dat = fmri_data(mask, mask);
datp = dat;
datp.dat = med_res.gen_Bp
dat.dat =med_res.gen_Bpath
%datp=apply_mask(datp, which('gray_matter_mask.nii'));
%dat=apply_mask(dat, which('gray_matter_mask.nii'));
dat.dat = dat.dat.*double(datp.dat <0.05);
dat.dat = dat.dat.*double(datp.dat <FDR(datp.dat,0.05));
orthviews(dat);

%% runing two-path mediaion 

[ cue_Apath,cue_Bpath,cue_Cppath,cue_Cpath,cue_ABpath, ...
    cue_Ap,cue_Bp,cue_Cpp,cue_Cp,cue_ABp] = deal(NaN(size(med_vars.m1,2),1));
tic;
parfor vox_i = 1:size(med_vars.m1,2)
    temp_gend = []; 
    [temp_gend2, temp_gend] = mediation(med_vars.X2, med_vars.Y, med_vars.m1{vox_i}', 'names', {'STIM' 'ratings' 'brain' },'boottop','bootsamples',10000,'covs', med_vars.X1);
    %[temp_gend2, temp_gend] = mediation(med_vars.X1, med_vars.Y, tas{vox_i}', 'names', {'STIM' 'ratings' 'brain' },'boottop','bootsamples',10000,'covs', med_vars.X2);
    
        cue_Ap(vox_i,:)=temp_gend.p(1);
    cue_Bp(vox_i,:)=temp_gend.p(2);
    cue_Cpp(vox_i,:)=temp_gend.p(3);
     cue_Cp(vox_i,:)=temp_gend.p(4);
    cue_ABp(vox_i,:)=temp_gend.p(5);
    
         cue_Apath(vox_i,:)=temp_gend.beta(1);
    cue_Bpath(vox_i,:)=temp_gend.beta(2);
    cue_Cppath(vox_i,:)=temp_gend.beta(3);
     cue_Cpath(vox_i,:)=temp_gend.beta(4);
    cue_ABpath(vox_i,:)=temp_gend.beta(5);
   disp(vox_i); 
end
toc;
med_res2.cue_Ap = cue_Ap; %P value of path A
med_res2.cue_Bp = cue_Bp; %P value of path A
med_res2.cue_Cpp = cue_Cpp; %P value of path A
med_res2.cue_Cp = cue_Cp; %P value of path A
med_res2.cue_ABp = cue_ABp; %P value of path A

med_res2.cue_Apath = cue_Apath; %P value of path A
med_res2.cue_Bpath = cue_Bpath; %P value of path A
med_res2.cue_Cppath = cue_Cppath; %P value of path A
med_res2.cue_Cpath = cue_Cpath; %P value of path A
med_res2.cue_ABpath = cue_ABpath; %P value of path A
med_res2.maks = fmri_data(mask,mask);

%%

%save(fullfile(basedir,'results','SEMIC_2mm_mediation_parfor_xCUE.mat'),'med_res2');
save(fullfile(basedir,'results','SEMIC_2mm_mediation_parfor_OTHER_xCUE.mat'),'med_res2');
%%
fdrval2 = FDR([med_res2.cue_ABp med_res2.cue_Bp med_res2.cue_Ap],0.05)

fdrval3 = FDR([med_res2.cue_ABp med_res2.cue_Bp med_res2.cue_Ap med_res.gen_ABp med_res.gen_Bp med_res.gen_Ap],0.05)
%%
% path ab
dat = fmri_data(mask, mask);
datp = dat;

datp.dat = med_res2.cue_ABp;
dat.dat = med_res2.cue_ABpath;
% datp = apply_mask(datp, which('gray_matter_mask.nii'));
% dat = apply_mask(dat, which('gray_matter_mask.nii'));
dat.dat = dat.dat.*double(datp.dat <0.05);
%dat.dat = dat.dat.*double(datp.dat <FDR(datp.dat,0.05));
dat.dat = dat.dat.*double(datp.dat <fdrval3);
orthviews(dat);

%% path b1 and b2
dat = fmri_data(mask, mask);
datp = dat;
datp.dat = med_res2.cue_Bp
dat.dat =med_res2.cue_Bpath
%datp=apply_mask(datp, which('gray_matter_mask.nii'));
%dat=apply_mask(dat, which('gray_matter_mask.nii'));
dat.dat = dat.dat.*double(datp.dat <0.05);
dat.dat = dat.dat.*double(datp.dat <FDR(datp.dat,0.05));
orthviews(dat);

%% PRUNING 

clear stat_img 
stat_img = statistic_image;
stat_img.dat = med_res2.cue_ABpath;
stat_img.p = med_res2.cue_ABp;
stat_img.volInfo = mask.volInfo;
%cue_tstat = threshold(stat_img,0.05,'fdr','k',10);
cue_tstat = threshold(stat_img,fdrval3,'unc','k',5);
prun_cue = pruning_img(stat_img, [0.05 0.01 fdrval3],[1 1 10]);