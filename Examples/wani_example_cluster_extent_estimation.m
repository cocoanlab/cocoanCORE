clear;

% outputdir = '/Volumes/RAID1/labdata/current/BMRK3/Imaging/Contrasts_Model15b';
outputdir = '/Volumes/current/BMRK3/Imaging/Contrasts_Model15b';
cd(outputdir);

% contrast results directory
d = dir;
for i = 1:length(d)
    names{i} = d(i).name;
end
names = names(2); % I'm only interested in the first contrast

% input variables
corrected_p = 0.05;
% prim_p = [.001 .0005 .0001 .00005 .00001];
prim_p = [.001 .0005 .0001 .00005 .00001 .000005 .000001];
voxelsize_mm = [2 2 2];
alphasim_dir = '/Users/clinpsywoo/abin/macosx_10.6_Intel_64'; % you have to give a directory that has alphasim (see help of get_cl_ext_FSL_Alphasim_SPM)

cd(fullfile(outputdir, names{1}));
load('data_obj.mat', 'imgs'); % load contrast images name
con_files = imgs;
con_files(:,9:22)=[];
% run
% [cl_ext, fwhm] = get_cl_ext_FSL_Alphasim_SPM(corrected_p, prim_p, rmm, con_files, outputdir);

[cl_ext, fwhm] = get_cl_ext_FSL_Alphasim_SPM(corrected_p, prim_p, rmm, voxelsize_mm, alphasim_dir, con_files, outputdir, 'twotail'); % do not run alphasim part
% [cl_ext, fwhm] = get_cl_ext_FSL_Alphasim_SPM(corrected_p, prim_p, rmm, con_files, outputdir, 'donotalphasim'); % do not run alphasim part

%% mediation analysis example

clear;
base_condir = '/Volumes/RAID1/labdata/current/BMRK3/Imaging/Imagination_mediation/mediation_Model_a1_boot_imgupdown_vw';
cd(base_condir);

con_files = filenames('X-M-Y_indiv_effect.img', 'absolute', 'char');
con_files = expand_4d_filenames(con_files, 33);

corrected_p = 0.05;
prim_p = [.001 0.0005 0.0001];
rmm = 2; 
voxelsize_mm = [2 2 2];
alphasim_dir = '/Users/clinpsywoo/abin/macosx_10.6_Intel_64'; % you have to give a directory that has alphasim (see help of get_cl_ext_FSL_Alphasim_SPM)

[cl_ext, fwhm] = get_cl_ext_FSL_Alphasim_SPM(corrected_p, prim_p, rmm, voxelsize_mm, alphasim_dir, con_files, base_condir, 'twotail');

%% cluster extent using smoothness from statistic image (not residual-based smoothness)
% this method is not recommended when there are residual images. However,
% sometimes we don't have residual images. In that case, we could use z or
% t statistic image to estimate smoothness level. 

clear;
% basedir = '/Volumes/RAID1/labdata/current/BMRK3';
basedir = '/Volumes/current/BMRK3';
resdir = fullfile(basedir, 'Imaging/threeway_wholebrain');

modeldir{1} = fullfile(resdir, 'model3_app_LNAC_WB_report');
modeldir{2} = fullfile(resdir, 'model2_app_WB_MPFC_report');
% modeldir{3} = fullfile(resdir, 'model1_app_BNAC_WB_report');

for i = 1:2
    cd(modeldir{i});

    path_you_want_to_display = 'X_M1_M2_Y';
    
%     p_image = [path_you_want_to_display '_pvals.img'];
%     beta_image = [path_you_want_to_display '_effect.img'];
    z_image = [path_you_want_to_display '_z.img'];
    
%     dat = fmri_data(p_image, p_image);
%     z_dat = spm_u(dat.dat/2, 32, 'Z'); % two-tailed
%     beta_dat = fmri_data(beta_image, beta_image);
%     z_sign = (beta_dat.dat > 0) + -1 * (beta_dat.dat < 0);
%     dat.dat = z_dat .* z_sign;
%     
    dat.fullpath = fullfile(modeldir{i}, z_image);
%     write(dat);

    mask = which('scalped_avg152T1_graymatter_smoothed.img');
    [dLh,resels,FWHM_mm, nVoxels, FWHM] = y_Smoothest_wani(dat.fullpath, mask, 32);
    
    mask = which('scalped_avg152T1_graymatter_smoothed.hdr');
    voxelsize_mm = [2 2 2];
    prim_p = [.005 .001 .0005 .0001 .00005 .00001 .000005 .000001];
    ClustSim_dir = '/Users/clinpsywoo/abin/macosx_10.6_Intel_64'; % you have to give a directory that has alphasim (see help of get_cl_ext_FSL_Alphasim_SPM)
    iter = 10000;
    corrected_p = .05;

    [fwhm{i}, cl_3dClustSim{i}] = cl_ext_3dClustSim(corrected_p, prim_p, dat.fullpath, voxelsize_mm, ClustSim_dir, 'iter', iter, 'fwhm', FWHM_mm', 'twotail', 'mask', mask); % two-tail - prim_p/2

end
