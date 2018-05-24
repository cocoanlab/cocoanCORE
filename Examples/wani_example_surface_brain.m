close all;
maskdir = '/Users/clinpsywoo/Dropbox/Resources/Wagerlab_Single_Trial_Pain_Datasets/wani_results/v11_controlfor_temp_nps_scale_pattern_later/ttest_v11_4_wo_scebl';
cl = region(fullfile(maskdir, 'nonnoc_v11_4_137subjmap_fdr05_unc0025_wttest.nii'));

poscm = colormap_tor([0.96 0.41 0], [1 1 0]);  % warm
negcm = colormap_tor([0.11 0.46 1], [.23 1 1]);  % cools

% h = add_surface(which('surf_BrainMesh_ICBM152.mat'));
% set(h, 'FaceColor', [.7 .7 .7], 'FaceAlpha', .1);


cluster_surf(cl ,'fsavg_right', 2, 'heatmap', 'colormaps', poscm, negcm)
% cluster_surf(cl ,which('surf_freesurf_inflated_Right.mat'), 2, 'heatmap', 'colormaps', poscm, negcm)
% cluster_surf(cl ,which('surf_workbench_very_inflated_32k_Left.mat'), 2, 'heatmap', 'colormaps', poscm, negcm)
%cluster_surf(cl ,which('surf_workbench_very_inflated_32k_Right.mat'), 2, 'heatmap', 'colormaps', poscm, negcm)
%cluster_surf(cl ,which('surf_workbench_inflated_32k_Left.mat'), 2, 'heatmap', 'colormaps', poscm, negcm)
%cluster_surf(cl ,which('surf_workbench_inflated_32k_Right.mat'), 2, 'heatmap', 'colormaps', poscm, negcm)

%cluster_surf(cl ,which('surf_BrainMesh_ICBM152Right_smoothed.mat'), 2, 'heatmap', 'colormaps', poscm, negcm)
%cluster_surf(clusters1,'surf_parahippo_havardoxford_20_l.mat', 2, 'heatmap', 'colormaps', poscm, negcm)

h = get(gca, 'children');
% set(h(2), 'FaceAlpha', .5);
% set(h(2), 'FaceAlpha', 1);
set(h(2), 'AmbientStrength', .5)
% set(h(2), 'AmbientStrength', 1)

axis vis3d;


%% Make a new mesh using atlas: parahippo

maskdir = '/Users/clinpsywoo/Documents/Workspace/Wagerlab_Single_Trial_Pain_Datasets/wani_results/v11_controlfor_temp_nps_scale_pattern_later/ttest_v11_4_wo_scebl';
cd(maskdir);
P = fullfile(maskdir, 'parahippo_havardoxford_20.nii');

close all;
create_figure('Brain Surface'); 

[hPatch,outP,FV, cl, myLight] = mask2surface(P);

v1 = [FV(1).vertices];
f1 = [FV(1).faces];

[sv1, sf1] = smoothMesh(v1, f1, 3);

vertices = sv1;
faces = sf1;

save('surf_parahippo_havardoxford_20_r.mat', 'vertices', 'faces');


v2 = [FV(2).vertices];
f2 = [FV(2).faces];

[sv2, sf2] = smoothMesh(v2, f2, 3);

vertices = sv2;
faces = sf2;

save('surf_parahippo_havardoxford_20_l.mat', 'vertices', 'faces');

%% Draw the anatomy: NAC

maskdir = '/Users/clinpsywoo/Documents/Workspace/Wagerlab_Single_Trial_Pain_Datasets/wani_results/v11_controlfor_temp_nps_scale_pattern_later/ttest_v11_4_wo_scebl';
cd(maskdir);
P = fullfile(maskdir, 'accumbens_havardoxford_20.nii');

close all;
create_figure('Brain Surface'); 

[hPatch, outP, FV, cl, myLight] = mask2surface(P);

v1 = [FV(1).vertices];
f1 = [FV(1).faces];

[sv1, sf1] = smoothMesh(v1, f1, 3);

vertices = sv1;
faces = sf1;

save('surf_accumbens_havardoxford_20_r.mat', 'vertices', 'faces');


v2 = [FV(2).vertices];
f2 = [FV(2).faces];

[sv2, sf2] = smoothMesh(v2, f2, 3);

vertices = sv2;
faces = sf2;

save('surf_accumbens_havardoxford_20_l.mat', 'vertices', 'faces');

%% DRAW
close all;
maskdir = '/Users/clinpsywoo/Documents/Workspace/Wagerlab_Single_Trial_Pain_Datasets/wani_results/v11_controlfor_temp_nps_scale_pattern_later/ttest_v11_4_wo_scebl';
clusters1 = region(fullfile(maskdir, 'parahippo_havardoxford_20_SIIPS1.nii'));
clusters2 = region(fullfile(maskdir, 'accumbens_havardoxford_20_SIIPS1.nii'));

poscm = colormap_tor([0.96 0.41 0], [1 1 0]);  %slate to orange to yellow
negcm = colormap_tor([0.11 0.46 1], [.23 1 1]);  % light blue to dark blue

% h = add_surface(which('surf_BrainMesh_ICBM152.mat'));
% set(h, 'FaceColor', [.7 .7 .7], 'FaceAlpha', .1);

cluster_surf(clusters1,'surf_parahippo_havardoxford_20_r.mat', 2, 'heatmap', 'colormaps', poscm, negcm)
cluster_surf(clusters1,'surf_parahippo_havardoxford_20_l.mat', 2, 'heatmap', 'colormaps', poscm, negcm)

cluster_surf(clusters2,'surf_accumbens_havardoxford_20_r.mat', 2, 'heatmap', 'colormaps', poscm, negcm)
cluster_surf(clusters2,'surf_accumbens_havardoxford_20_l.mat', 2, 'heatmap', 'colormaps', poscm, negcm)

lighting gouraud;
material shiny
% camlight right
%camlight headlight
view(-197, 33);

a = get(gca, 'children');
for i = 1:numel(a)
    try
        set(a(i), 'ambientStrength', .5);
    catch

    end
end

% %% save freesurfer surfs
% 
% [vertices, faces] = read_surf('/Users/clinpsywoo/Dropbox/MATLAB/Connectivity_toolboxes/conn/utils/surf/lh.inflated.surf');
% faces = faces + 1; 
% vertices = vertices + 1;
% save('/Users/clinpsywoo/Nas/Resources/github_nas/cocoanlab/cocoanCORE/Canonical_brains/surf_freesurf_inflated_Left.mat', 'vertices', 'faces');
% 
% [vertices, faces] = read_surf('/Users/clinpsywoo/Dropbox/MATLAB/Connectivity_toolboxes/conn/utils/surf/rh.inflated.surf');
% faces = faces + 1; 
% vertices = vertices + 1;
% save('/Users/clinpsywoo/Nas/Resources/github_nas/cocoanlab/cocoanCORE/Canonical_brains/surf_freesurf_inflated_Right.mat', 'vertices', 'faces');
