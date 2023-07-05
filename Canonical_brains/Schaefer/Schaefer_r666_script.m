%% Schaefer666 - combine with subcortex r66
gitdir = '/Volumes/homeo/github';
spmdir = '/Volumes/homeo/dropbox/resources/spm12';
addpath(genpath(spmdir));

atlas_subcortex = spm_read_vols(spm_vol(fullfile(gitdir, 'cocoanlab/cocoanCORE/Canonical_brains/subcortex_comb/subcortex_comb_r66.nii')));
atlas_cortex = spm_read_vols(spm_vol(fullfile(gitdir, 'cocoanlab/cocoanCORE/Canonical_brains/Schaefer/Schaefer2018_600Parcels_7Networks_order_FSLMNI152_2mm.nii')));

atlas_combined = atlas_cortex;
atlas_combined(atlas_subcortex > 0) = ...
    atlas_subcortex(atlas_subcortex > 0) + max(atlas_combined, [], 1:3);

atlas_combined_vol = spm_vol(fullfile(gitdir, 'cocoanlab/cocoanCORE/Canonical_brains/subcortex_comb/subcortex_comb_r66.nii'));
atlas_combined_vol.fname = fullfile(gitdir, 'cocoanlab/cocoanCORE/Canonical_brains/Schaefer/Schaefer_666_combined_2mm.nii');
spm_write_vol(atlas_combined_vol, atlas_combined);

rmpath(genpath(spmdir));

%% Schaefer666 - metadata

load(fullfile(gitdir, 'cocoanlab/cocoanCORE/Canonical_brains/subcortex_comb/subcortex_comb_r66_labels.mat'));
Schaefer_label_dat = importdata(fullfile(gitdir, 'cocoanlab/cocoanCORE/Canonical_brains/Schaefer/Schaefer2018_600Parcels_7Networks_order.txt'));
Schaefer_network_match = regexp(Schaefer_label_dat.textdata(:,2), '(?<=H_)[a-zA-Z]*(?=_)', 'match');
Schaefer_network_groups = {'Vis', 'SomMot', 'DorsAttn', 'SalVentAttn', 'Limbic', 'Cont', 'Default'};
Schaefer_network = [];
for reg_i = 1:numel(Schaefer_network_match)
    Schaefer_network(reg_i, 1) = find(strcmp(Schaefer_network_match{reg_i}{1}, Schaefer_network_groups));
end
Schaefer_network = [Schaefer_network; subcortex_comb_r66.dat(:,2) + max(Schaefer_network)];

Schaefer_Net_Labels.dat(:,1) = 1:numel(Schaefer_network);
Schaefer_Net_Labels.dat(:,2) = Schaefer_network;
Schaefer_Net_Labels.names = [strrep(Schaefer_label_dat.textdata(:,2), '7Networks_', ''); subcortex_comb_r66.names]; 
Schaefer_Net_Labels.dat(:,3) = zeros(numel(Schaefer_network),1);
Schaefer_Net_Labels.dat(contains(Schaefer_Net_Labels.names, {'LH', 'Left', '_L_'}), 3) = -1;
Schaefer_Net_Labels.dat(contains(Schaefer_Net_Labels.names, {'RH', 'Right', '_R_'}), 3) = 1;
Schaefer_Net_Labels.descrip = {'Column 1: Original mask value, Column 2: Network value, Column 3: Laterality', ...
    '1: visual (occipital), 2: somatomotor, 3: dAttention, 4: vAttention, 5: Limbic, 6: Frontoparietal 7: Default, 8: Thalamus, 9: Hippocampus/Amygdala, 10: Basal ganglia, 11: Cerebellum, 12: Brainstem', ...
    'Laterality. -1: Left, 1: Right, 0: No laterality'};

save(fullfile(gitdir, 'cocoanlab/cocoanCORE/Canonical_brains/Schaefer/Schaefer_Net_Labels_r666.mat'), 'Schaefer_Net_Labels');
