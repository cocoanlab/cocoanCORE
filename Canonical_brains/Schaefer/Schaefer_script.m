[rootdir, basedir, gitdir] = set_path_env('ncpu', 1);

%% Schaefer265 - add subcortex, brainstem, cerebellum

Schaefer_img = spm_vol(which('Schaefer2018_200Parcels_7Networks_order_FSLMNI152_2mm.nii'));
Schaefer_dat = spm_read_vols(Schaefer_img);
Brainnetome_img = spm_vol(which('BN_Atlas_275_combined_2mm_MNI_brainstem.nii'));
Brainnetome_dat = spm_read_vols(Brainnetome_img);

Schaefer_new_img = Schaefer_img;
Schaefer_new_img.fname = strrep(Schaefer_img.fname, 'Schaefer2018_200Parcels_7Networks_order_FSLMNI152_2mm.nii', 'Schaefer_265_combined_2mm.nii');

Schaefer_new_dat = Schaefer_dat;
Schaefer_new_idx = max(Schaefer_dat(:));
for reg_i = 211:275
    Schaefer_new_idx = Schaefer_new_idx + 1;
    wh_add = Brainnetome_dat == reg_i & Schaefer_new_dat == 0;
    Schaefer_new_dat(wh_add) = Schaefer_new_idx;
end

spm_write_vol(Schaefer_new_img, Schaefer_new_dat);

%% Schaefer265 - resample to 4mm

system(['export FSLOUTPUTTYPE=NIFTI;' ...
    ...
    'flirt' ...
    ' -in ' Schaefer_new_img.fname ...
    ' -ref ' Schaefer_new_img.fname ...
    ' -applyisoxfm 4' ...
    ' -interp nearestneighbour' ...
    ' -out ' strrep(Schaefer_new_img.fname, 'Schaefer_265_combined_2mm.nii', 'Schaefer_265_combined_4mm.nii')]);

%% Schaefer265 - metadata

load(which('cluster_Fan_Net_r275.mat'));
Schaefer_label_dat = importdata(which('Schaefer2018_200Parcels_7Networks_order.txt'));
Schaefer_network_match = regexp(Schaefer_label_dat.textdata(:,2), '(?<=H_)[a-zA-Z]*(?=_)', 'match');
Schaefer_network_groups = {'Vis', 'SomMot', 'DorsAttn', 'SalVentAttn', 'Limbic', 'Cont', 'Default'};
for reg_i = 1:numel(Schaefer_network_match)
    Schaefer_network(reg_i, 1) = find(strcmp(Schaefer_network_match{reg_i}{1}, Schaefer_network_groups));
end
Schaefer_network = [Schaefer_network; cluster_Fan_Net.dat(211:275,9)];

Schaefer_Net_Labels.dat(:,1) = 1:265;
Schaefer_Net_Labels.dat(:,2) = Schaefer_network;
Schaefer_Net_Labels.dat(:,3) = zeros(265,1);
Schaefer_Net_Labels.dat(contains(Schaefer_Net_Labels.names, {'LH', 'Left', '_L_'}), 3) = -1;
Schaefer_Net_Labels.dat(contains(Schaefer_Net_Labels.names, {'RH', 'Right', '_R_'}), 3) = 1;
Schaefer_Net_Labels.descrip = {'Column 1: Original mask value, Column 2: Network value', ...
    '1: visual (occipital), 2: somatomotor, 3: dAttention, 4: vAttention, 5: Limbic, 6: Frontoparietal 7: Default, 8: Subcortical, 9:Cerebellum , 10:Brainstem', ...
    'Laterality. -1: Left, 1: Right, 0: No laterality'};
Schaefer_Net_Labels.ten_network_col = cluster_Fan_Net.ten_network_col;
Schaefer_Net_Labels.r_2mm = region(which('Schaefer_265_combined_2mm.nii'), 'unique_mask_values');
Schaefer_Net_Labels.names = [strrep(Schaefer_label_dat.textdata(:,2), '7Networks_', ''); cluster_Fan_Net.names(211:275)]; 

save(strrep(which('Schaefer_265_combined_2mm.nii'), 'Schaefer_265_combined_2mm.nii', 'Schaefer_Net_Labels_r265.mat'), 'Schaefer_Net_Labels');

%% Schaefer100 - metadata

load(which('cluster_Fan_Net_r275.mat'));
Schaefer_label_dat = importdata(which('Schaefer2018_100Parcels_7Networks_order.txt'));
Schaefer_network_match = regexp(Schaefer_label_dat.textdata(:,2), '(?<=H_)[a-zA-Z]*(?=_)', 'match');
Schaefer_network_groups = {'Vis', 'SomMot', 'DorsAttn', 'SalVentAttn', 'Limbic', 'Cont', 'Default'};
for reg_i = 1:numel(Schaefer_network_match)
    Schaefer_network(reg_i, 1) = find(strcmp(Schaefer_network_match{reg_i}{1}, Schaefer_network_groups));
end

Schaefer_Net_Labels.dat(:,1) = 1:100;
Schaefer_Net_Labels.dat(:,2) = Schaefer_network;
Schaefer_Net_Labels.names = strrep(Schaefer_label_dat.textdata(:,2), '7Networks_', '');
Schaefer_Net_Labels.dat(:,3) = zeros(100,1);
Schaefer_Net_Labels.dat(contains(Schaefer_Net_Labels.names, {'LH', 'Left', '_L_'}), 3) = -1;
Schaefer_Net_Labels.dat(contains(Schaefer_Net_Labels.names, {'RH', 'Right', '_R_'}), 3) = 1;
Schaefer_Net_Labels.descrip = {'Column 1: Original mask value, Column 2: Network value', ...
    '1: visual (occipital), 2: somatomotor, 3: dAttention, 4: vAttention, 5: Limbic, 6: Frontoparietal 7: Default, 8: Subcortical, 9:Cerebellum , 10:Brainstem', ...
    'Laterality. -1: Left, 1: Right, 0: No laterality'};
Schaefer_Net_Labels.seven_network_col = cluster_Fan_Net.ten_network_col(1:7,:);
Schaefer_Net_Labels.r_2mm = region(which('Schaefer2018_100Parcels_7Networks_order_FSLMNI152_2mm.nii'), 'unique_mask_values');

save(strrep(which('Schaefer_Net_Labels_r265.mat'), 'Schaefer_Net_Labels_r265.mat', 'Schaefer_Net_Labels_r100.mat'), 'Schaefer_Net_Labels');


