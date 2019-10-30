%% make a mask for MPFC

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
cloutdat = region2imagevec(clout);
cloutdat.dat(cloutdat.dat~=0) = ones(sum(cloutdat.dat~=0),1);
outputdir = '/Volumes/RAID1/labdata/current/Metaanalysis_Anjali/Anjali_MPFC_subcortical_connectivity/data/ROI_masks';
cloutdat.fullpath = fullfile(outputdir, 'MPFC_mask.nii');
write(cloutdat)

%% make a mask for OFC
cluster_orthviews(cl, 'unique');

clout = [];

% choose OFC clusters
[clout,cl] = cluster_graphic_select(cl,clout);

% display
cluster_orthviews(clout, 'unique');

for i =7:length(clout)
    k = (clout(i).XYZmm(3,:) > -13.5);
    clout(i).XYZ(:,k) = [];
    clout(i).XYZmm(:,k) = [];
    clout(i).val(k) = [];
    clout(i).Z(k) = [];
    clout(i).numVox = length(clout(i).Z);
end

% save
outputdir = '/Volumes/RAID1/labdata/current/Metaanalysis_Anjali/Anjali_MPFC_subcortical_connectivity/data/ROI_masks';
filename = 'OFC_mask.img';
wani_make_mask(outputdir, filename, clout);

%% make a mask for VLFC
cluster_orthviews(cl, 'unique');

clout = [];

% choose VLPFC clusters
[clout,cl] = cluster_graphic_select(cl,clout);

% display
cluster_orthviews(clout, 'unique');

% save
outputdir = '/Volumes/RAID1/labdata/current/Metaanalysis_Anjali/Anjali_MPFC_subcortical_connectivity/data/ROI_masks';
filename = 'VLPFC_mask.img';
wani_make_mask(outputdir, filename, clout);

%% make a mask for subcortex
cl = region(which('atlas_labels_combined.img'), 'unique_mask_values'); %anat_lbpa_thal.img'); %% 'lpba40.spm5.avg152T1.label.nii');

clout = cl([1:7 56:91]);

% display
cluster_orthviews(clout, 'unique');

% save
outputdir = '/Volumes/RAID1/labdata/current/Metaanalysis_Anjali/Anjali_MPFC_subcortical_connectivity/data/ROI_masks';
filename = 'subcortex_mask.img';
wani_make_mask(outputdir, filename, clout);


%% make a mask for cortex
cl = region(which('atlas_labels_combined.img'), 'unique_mask_values'); %anat_lbpa_thal.img'); %% 'lpba40.spm5.avg152T1.label.nii');

clout = cl([8:55]);

% display
cluster_orthviews(clout, 'unique');

% save
outputdir = '/Volumes/RAID1/labdata/current/Metaanalysis_Anjali/Anjali_MPFC_subcortical_connectivity/data/ROI_masks';
filename = 'cortex_mask.img';
wani_make_mask(outputdir, filename, clout);

%% make a frontal mask
cl = region(which('atlas_labels_combined.img'), 'unique_mask_values'); %anat_lbpa_thal.img'); %% 'lpba40.spm5.avg152T1.label.nii');
clout = cl(14:15);
cl_frontal = region ('/Volumes/RAID1/labdata/current/Metaanalysis_Anjali/Anjali_MPFC_subcortical_connectivity/frontal_region/frontal_mask_image.img');

cl_frontal(2:3) = clout(1:2);

outputdir = '/Volumes/RAID1/labdata/current/Metaanalysis_Anjali/Anjali_MPFC_subcortical_connectivity/data/ROI_masks';
filename = 'frontal_mask.img';
wani_make_mask(outputdir, filename, cl_frontal);


%% make mask data

clear;
outputdir = '/Volumes/RAID1/labdata/current/Metaanalysis_Anjali/Anjali_MPFC_subcortical_connectivity/data/ROI_masks';

cd(outputdir);

masks = filenames('*.img');

for i = 4 %1:length(masks)
    [~, names{i}, ~] = fileparts(masks{i});

    eval(['dat = fmri_data(''' names{i} '.img'');']);
    data = zeros(size(dat.volInfo.image_indx),1);
    data(dat.volInfo.wh_inmask) = dat.dat;

    fileID = fopen([names{i} '_wholebrain.txt'], 'w');
    fprintf(fileID, '%d\n', data);
    fclose(fileID);
end


%% get study by roi (I don't need this)
clear;
outputdir = '/Volumes/RAID1/labdata/current/Metaanalysis_Anjali/Anjali_MPFC_subcortical_connectivity/data/ROI_masks';

cd(outputdir);

masks = filenames('*.img');

for i = 1:3
    cl{i} = region(masks{i}, 'unique_mask_values');
    [~, names{i}, ~] = fileparts(masks{i});
end

dosave = 0;

get_studies_using_cl(cl, names, outputdir, dosave)


%% write a text file

load studybycl;

for i = 1:3
    fileID = fopen([names{i} '.txt'], 'w');
    fprintf(fileID, '%d\n', studybyroi(:,i));
    fclose(fileID);
end