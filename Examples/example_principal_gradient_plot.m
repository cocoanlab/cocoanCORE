% Exmaple scripts for visualization of principal_gradient_axis 
% using "principal_gradient_plot" 
%
% Suhwan Gim 
% 2020. Oct. 20

%% 0) path
addpath(genpath('/Volumes/sein/github/cocoanlab/cocoanCORE'));

%% 1) visualization of first gradient map of the principal gradient axis made by BrainSpace toolbox

%gradient_mask = fmri_data(which('ten_prinicpal_gradients_volumn_DE.nii')); % in cocoanCORE/Canocial_brains (n=56)
gradient_mask = fmri_data( which('Volumetric_hcp_gradients_GSP_90_DE_wholebrain_2mm.nii')); % in cocoanCORE/Canocial_brains (n=1000)
first_gradient_mask = gradient_mask.get_wh_image(1); %1~10
first_gradient_mask.dat = first_gradient_mask.dat - min(first_gradient_mask.dat); % just for visualization 
% color
cmap = csvread(which('colormap_rainbow_gradient.csv')); 
cmap = cmap(:,1:3);
% visualization 
depth = 3;
surface_name_L = 'fsavg_left';
create_figure;
cluster_surf(region(first_gradient_mask), depth, surface_name_L,'colormaps', flip(cmap), [], 'heatmap'); % canlabCore
out.h = get(gca, 'children');
set(out.h(2), 'BackFaceLighting', 'lit')
camlight(-90,-20);
axis vis3d;
view(-90, 0);

%% 2) analyzing own threhsed images on the first gradient map
help principal_gradients_plot;
thresh_img = {'thresh_img_005.nii'};  % example
results = principal_gradients_plot(thresh_img);
disp(results)
% 2-1) save figure
pagesetup(gcf);
savenames = 'temp_map.pdf';
saveas(gcf, savenames);

%% 3) with several images 
thresh_img = {'thresh_img_005.nii','thresh_img_005.nii'};  % example
results = principal_gradients_plot(thresh_img);
%% 4) existing figures
results = principal_gradients_plot(thresh_img,'samefigure');
%% 5) Use other component with different number of bins
numComp = 2; % second components (1~10)
results = principal_gradients_plot(thresh_img,'other_comp',numComp,'numBins',10); 



