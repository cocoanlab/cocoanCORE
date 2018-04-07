function cluster_table_from_image(effectimg, prim_p, k, subj_n, varargin)

% function cluster_table_from_image(beta_con_t_img, prim_p, k, subj_n, varargin)
%
% This function will give provide a cluster table that consists of
% anatomical labels (using canlab atlas and AAL atlas), XYZmm, direction
% (activation or deactivation), max statistics (including T/Z value,-logP), 
% the number of voxels, and the volume(mm^3) of clusters.
% 
% To use this, you need AAL atlas data (you can download from the 
% website - http://www.cyceron.fr/web/aal__anatomical_automatic_labeling.html) 
% and canlab atlas data. Also you need a function "cluster_table_wani.m",
% which is modified from the function, "cluster_table.m". 
% 
% input: 
%   beta_t_con_img: beta or t or contrast image. If you use beta or contrast
%       image, you need to add p image in varargin.  
%   prim_p and k: parameters for cluster-extent thresholding. 
%       prim_p: a primary p value - uncorrected. 
%       k: cluster extent size
%     if you don't want to use extent-based thresholding, k = 1. 
%   subj_n: the number of subjects
%   varargin:
%       'subcluster': provide subcluster info
%       'dosave' or 'save': save text file. you need to provide a name of savefile.
%       'mask': use mask. you need to provide the name of mask image. 
%           default mask is 'brainmask.nii', which is under spm8/apriori
%           directory
%       'pimg': If you are using beta or contast image, you need to give
%           the function p image as well. {'pimg', pimg}
%
% Example:
% 
% % Example for robust regression output files
% beta_or_con_img = 'rob_beta_0001.img'; 
% pimg = 'rob_p_0001.img'; 
% mask = which('mask_colin27.img');
% savefile = fullfile(outputdir, 'test.txt');
% prim_p = .005; k = 75; subj_n = 33;
%
% %1. use beta or contrast image, do not save, do not use subcluster
%     cluster_table_from_image(beta_or_con_img, prim_p, k, subj_n, 'mask', mask, 'pimg', pimg)
% %2. use beta or contrast image, save cluster table 
%     cluster_table_from_image(beta_or_con_img, prim_p, k, subj_n, 'mask', mask, 'pimg', primg, 'save', savefile)
%
% % Example for SPM output files
% t_img = 'spmT_0001.img';
% mask = which('mask_colin27.img');
% savefile = fullfile(outputdir, 'test.txt');
% prim_p = .005; k = 75; subj_n = 33;
%
% %3. use t image, save cluster table and report subclusters
%     cluster_table_from_image(t_img, prim_p, k, subj_n, 'subcluster', 'mask', mask, 'save', savefile)
%

subcluster = 0;
dosave = 0;

try 
    mask = which('brainmask.nii');
catch
    mask = effectimg;
end

use_timg = 1;

if ~isempty(varargin)
    for i = 1:length(varargin)
        if ischar(varargin{i})
            switch varargin{i}
                case 'subcluster', subcluster = 1;
                case {'dosave', 'save'}, dosave = 1; savefilename = varargin{i+1};
                case 'mask', mask = varargin{i+1};
                case 'pimg', use_timg = 0; pimg = varargin{i+1};
            end
        end
    end
end

data = fmri_data(mask, effectimg);

dat = statistic_image('image_names', effectimg, 'type', 'beta');

if use_timg
    datp = dat.dat;
    datp = tcdf(-abs(datp), subj_n-1);
else
    [datp, volInfo] = iimg_threshold(pimg, 'mask', effectimg);
end 
    
datp(datp == 0 & isnan(dat.dat)) = NaN;
dat.p = datp;
dat = threshold(dat, 0.01, 'unc', 'k', 1);

% mask 
dat.dat(~logical(data.dat)) = 0;
dat.sig(~logical(data.dat)) = 0;
dat.p(~logical(data.dat)) = 1;

% threshold 
dat = replace_empty(dat);
dat = threshold(dat, prim_p, 'unc', 'k', k);
dat1 = dat;
dat1.dat = -log10(dat1.p) .* dat1.sig .* (-(dat1.dat < 0) + (dat1.dat > 0));
dat1.fullpath = fullfile(pwd, 'temp_img_for_table.img');
write(dat1);
cl = region(dat1.fullpath);

if dosave
    cluster_table_wani(cl, subcluster, 0, [1 subj_n-1], 'writefile', savefilename);
else
    cluster_table_wani(cl, subcluster, 0, [1 subj_n-1]);
end

delete(dat1.fullpath);

end