function dat = wani_make_mask(fname, cl, varargin)

% function wani_make_mask(fname, cl, varargin)
%
% This is a function to make a mask file. This function is using
% clusters2mask.m and 'brainmask.nii' as a overlay. 
% input: 
%    fname: filename and path
%       e.g.) outputdir = 'Volume/current/BMRK3';
%             fname = fullfile(outputdir, 'mask_DMPFC.img');
% varargin options:
%    'unique_mask_values': will save cl with different numbers
%    'original_mask_values': will saven cl with original values
%
%

dounique = 0;
doorigi = 0;
dobinary = 1;
domask = 0;
dowrite = 1;

for i = 1:numel(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            case 'unique_mask_values'
                dounique = 1;
                dobinary = 0;
            case 'original_mask_values'
                doorigi = 1;
                dobinary = 0;
            case 'mask'
                domask = 1;
                mask = varargin{i+1};
            case 'nowrite'
                dowrite = 0;
        end
    end
end

fprintf('\nStart creating mask:\n');

if dobinary

    % new method using region object
%     try
%         cldat = region2imagevec(cl);
%     catch
%         for i = 1:numel(cl), cl(i) = cluster2region(cl(i)); end
%         cldat = region2imagevec(cl);
%     end
%     cldat.dat(cldat.dat ~= 0) = 1;
%     cldat.fullpath = fname;
%     write(cldat);

% old method: this is much faster
    P = which('brainmask.nii');
    % P = which('weights_NSF_grouppred_cvpcr.img');
    V = spm_vol(P);
    V = V(1);
    V(1).private = [];
    V(1).fname = fname;

    clusters2mask(cl,V,0,fname);

elseif dounique
    fprintf('You chose to make a mask with unique mask values.\n');
    if domask
        if isa(mask, 'image_vector')
           dat = mask; clear mask;
           dat.dat = zeros(size(dat.dat));
        else
           error('not implemented yet');
        end
        
        for i = 1:numel(cl)
            fprintf('Working on cluster %04d among %04d clusters... ', i, numel(cl));
            cldat = region2imagevec(cl(i));
            dat.dat(~cldat.removed_voxels(~dat.removed_voxels)) = i;
            fprintf('Done.\n');
        end
    
    else
        try
            cldat = region2imagevec(cl(1));
        catch
            cl(1) = cluster2region(cl(1));
            cldat = region2imagevec(cl(1));
        end
        cldat = replace_empty(cldat);
        dat = cldat;
        dat.dat = zeros(size(dat.dat));
        
        for i = 1:numel(cl)
            fprintf('Working on cluster %04d among %04d clusters... ', i, numel(cl));
            cldat = region2imagevec(cl(i));
            cldat.dat(cldat.dat ~= 0) = 1;
            cldat = replace_empty(cldat);
            dat.dat(cldat.dat ~=0) = i;
            fprintf('Done.\n');
        end
    end

    dat.fullpath = fname;
    
    if dowrite
        write(dat);
    end

elseif doorigi
    fprintf('You chose to make a mask with original mask values.\n');
    try
        cldat = region2imagevec(cl(1));
    catch
        cl(1) = cluster2region(cl(1));
        cldat = region2imagevec(cl(1));
    end
    cldat = replace_empty(cldat);
    dat = cldat;
    dat.dat = zeros(size(dat.dat));

    for i = 1:numel(cl)
        fprintf('Working on cluster %04d among %04d clusters... ', i, numel(cl));
        cldat = region2imagevec(cl(i));
        cldat.dat(cldat.dat ~= 0) = 1;
        cldat = replace_empty(cldat);
        if ~isempty(cl(i).val)
            dat.dat(cldat.dat ~=0) = cl(i).val;
        elseif ~isempty(cl(i).Z)
            dat.dat(cldat.dat ~=0) = cl(i).Z';
        else
            error('cl,val and cl.Z is empty');
        end
        fprintf('Done.\n');
    end

    dat.fullpath = fname;
    
    if dowrite
        write(dat);
    end
end

end

