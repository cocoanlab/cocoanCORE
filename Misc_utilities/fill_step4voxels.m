% When you have input data from 'step4voxels' searchlight method which has 
% 3297 rows, fill input data into an empty fmri_data(dat_template). 
%
% out_fmridata = fill_step4voxels(input,varargin)
%
% example :: fill_step4voxels(input, 'top', 30);
%
% the input should have 3297 rows.
%
% Copyright Byeol Kim & Wani Woo, 2019
%
function [out_fmridata] = fill_step4voxels(input, varargin)
sorting = false;
do_display = true;
use_sphere = false;     % default is cube
do_values = false;

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            case {'sort','sorting','top'}
                sorting = true;
                top_n = varargin{i+1};
            case {'no_display'}
                do_display = false;
            case {'sphere'}
                use_sphere = true;
            case {'values'}
                do_values = true;
        end
    end
end

if ~do_values
    input = zeros(3297,1);
    input(input) = 1;
elseif size(input,1) ~= 3297
    error('The row of input is not 3297. Check it out.')
end

load(which('gray_matter_step4voxels.mat'))

if sorting
    % only input has a column
    if size(input,2) == 1
        out_fmridata = dat_template;                    % 211339
        [~,idx] = sort(input, 'descend');               % idx = 3297x1
        input_sort = input(idx(1:top_n,:));
        if use_sphere
            idx_sort = fillsphere_idx_211339(:,idx(1:top_n,:));
        else
            idx_sort = fillcube_idx_211339(:,idx(1:top_n,:));
        end
        out_fmridata.dat = idx_sort * input_sort;
    else
        warning('You used sorting option but the column number of input is not 1.')
    end
    
else
    out_fmridata = dat_template;                  % 211339
    if use_sphere
        out_fmridata.dat = fillsphere_idx_211339 * input;
    else
        out_fmridata.dat = fillcube_idx_211339 * input;
    end
    
end

if do_display && size(out_fmridata.dat,2) == 1
    orthviews(out_fmridata);
end

end