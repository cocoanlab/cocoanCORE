function obj = global_norm_spm_style(obj)

% This function implements the Global Normalization of the SPM manual 8.7. 
% with the 'None' option: values x 100/gs, which is actually L1norm with
% making the average L1 distance as 100

[r,c] = size(obj.dat);
obj.dat = obj.dat .* (100 .* r .* c ./ sum(sum(abs(obj.dat))));
    
end