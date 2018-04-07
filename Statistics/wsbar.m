function bar = wsbar(data)

% bar = wsbar(data)
% 
% A function for the calculation of within-subject error bar by Matt Jones.
% This function has been tested against barplot_get_within_ste 
% with data from Loftus & Masson, Table 2
% 
% each row should be a different subject, and columns should be repeated
% measures. 
%

resid = data - repmat(nanmean(data),size(data,1),1) - ...
    repmat(nanmean(data,2),1,size(data,2)) + mean(nanmean(data));

% replacing nans with the subject-wise mean resid values
interp_resid = isnan(resid).*repmat(nanmean(resid,2), 1, size(resid,2));
resid(isnan(resid))=0;
resid = resid+interp_resid;

ms = sum(sum(resid.^2))/prod(size(data)-1);

bar = sqrt(ms/size(data,1));

end

