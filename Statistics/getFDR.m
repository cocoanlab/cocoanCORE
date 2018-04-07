function [pFDR] = getFDR(p,q)
% function [pFDR] = getFDR(p,q)
%   p: a matrix of p-values across each voxel (eg: p=t.p; %in case of t-image)
%   q: q-value for FDR (eg: q=0.05)
%   V: the number of voxels tested (eg: V=90347)
%   by ChoongWan Woo

pFDR = 1;
V = length(p);
P=sort(p);
for i=1:V
    if P(i) <= i*q/V
    else
        if i~=1
            pFDR = P(i-1);
        else
            pFDR = P(1);
        end
        break
    end
end

end