function [pFDR] = getFDR(p,q)
% function [pFDR] = getFDR(p,q)
%   p: a matrix of p-values across each voxel (eg: p=t.p; %in case of t-image)
%   q: q-value for FDR (eg: q=0.05)
%   V: the number of voxels tested (eg: V=90347)
%   by ChoongWan Woo

pFDR = 1;
p = reshape(p, numel(p), 1);
p(isnan(p)) = []; % Suhwan added (2019.09.09)
V = length(p);
P = sort(p);
for i=1:V
    if P(i) <= i*q/V
    else
        if i~=1
            pFDR = P(i-1);
        else
            pFDR = P(1);
            % print_error(P, V, q); % Suhwan added (2019.09.06)
        end
        break
    end
end

end



%% ---------------------------------------------- %%
%                 SUB-FUNCTION                     %
% ----------------------------------------------- %%
function print_error(P, V, q)
warning([10 '** The output is not a version of corrected for FDR **', 10 ... 
    10 'We failed to find corrected p-value based on your input(q).', 10 ... 
    'Instead, we returned the most lowest p-value of the matrix of p-values. ' 10 10, ...;    
    'Please be careful when using value of output: (SEE PLOT)' 10, ...
    '(You can try to change the q-value)']);
d = []; plot(P); hold on;
for ii=1:V, d(ii) = ii*q/V'; end
plot(d); hold off; legend('the sorted matrix of p-value','corrected q-line');
xlim_max = fix(V*0.01); if xlim_max < 50, xlim_max = 50; end
set(gca,'Xlim', [1 xlim_max]); % see only 5 percent of x 
end
