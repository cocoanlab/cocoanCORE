function idx = quartile_idx(x, dotertile)

% idx = quartile_idx(x, dotertile)
%
% returns 1-4 for quartiles
% returns 1-3 for tertiles if dotertile is true

idx = zeros(size(x));

if ~dotertile
    idx(x < prctile(x,25)) = 1;
    idx(x >= prctile(x,25) & x < prctile(x,50)) = 2;
    idx(x >= prctile(x,50) & x < prctile(x,75)) = 3;
    idx(x >= prctile(x,75)) = 4;
else
    idx(x < prctile(x,33.34)) = 1;
    idx(x >= prctile(x,33.34) & x < prctile(x,66.67)) = 2;
    idx(x >= prctile(x,66.67)) = 3;
end
    

end