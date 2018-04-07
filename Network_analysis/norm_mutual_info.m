function nmi = norm_mutual_info(za, zb)

% function nmi = norm_mutual_info(za, zb)
%
% feature: This function calculates the normalized mutual information,
%          which tells us how much information we learn about Community 
%          structure A (or B) if we know B (or A). If A and B are 
%          identical, we learn everything about A from B (or B from A). 
%          In this case, nmi returns 1. If they are entirely uncorrelated, 
%          we learn nothing. In this case, nmi returns 0. 
%
% input: za      membership vector for community structure A
%        zb      membership vector for community structure B
%
% output: nmi    normalized mutual information
%
% Reference:
%   *Eq.(2) from Danon et al. (2005) Comparing community structure identification
%   Also see Karrer et al. (2008) Robustness of community structure in networks

ua = unique(za);
ub = unique(zb);
N = length(za);

for i = 1:numel(ua)
    for j = 1:numel(ub)
        n(i,j) = sum(za == ua(i) & zb == ub(j));
    end
end

ni = sum(n,2); % row sum
nj = sum(n,1); % column sum

% numerator of Eq. (2) of Danon et al. (2005)
numerator = 0;
for i = 1:numel(ua)
    for j = 1:numel(ub)
        if n(i,j) ~= 0
            numerator = numerator + n(i,j) .* log((n(i,j).*N)./(ni(i).*nj(j)));
        end
    end
end
numerator = -2 .* numerator;

% denominator of Eq. (2) of Danon et al. (2005)
den1 = 0;
for i = 1:numel(ua)
    den1 = den1 + ni(i) .* log(ni(i)./N);
end

den2 = 0;
for j = 1:numel(ub)
    den2 = den2 + nj(j) .* log(nj(j)./N);
end

% Eq. (2) of Danon et al. (2005)
nmi = numerator ./ (den1 + den2);

end