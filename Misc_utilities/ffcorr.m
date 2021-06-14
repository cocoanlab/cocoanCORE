function y = ffcorr(x)
% Fast & Flattened correlation
                
[n,p] = size(x);
zx = (x - mean(x)) ./ (std(x) * (n-1).^0.5); % z-scoring divided by (n-1).^0.5; do not use 'zscore.m' that ignores zero array... sad.
y = subsref(zx.' * zx, struct('type', {'()'}, 'subs', {{(1:p) > (1:p).'}}));
wh_bad = abs(y)>1;
y(wh_bad) = sign(y(wh_bad));

end