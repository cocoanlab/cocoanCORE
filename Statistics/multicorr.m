function corrvals = multicorr(x, y)

if size(x) ~= size(y)
    error('Size of x and y does not match.');
end

[obsnum, subjnum] = size(x);
corrvals = sum(zscore(x).*zscore(y)) ./ (obsnum-1);

end