function z = f2z(f, ndf, ddf)
%f2z   Convert F variates to z variates (selector algorithm)
%
%       z = f2z(f, ndf, ddf)
%
% The f2z function converts an array of F variates with ndf degrees of
% freedom in the numerator and ddf degrees of freedom in the denominator
% to an array of standard normal variates.  The calculation is arranged 
% to provide accurate answers for both very small and very large F values.
% See the article "Accurate Computation of the F-to-z and t-to-z
% Transforms for Large Arguments" by Paul Hughett for more details.

% See also f2z_bloc.m for a version that computes by blocks to manage cache;
% this can provide about a 30% speed improvement, but only if you get the
% blocking factor correctly matched to the computer's cache size.

% Verify consistency of input arguments
felem = prod(size(f));
nelem = prod(size(ndf));
delem = prod(size(ddf));
if felem > 1
    rsize  = size(f);
    ntotal = felem;
    if (nelem ~= 1) && (nelem ~= ntotal)
        error('f2z:badopt', 'Argument sizes not compatible');
    end
    if (delem ~= 1) && (delem ~= ntotal)
        error('f2z:badopt', 'Argument sizes not compatible');
    end
elseif nelem > 1
    rsize  = size(ndf);
    ntotal = nelem;
    if (delem ~= 1) && (delem ~= ntotal)
        error('f2z:badopt', 'Argument sizes not compatible');
    end
else
    rsize  = size(ddf);
    ntotal = delem;
end

% Allocate array for results
z = zeros(rsize);

% Compute crossover point; build selectors
% NOTE: Always setting the crossover to 1.0 makes this function about
% 4x faster if ndf or ddf is a vector, but this costs accuracy if ndf
% and ddf are very different.
if 1
    c  = finv(0.5, ndf, ddf);
else
    c = 1.0;
end
k1 = (f <= c);
k2 = ~k1;

% Extract blocks for small and large F values
if felem > 1
    f1 = f(k1);
    f2 = f(k2);
else
    f1 = f;
    f2 = f;
end
if nelem > 1
    ndf1 = ndf(k1);
    ndf2 = ndf(k2);
else
    ndf1 = ndf;
    ndf2 = ndf;
end
if delem > 1
    ddf1 = ddf(k1);
    ddf2 = ddf(k2);
else
    ddf1 = ddf;
    ddf2 = ddf;
end

% Compute for small F values
p = fcdf(f1, ndf1, ddf1);
z(k1) = norminv(p);

% Compute for large F values
q = fcdf(1.0 ./ f2, ddf2, ndf2);
z(k2) = -norminv(q);

