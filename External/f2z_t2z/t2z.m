function z = t2z(t, df)
%t2z   Convert t variates to z variates
%
%       z = t2z(t, df)
%
% The t2z function converts an array of t variates with df degrees of 
% freedom to an array of standard normal variates.  The calculation is 
% arranged to provide accurate answers for both very small and very large
% t values.  See the article "Accurate Computation of the  F-to-z and 
% t-to-z Transforms for Large Arguments" by Paul Hughett for more details.

% See also t2z_bloc.m for a version that computes by blocks to manage cache;
% this can provide about a 30% speed improvement, but only if you get the
% blocking factor correctly matched to the computer's cache size.

% Verify consistency of input arguments
telem = prod(size(t));
delem = prod(size(df));
if telem > 1
    rsize  = size(t);
    ntotal = telem;
    if (delem ~= 1) && (delem ~= ntotal)
        error('t2z:badopt', 'Argument sizes not compatible');
    end
else
    rsize  = size(df);
    ntotal = delem;
end

% Allocate array for results
z = zeros(rsize);

% Build selectors for large and small values
c  = zeros(size(df));
k1 = (t <= c);
k2 = ~k1;

% Extract blocks for small and large t values
if telem > 1
    t1 = t(k1);
    t2 = t(k2);
else
    t1 = t;
    t2 = t;
end
if delem > 1
    df1 = df(k1);
    df2 = df(k2);
else
    df1 = df;
    df2 = df;
end

% Compute accurate result for small t values
p = tcdf(t1, df1);
z(k1) = norminv(p);

% Compute accurate result for large t values
q = tcdf(-t2, df2);
z(k2) = -norminv(q);

