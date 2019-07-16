function z = t2z_bloc(t, df)
%t2z_bloc   Convert t variates to z variates
%
%       z = t2z_bloc(t, df)
%
% The t2z_bloc function converts an array of t variates with df degrees of 
% freedom to an array of standard normal variates.  The calculation is 
% arranged to provide accurate answers for both very small and very large
% t values.  See the article "Accurate Computation of the F-to-z and
% t-to-z Transforms for Large Arguments" by Paul Hughett for more details.

% See also t2z.m for a unblocked (and thus simpler) version.  This
% blocked version can provide about a 30% speed improvement, but only if 
% you get the blocking factor correctly matched to the computer's cache size.

% Specify the block size; must be tuned to the particular computer.
nblock = 10000;

% Verify consistency of input arguments
telem = prod(size(t));
delem = prod(size(df));
if telem > 1
    rsize  = size(t);
    ntotal = telem;
    if (delem ~= 1) && (delem ~= ntotal)
        error('t2z_bloc:badopt', 'Argument sizes not compatible');
    end
else
    rsize  = size(df);
    ntotal = delem;
end

% Create result array; extract block of scalar arguments
z = zeros(rsize);
if telem == 1
    t1 = t;  kt1 = 1;  kt2 = 1;
end
if delem == 1
    df1 = df;  kd1 = 1;  kd2 = 1;
end

% Process in blocks to minimize cache reloads
for n1 = 1:nblock:ntotal

    % Choose limits of the current block
    n2 = n1 + nblock - 1;
    nb = nblock;
    if n2 > ntotal
        n2 = ntotal;
        nb = n2 - n1 + 1;
    end
    [n1 n2]    % DEBUG

    % Extract current input blocks
    if telem > 1
        t1 = t(n1:n2);
    end
    if delem > 1
        df1 = df(n1:n2);
    end

    % Allocate output block; build selectors
    if length(t1) > 1
        z1 = zeros(size(t1));
        c1 = zeros(size(t1));
    else
        z1 = zeros(size(df1));
        c1 = zeros(size(df1));
    end
    k1 = (t1 <= c1);
    k2 = ~k1;
    if telem > 1
        kt1 = k1;   kt2 = k2;
    end
    if delem > 1
        kd1 = k1;   kd2 = k2;
    end

    % Compute accurate result for small F
    p = tcdf(t1(kt1), df1(kd1));
    z1(k1) = norminv(p);

    % Compute accurate result for large F
    q = tcdf(-t1(kt2), df1(kd2));
    z1(k2) = -norminv(q);

    % Save a block of results
    z(n1:n2) = z1;

end

