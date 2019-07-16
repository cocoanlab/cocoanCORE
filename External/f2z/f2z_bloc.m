function z = f2z_bloc(f, ndf, ddf)
%f2z_bloc   Convert F variates to z variates  (blocked version)
%
%       z = f2z_bloc(f, ndf, ddf)
%
% The f2z_bloc function converts an array of F variates with ndf degrees of
% freedom in the numerator and ddf degrees of freedom in the denominator
% to an array of standard normal variates.  The calculation is arranged 
% to provide accurate answers for both very small and very large F values.
% See the article "Accurate Computation of the F-to-z and t-to-z
% Transforms for Large Arguments" by Paul Hughett for more details.

% NOTE: This code might be inefficient (and is certainly ugly) for the
% case that both ndf and ddf are scalars.

% See also f2z.m for a unblocked (and thus simpler) version.  This
% blocked version can provide about a 30% speed improvement, but only if 
% you get the blocking factor correctly matched to the computer's cache size.

% Specify the block size; must be tuned to the particular computer.
nblock = 10000;

% Verify consistency of input arguments
felem = prod(size(f));
nelem = prod(size(ndf));
delem = prod(size(ddf));
if felem > 1
    rsize  = size(f);
    ntotal = felem;
    if (nelem ~= 1) && (nelem ~= ntotal)
        error('f2z_bloc:badopt', 'Argument sizes not compatible');
    end
    if (delem ~= 1) && (delem ~= ntotal)
        error('f2z_bloc:badopt', 'Argument sizes not compatible');
    end
elseif nelem > 1
    rsize  = size(ndf);
    ntotal = nelem;
    if (delem ~= 1) && (delem ~= ntotal)
        error('f2z_bloc:badopt', 'Argument sizes not compatible');
    end
else
    rsize  = size(ddf);
    ntotal = delem;
end

% Create result array; extract block of scalar arguments
z = zeros(rsize);
if felem == 1
    f1 = f;  kf1 = 1;  kf2 = 1;
end
if nelem == 1
    ndf1 = ndf;  kn1 = 1;  kn2 = 1;
end
if delem == 1
    ddf1 = ddf;  kd1 = 1;  kd2 = 1;
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

    % Extract current input blocks
    if felem > 1
        f1 = f(n1:n2);
    end
    if nelem > 1
        ndf1 = ndf(n1:n2);
    end
    if delem > 1
        ddf1 = ddf(n1:n2);
    end

    % Allocate output block
    z1 = zeros(1,nb);

    % Compute crossover point
    c1 = finv(0.5, ndf1, ddf1);
    k1 = (f1 <= c1);
    k2 = ~k1;

    % Build selectors
    if felem > 1
        kf1 = k1;   kf2 = k2;
    end

    if nelem > 1
        kn1 = k1;   kn2 = k2;
    end

    if delem > 1
        kd1 = k1;   kd2 = k2;
    end

    % Compute accurate result for small F
    p = fcdf(f1(kf1), ndf1(kn1), ddf1(kd1));
    z1(k1) = norminv(p);

    % Compute accurate result for large F
    q = fcdf(1.0 ./ f1(kf2), ddf1(kd2), ndf1(kn2));
    z1(k2) = -norminv(q);

    % Save a block of results
    z(n1:n2) = z1;

end

