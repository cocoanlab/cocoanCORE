function [V_b, D_b, U_b, rank_b] = bootSVD(Y, P_b_idx, V, D, U, do_flip)

% [V_b, D_b, U_b] = bootSVD(Y, P_b_idx, V, D, U)
%
% Y_b = Y * P_b
% = (V * D * U.') * P_b      //  [V, D, U] = svd(Y - mean(Y,2), 'econ');
% = V * (D * U.' * P_b)
% = V * DUtPb                //  DUtPb = D * U.' * P_b;
% = V * (A_b * S_b * R_b.')  //  [A_b, S_b, R_b] = svd(DUtPb - mean(DUtPb,2), 'econ');
% = (V * A_b) * S_b * R_b.'
% = V_b * D_b * U_b.'        //  V_b = V * A_b; D_b = S_b; U_b = R_b;
%
% Adapted from
% "Fast, Exact Bootstrap Principal Component Analysis for p > 1 Million",
% A Fisher, B Caffo, B Schwartz, V Zipunnikov,
% Journal of the American Statistical Association 111 (514), 846-860
%
% Jae-Joong Lee, 2023. 10. 26
%
%

if nargin < 5
    [V, D, U] = svd(Y - mean(Y,2), 'econ');
end

if nargin < 6
    do_flip = false;
end

DUt = D * U.';

DUtPb = DUt(:,P_b_idx);
% assert(isequal(DUtPb, D * U.' * double(repmat(1:n,n,1).' == P_b_idx(:).')))

[A_b, S_b, R_b] = svd(DUtPb - mean(DUtPb,2), 'econ');

V_b = V * A_b;
D_b = S_b;
U_b = R_b;

rank_b = sum(diag(D_b) > max(size(Y)) * eps(norm(diag(D_b),inf)));

if do_flip
    % sign of V and Vb should ideally be same; V.' * V_b = V.' * V * A_b = A_b;
    % flip the vectors with negative A_b!
    wh_flip = sign(diag(A_b)).';
    V_b = V_b .* wh_flip;
    U_b = U_b .* wh_flip;
end

% [V_b, D_b, U_b] = svd(Y(:,P_b_idx) - mean(Y(:,P_b_idx), 2), 'econ');
% VA_b = V * A_b; VA_b = VA_b .* sign(diag(VA_b.' * V_b)).';
% histogram(V_b(:,1:rank_b) - VA_b(:,1:rank_b));

end