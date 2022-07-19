function [varargout] = modularity(A, g)
% [Q] = modularity(A, group) computes the modularity Q of the network
% defined by the adjecency matrix A and partitioned as indicated in vector
% g.
% A is a square matrix of size N. g is a vector of length N; each entry
% g(i) of g corresponds to the group to which node i belongs; e.g., if g =
% [1 2 1], then nodes 1 and 3 belong to group 1, while node 2 belongs to
% group 2.
% If A is not symmetric, the modularity for directed networks is computed.
% Also, if A is a weighted adjacency matrix, the modularity for weighted
% networks is computed. The unweighted modularity is computed when all
% non-zero values in A are equal to 1.
% [Q, Qv] = modularity(A, group) computes the modularity Q and the vector
% Qv of the contributions of each community to the total modularity, so
% that sum(Qv) = Q
% -------------------------------------------------------------------------
% For references, see:
% https://doi.org/10.1103/PhysRevE.102.043109
% https://doi.org/10.1103/PhysRevE.70.056131
% https://doi.org/10.1103/PhysRevLett.100.118703
% https://doi.org/10.1088/1742-5468/2009/03/P03024
% -------------------------------------------------------------------------
% Code by: Davide Perrone [davide.perrone@polito.it]
nGroups = numel(unique(g));
nNodes = length(g);
nargoutchk(0,2)
e = zeros(nGroups); 
for i = 1:nNodes
    for j = 1:nNodes
        e(g(i), g(j)) = e(g(i), g(j)) + A(i, j);
    end
end
nLinks = sum(A(:));
a_out = sum(e, 2);
a_in = (sum(e, 1))';
a = a_in.*a_out/nLinks^2;
Q = trace(e)/nLinks - sum(a);
Qv = diag(e)/nLinks - a;
if nargout <= 1
    varargout{1} = Q;
elseif nargout == 2
    varargout{1} = Q;
    varargout{2} = Qv;
end
    