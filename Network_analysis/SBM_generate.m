function A = SBM_generate(M,z, varargin)

% usage: A = SBM_generate(M,z, varargin)
%
% feature: generate a simple Stochastic Block Model (SBM) (undirected, 
%          unweighted, no multi- or self-loop) 
%
% input:   M      a k x k stochastic block matrix, where Muv gives the
%                 probability that a vertex in group u is connected to 
%                 a vertex in group v
%
%          z      a n x 1 vector where z(i) gives the group index of 
%                 vertex i. z(i) can have {1,2,...,k} values. 
% 
% output:  A      adjacency matrix
%
% optional input:   
%          'sparse'   make adjacency matrix as a sparse matrix
% 
% % example:
% M = [.2 .05; .05 .2];
% z = [repmat(1,1,30) repmat(2,1,20)];
% A = SBM_generate(M,z);
%
% All calculations are based on the lecture note of Aaron Clauset's 
% Network analysis and modeling class (Fall 2014).
% see  http://tuvalu.santafe.edu/~aaronc/courses/5352/

dosparse = true; % sparse default

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            % functional commands
            case {'sparse'}
                dosparse = true;
        end
    end
end

k = unique(z);
if ~isequal(length(k), size(M,1))
    error('length(unique(z)) and size(M,1) are different. Please check.');
end

n = length(z);

A = zeros(n,n);

for i = 1:numel(k)
    for j = i:numel(k)
        
        if i == j
            A_temp = triu(rand(sum(z == k(i)), sum(z == k(j))) < M(i,j),1);
        else
            A_temp = rand(sum(z == k(i)), sum(z == k(j))) < M(i,j);
        end
        
        A(z == k(i), z == k(j)) = A_temp;
        
    end
end

if dosparse, A = sparse(A); end

A = A+A';

end