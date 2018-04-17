function A = build_adjacency_from_list(ij, varargin)

% function A = build_adjacency_from_list(ij, [optional 'sparse' OR 'weights', w])
%
% feature: This function builds an adjacency matrix from the adjacency
%          list.
%
% input:   ij   - a list of the row and column numbers that contain edges
%
% optional input:   'sparse'   - make A into sparse matrix
%                   'weights'  - use w to make adjacency matrix
%
% output:  A    - adjacency matrix
%
% All calculations are based on the lecture note of Aaron Clauset's Network
% analysis and modeling class (Fall 2014).
% see  http://tuvalu.santafe.edu/~aaronc/courses/5352/

dosparse = false;
doweight = false;

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            % functional commands
            case {'sparse'}
                dosparse = true;
            case {'weights'}
                doweight = true;
                w = varargin{i+1};
                if ~isequal(length(ij), length(w))
                	error('The length of ij should be same with length of w.');
                end
        end
    end
end

A = zeros(max(max(ij)),max(max(ij)));

for i = 1:length(ij)
    if doweight
        A(ij(i,1), ij(i,2)) = w(i);
    else
        A(ij(i,1), ij(i,2)) = 1;
    end
end

if dosparse
    A = sparse(A);
end

if ~issymmetric(A)
    A = A + A';
end

end