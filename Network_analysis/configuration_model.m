function A = configuration_model(ki, varargin)

% function A = configuration_model(k, varargin)
% 
% feature: This function generates configuration model for a specific 
%          degree sequence, ki. 
%
% input:    ki   degree sequence
%
% optional inputs: 'sparse'
%
% output:   A    adjacency matrix
%
% All calculations are based on the lecture note of Aaron Clauset's Network
% analysis and modeling class (Fall 2014).
% see  http://tuvalu.santafe.edu/~aaronc/courses/5352/

rng 'shuffle'; %pick random seed for randomization

m2 = sum(ki);                    % 2m
N = length(ki);                  % the network size
ki_exp = zeros(1, m2);           % the stub list, 1 x 2m vector
for i = 1:numel(ki)
    if i == 1
        ki_1 = 0;
    else
        ki_1 = sum(ki(1:i-1));
    end
    ki_exp(1, ki_1+1:ki_1+ki(i)) = repmat(i,1,ki(i));
end

ki_exp = ki_exp(randperm(m2));   % random permutation the stub list
    
% re-build the adjacent matrix
ij = reshape(ki_exp,2,m2./2)'; 
A = build_adjacency_from_list(ij, varargin);

end