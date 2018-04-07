function q = modularity_wani(A,z,varargin)

% function q = modularity_wani(A,z, [optional 'scalar' or 'degree'])
%
% feature: This function calculates the modularity (same with 
%          assortativity), q. It works for categorical and continuous
%          (needs optional input, 'scalar') attributes. 
% 
% input:   z   attributes or category membership
%          A   adjacency matrix
%
% optional input:
%      'scalar'   calculate the assortativity coefficient (which is a
%                 network-based generalization of the Pearson correlation
%                 coefficient). 
%      'degree'   calculate the degree assortativity coefficient (a special
%                 case of scalar assortativity coefficient. In this case,
%                 you don't need z. eg) q = modularity_wani(A,[],'degree')
%          
% output:  q   modularity or assortativity coefficient
% 
% All calculations are based on the lecture note of Aaron Clauset's Network
% analysis and modeling class (Fall 2014).
% see  http://tuvalu.santafe.edu/~aaronc/courses/5352/

doscalar = false;
dodegree = false;
docateg = true;

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            % functional commands
            case {'scalar'}
                doscalar = true;
                docateg = false;
            case {'degree'}
                dodegree = true;
                docateg = false;
        end
    end
end

% examine the data
if ~dodegree
    if length(z) ~= size(A,1)
        error('The length of z should be same with the network size. Check the data.');
    end
end

m2 = sum(sum(A)); % 2m

% 1. categorical attributes
if docateg

    u = unique(z);
    q = 0;
    
    for i = 1:numel(u)
        % Eq. (4) of Lecture note 5
        q = q + sum(sum(A(z==u(i),z==u(i))))./m2 - (sum(sum(A(z==u(i),:)))./m2)^2;
    end
    
% 2. scalar attributes
elseif doscalar
    
    k = sum(A); 
    kikj = k'*k; 
    
    if size(z,2)<size(z,1), z = z'; end
    zizj = z'*z; 
    
    % Eq. (6) of Lecture note 5
    q = sum(sum((A - kikj ./ m2) .* zizj)) ./ ...
        sum(sum(((repmat(k,length(k),1) .* eye(length(k))) - (kikj ./ m2)) .* zizj)); 

% 3. degree assortativity coefficient
elseif dodegree
    
    k = sum(A); 
    kikj = k'*k; 
    q = sum(sum((A - kikj ./ m2) .* kikj)) ./ ...
        sum(sum(((repmat(k,length(k),1) .* eye(length(k))) - (kikj ./ m2)) .* kikj)); 

end

end