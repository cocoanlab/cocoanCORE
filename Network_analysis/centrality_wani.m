function cent = centrality_wani(A, varargin)

% function cent = centrality_wani(A, varargin)
%
% feature: This function calculate four different centrality measures, 
%          including degree, eigenvector, harmonics, and betweenness 
%          centrality measures. This only works for undirected/unweighted 
%          graph for now. 
%
% input:   A      adjacency matrix
% 
% options: 'all'          returns all four centrality measures
%          'degree'       returns a degree centrality
%          'eigenvector'  returns a eigenvector centrality
%          'harmonics'    returns a harmonics centrality
%          'betweenness'  returns a harmonics centrality
%          'add', cent    add centrality values into an existing cent
% 
% example: 
%       A = [0 1 1; 1 0 0; 1 0 0];
%       cent = centrality_wani(A, 'degree')  % returns degree centrality
%       cent = centrality_wani(A, 'betweenness', 'add', cent)  % add betweenness centrality into cent
%       cent = centrality_wani(A)     % returns all centrality measures
% 
% All calculations are based on the lecture note of Aaron Clauset's Network
% analysis and modeling class (Fall 2014).
% see  http://tuvalu.santafe.edu/~aaronc/courses/5352/

doall = true;
dodeg = false;
doeig = false;
dohar = false;
dobet = false;

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            % functional commands
            case {'all'}
                doall = true;
            case {'degree', 'deg'}
                dodeg = true; doall = false;
            case {'eigenvector', 'eig', 'eigen'}
                doeig = true; doall = false;
            case {'harmonics', 'har'}
                dohar = true; doall = false;
            case {'betweenness', 'bet', 'between'}
                dobet = true; doall = false;
            case {'add'}
                cent = varargin{i+1};
        end
    end
end

N = size(A,1);

%% 1. degree centrality
if dodeg || doall
    cent.degree = sum(A,2)./(sum(sum(A))/2);  % get # of degree / m
end

%% 2. eigenvector centrality
if doeig || doall
    x = ones(N,1);          % the initial condition
    for t = 1:100000
        x(:,t+1) = (A*x(:,t))./norm(A*x(:,t),1);
                            % eq. 1 of the lecture note (Ch.3) p.4 with normalization
        if sum(abs(x(:,t+1) - x(:,t))) < .00001, break, end
                            % repeat the previous line until the change of values becomes minimal
    end
    cent.eigenvector = x(:,end);
end

%% 3. harmonic centrality
if dohar || doall
    if issparse(A)
        D = all_shortest_paths(A);          % calculate the shortest paths (D) using matlab_bgl toolbox
    else
        D = all_shortest_paths(sparse(A));  % calculate the shortest paths (D) using matlab_bgl toolbox
    end
    cent.harmonics = sum(triu(1./D,1) + tril(1./D,-1),2)./(N-1);  
                                            % eq. 8 of the lecture note (Ch.3) p.10
end

%% 4. betweenness centrality
if dobet || doall
    between_cent = zeros(N,1);              % init betweenness_centrtrality, N = # of vertices
    
    if issparse(A)
        D = all_shortest_paths(A);          % calculate the shortest paths (D) using matlab_bgl toolbox
    else
        D = all_shortest_paths(sparse(A));  % calculate the shortest paths (D) using matlab_bgl toolbox
    end
    
    for j = 1:N
        for k = 1:N
            i = find(D(j,:)+D(k,:)==D(j,k));        % finding i mediating the geodesic paths from j to k
            between_cent(i) = between_cent(i) + (1/numel(i))./(N^2);  % add betweenness coefficients to vertex i
        end
    end
    cent.betweenness = between_cent;
end

return