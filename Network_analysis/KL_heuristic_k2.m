function [bestL, bestP, info] = KL_heuristic_k2(A, varargin)

% usage: [bestL, bestP, info] = KL_heuristic_k2(A, varargin)
%
% feature: use the Kernighan-Lin (KL) heuristic to optimize any partition 
%          score function, e.g., modularity Q or stochastic block model's
%          likelihood function. This works only for 2 partitioning problem.
%
% input:   A        adjacency matrix
%
% output:  bestL    maximum likelihood over all j, n-j paritioning.
%          bestP    the best partitioning that provides maximum likelihood
%          info.L   likelihood over t rounds. 
%
% optional input:   
%    'modularity'   This option uses modularity Q as a partitioing 
%                   score function. 
%    'fixed'        This option assigns a fixed partitioning membership 
%                   to certain vertices (e.g., known vertices). 
%                'fixed' should be followed by two vectors. One is the 
%                   vertex index/number that to be fixed. Also you need to
%                   provide the partitioning. 
%                   e.g., 'fixed', fixed, z (see the below example)
%    'dcsbm'        This option uses degree-corrected SBM likelihood as 
%                   a partitioing score function. 
%
% ** TERMINATION OPTIONS **
%
%    'run_to_n'{default}  The algorithm is repeated n times.
%    'early_stop'         The algorithm is terminated at the point where 
%                         the score (e.g., log-likelihood) doesn't increase
%    'lag'                The algorithm is terminated when the score does 
%                         not increase lagn (which should be specified) 
%                         times. e.g., KL_heuristic_k2(A, 'lag', 5)
% 
% % example
% A = [0 1 1 0 0 0; 1 0 1 0 0 0; 1 1 0 1 0 0; 0 0 1 0 1 1; ...
%      0 0 0 1 0 1; 0 0 0 1 1 0];
% fixed = [1,3];
% z = [1,2];
% [bestL, bestP] = KL_heuristic_k2(A, 'fixed', fixed, z)
%
% All calculations are based on the lecture note of Aaron Clauset's 
% Network analysis and modeling class (Fall 2014).
% see  http://tuvalu.santafe.edu/~aaronc/courses/5352/

n = size(A,1);
P0 = ones(1,n);
fixed = false(n,1);
scoref = @(A,P) SBM_likelihood(A,P); % SBM is a default
run_to_n = true; % default
do_earlystop = false;
do_lag = false;

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            % score function
            case {'modularity'}
                scoref = @(A,P) modularity_wani(A,P);
            case {'fixed'}
                fixed(varargin{i+1}) = true;       % fixed some vertices
                P0(varargin{i+1}) = varargin{i+2}; % assign groups
            case {'dcsbm'}
                scoref = @(A,P) SBM_likelihood(A,P, 'dcsbm'); 
            % termination options
            case {'run_to_n'}  
                % default; do nothing
            case {'early_stop'}  
                            % stop at the point where L does not increase
                do_earlystop = true; run_to_n = false;
            case {'lag'}    % stop if L does not increase lagn times. 
                do_lag = true; lagn = varargin{i+1}; 
                run_to_n = false;
        end
    end
end

for j = 150:((n-sum(fixed))/2)
    fprintf('Working on %d/%d\n', j, ((n-sum(fixed))/2));
    % t = 1
    P{j,1} = rand_partition(n,j,fixed,P0); % random (j,n-j) partition
    L(j,1) = scoref(A,P{j,1});             % initial log-likelihood of P
    
    % t = 2,3,...
    t = 1;
    continuing = true;
    testn = [];
    testm = -Inf;
    while continuing
        P1 = P{j,t};
        t = t + 1;
        [P{j,t}, L(j,t)] = swapping(A, P1, fixed, scoref);
        
        if do_earlystop  % stop at the point where L does not increase
            continuing = L(j,t) > L(j,t-1);
        elseif do_lag    % stop if L does not increase lagn times. 
            testm = max(testm, L(j,t));
            if L(j,t) < testm || L(j,t) == L(j,t-1)
                testn(end+1) = 1;
            else
                testn = [];
            end
            continuing = (sum(testn) < lagn) && (t < n); % upper limit is n
        elseif run_to_n 
            continuing = t < n;
        end

    end
    
end

L(L==0) = -Inf;

bestL = max(max(L));
[i,j] = find(L == bestL);
bestP = P{i,j};

% save P and L for t (the number of rounds)
% info.P = P;
info.L = L;

end

% =============== SUBFUNCTIONS ===============

function [P, P_idx] = rand_partition(n,j,fixed,P)

% CHOOSE a random (j, n-j)-partition

fn = sum(fixed);
temp_P_idx = datasample(1:(n-fn),j,'Replace',false);
P_idx = find(~fixed);
P_idx = P_idx(temp_P_idx);
P(P_idx) = 2;

end

function [P, L] = swapping(A, P0, fixed_s, scoref)

% GET possible pairs that can be swapped.

[~, min_idx] = min([sum(P0(~fixed_s)==1), sum(P0(~fixed_s)==2)]);

for i = 1:sum(P0(~fixed_s) == min_idx)

    nonfix_idx = find(~fixed_s);
    eta_i = datasample(nonfix_idx(P0(~fixed_s)==1), 1);
    eta_j = datasample(nonfix_idx(P0(~fixed_s)==2), 1);
    
    P0(eta_i)=2;
    P0(eta_j)=1;
    
    fixed_s(eta_i) = true;
    fixed_s(eta_j) = true;
    
    P_temp(i,:) = P0;
    L_temp(i) = scoref(A, P_temp(i,:));
end

[~, max_idx] = max(L_temp);
L = L_temp(max_idx);
P = P_temp(max_idx,:);

end