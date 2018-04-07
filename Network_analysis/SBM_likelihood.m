function L = SBM_likelihood(A, z, varargin)

% usage: L = SBM_likelihood(A, z)
%
% feature: log-likelihood estimation for simple Stochastic Block Model (SBM)
%          (undirected, unweighted) based on Maximum Likelihood estimation. 
%
% input:   A      adjacency matrix (simple graph)
%          z      group membership
%
% optional_inputs:  'log'{default} returns log-likelihood
%                   'nolog'        returns likelihood
%                   'dcsbm'        returns log-likelihood for 
%                                  degree-corrected SBM
% 
% output:  L      likelihood for z
% 
% example: test against an example on p.10 of the Lecture note 6:
%
%  A = [0 1 1 0 0 0; 1 0 1 0 0 0; 1 1 0 1 0 0; 0 0 1 0 1 1; ...
%      0 0 0 1 0 1; 0 0 0 1 1 0];
%  z = [1 1 1 2 2 2];
%  L = SBM_likelihood(A, z); % L = 0.0433
%
%  z = [1 1 1 1 2 2];
%  L = SBM_likelihood(A, z); % L = 2.4414e-04
%
%  % Additionally this function has been tested against the Zachary karate 
%  % club example on p.12 of Lecture Note 6. 
%
% All calculations are based on the lecture note of Aaron Clauset's 
% Network analysis and modeling class (Fall 2014).
% see  http://tuvalu.santafe.edu/~aaronc/courses/5352/

dolog = true;
do_dcsbm = false;
for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            case {'log'} 
                % do nothing
            case {'nolog'} 
                dolog = false; % log-likelihood 
            case {'dcsbm'} 
                do_dcsbm = true; % degree-corrected SBM
        end
    end
end

z_uniq = unique(z);

[u,v] = find(tril(ones(length(z_uniq),length(z_uniq))));

if ~do_dcsbm
    
    % Standard SBM
    L = 1;
    for i = 1:length(u)
        
        % In this calculation, Nuv is the number of possible edges between
        % group u and v, and Nu (or Nv) is the number of vertices with labe
        % u (or v). Euv is the number of observed edges bewteen group u and
        % v (or the number of observed edges within a group u when u = v).
        
        if u(i) == v(i) % first, when u and v are in the same group
            Nu = sum(z == z_uniq(u(i)));
            if Nu == 1  % if Nu=1, then Nuu=1, Euu=0; (no self-loop)
                Nuv = 1; Euv = 0;
            else
                Nuv = nchoosek(Nu,2);
                Euv = sum(sum(A(z== z_uniq(u(i)),z==z_uniq(u(i)))))./2;
            end
        else
            Nu = sum(z == z_uniq(u(i)));
            Nv = sum(z == z_uniq(v(i)));
            
            Nuv = Nu .* Nv;
            Euv = sum(sum(A(z==z_uniq(u(i)),z==z_uniq(v(i)))));
        end
        
        L = L .* ((Euv/Nuv).^Euv.*(1-Euv/Nuv).^(Nuv-Euv)); % Eq. (3)
        
    end
    if dolog, L = log(L); end

else
    
    % Degree Corrected SBM; ref) p.11 of the Lecture Note 6, or Karrer & 
    % Newman (2011) Stochastic blockmodels and community structure in 
    % networks

    L = 0; % log L
    m2 = sum(sum(A));
    ki = sum(A);
    for i = 1:length(u)
        ku = sum(ki(z == z_uniq(u(i))));
        kv = sum(ki(z == z_uniq(v(i))));
        Euv = sum(sum(A(z == z_uniq(u(i)),z==z_uniq(v(i)))));
        if Euv == 0
            L = L + 0;
        else
            L = L + (Euv./m2) .* log((Euv.*m2)./(ku.*kv)); 
                                            % Eq. (11) of Lecture Note 6
        end
    end
    if ~dolog, L = exp(L); end
    
end

end