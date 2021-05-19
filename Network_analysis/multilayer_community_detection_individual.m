function N = multilayer_community_detection_individual(A, coupling_type, varargin)

% This function detects multi-layer community structure within individuals.
%
%
% :Usage:
% ::
%     N = multilayer_community_detection_individual(A, coupling_type, varargin)
%
%
% :Input:
% ::
%   - A                  Adjacency matrix. Data type should be cell.
%                        Rows indicate subjects.
%                        Columns indicate runs.
%                        (e.g., 19 subjects, 4 runs => 19 X 4 cell)
%   - coupling_type      Type of coupling. Either 'ord' or 'cat'.
%                        'ord': Runs are time-continuous.
%                        'cat': Runs are not time-continuous.
%
%
% :Optional Input:
%
%   - n_repeat           The number of repetition for getting consensus
%                        community. (default: 100)
%   - thresh_type        Type of thresholding for consensus community
%                        detection. See 'calc_allegiance' function.
%                        (default: 'max')
%   - gamma              Intra-layer resolution parameter. (default: 1)
%   - omega              Inter-layer coupling parameter. (default: 1)
%
% :Output:
% ::   
%   - N                  Multi-layer community structure within individuals.
%
%
% :Example:
% ::
%
%
%     Author and copyright information:
%
%     Copyright (C) May 2020  Jae-Joong Lee
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
%
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
%
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
% ..

n_repeat = 100;
thresh_type = 'max';
gamma = 1;
omega = 1;

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            case {'n_repeat'}
                n_repeat = varargin{i+1};
            case {'thresh_type'}
                thresh_type = varargin{i+1};
            case {'gamma'}
                gamma = varargin{i+1};
            case {'omega'}
                omega = varargin{i+1};
        end
    end
end

n_subj = size(A,1);
n_run = size(A,2);
n_node = size(A{1,1}, 1);

fprintf('Multilayer community detection ... \n');

for subj_i = 1:n_subj
    
    fprintf('Working on SUBJECT %.3d ... \n', subj_i);

    switch coupling_type
        case 'ord'
            [B,twomu] = multiord(A(subj_i,:),gamma,omega);
        case 'cat'
            [B,twomu] = multicat(A(subj_i,:),gamma,omega);
        otherwise
            error('Coupling type should be either ''ord'' or ''cat'' !!');
    end

    PP = @(S) postprocess_categorical_multilayer(S, n_run);

    for rep_i = 1:n_repeat
        [S, Q, n_iter] = iterated_genlouvain(B, 10000, 0, 1, 'moverandw', [], PP);
        N{subj_i}.multi_modQ{rep_i} = Q/twomu;
        N{subj_i}.multi_modmm{rep_i} = twomu;
        N{subj_i}.multi_niter{rep_i} = n_iter;
        N{subj_i}.multi_module{rep_i} = reshape(S, n_node, n_run);
    end


    for run_i = 1:n_run
        for cons_iter = 1:10

            if cons_iter == 1
                cons_comm = cat(3, N{subj_i}.multi_module{:});
                cons_comm = squeeze(cons_comm(:,run_i,:));
            end
            [allegiance_mat, perm_max(cons_iter, 1)] = calc_allegiance(cons_comm, 'threshold', thresh_type);

            wh_conn = any(allegiance_mat ~= 0);
            allegiance_mat = allegiance_mat(wh_conn, wh_conn);

            cons_comm = NaN(n_node, n_repeat);
            [B, ~] = modularity(allegiance_mat, gamma);
            allegiance_mat = [];
            for rep_i = 1:n_repeat
                [S, ~] = genlouvain(B, 10000, 0);
                cons_comm(wh_conn, rep_i) = S;
            end
            B = [];

            if all(all(diff(cons_comm(wh_conn,:), [], 2) == 0))
                N{subj_i}.multi_module_consensus(:,run_i) = cons_comm(:,1);
                N{subj_i}.multi_niter_consensus(:,run_i) = cons_iter;
                N{subj_i}.multi_permmax_consensus{run_i} = perm_max;
                perm_max = [];
                break;
            end

        end
    end
    
    N{subj_i}.multi_module_consensus_sorted = match_community_affiliation(N{subj_i}.multi_module_consensus);
    
end

end