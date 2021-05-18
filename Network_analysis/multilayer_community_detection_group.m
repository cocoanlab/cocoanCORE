function N_all = multilayer_community_detection_group(N, varargin)

% This function detects consensus community structure across individuals.
%
%
% :Usage:
% ::
%     N_all = multilayer_community_detection_group(N, varargin)
%
%
% :Input:
% ::
%   - N                  Multi-layer community structure of individuals,
%                        obtained from
%                        "multilayer_community_detection_individual"
%                        function.
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
%
% :Output:
% ::   
%   - N_all              Multi-layer community structure across individuals.
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

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            case {'n_repeat'}
                n_repeat = varargin{i+1};
            case {'thresh_type'}
                thresh_type = varargin{i+1};
            case {'gamma'}
                gamma = varargin{i+1};
        end
    end
end

n_subj = numel(N);
n_run = size(N{1}.multi_module_consensus, 2);
n_node = size(N{1}.multi_module_consensus, 1);

for run_i = 1:n_run
    for cons_iter = 1:10
        
        if cons_iter == 1
            cons_comm = NaN(n_node, n_subj);
            for subj_i = 1:n_subj
                cons_comm(:,subj_i) = N{subj_i}.multi_module_consensus(:,run_i);
            end
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
            N_all.multi_module_consensus(:,run_i) = cons_comm(:,1);
            N_all.multi_niter_consensus(:,run_i) = cons_iter;
            N_all.multi_permmax_consensus{run_i} = perm_max;
            perm_max = [];
            break;
        end
        
    end
end

comm_idx = N_all.multi_module_consensus;
comm_sorted_idx = comm_idx;
overlap_mat = [];
for run_i = 2:n_run
    u_comm_before_idx = unique(comm_sorted_idx(:,run_i-1));
    u_comm_before_idx(isnan(u_comm_before_idx)) = [];
    u_comm_after_idx = unique(comm_idx(:,run_i));
    u_comm_after_idx(isnan(u_comm_after_idx)) = [];
    for u_i = 1:numel(u_comm_after_idx)
        for u_j = 1:numel(u_comm_before_idx)
            overlap_mat(u_j, u_i) = sum(comm_sorted_idx(:,run_i-1) == u_comm_before_idx(u_j) & comm_idx(:,run_i) == u_comm_after_idx(u_i));
        end
    end
    [max_overlap_val, max_overlap_idx] = max(overlap_mat);
    target_duplicate = find(histcounts(max_overlap_idx) > 1);
    if ~ismepty(target_duplicate)
        fprintf('Divergence found between run %d and %d ... \n', run_i-1, run_i);
        source_duplicate = find(max_overlap_idx == target_duplicate);
        [~, max_duplicate_idx] = max(max_overlap_val(source_duplicate));
        source_duplicate(max_duplicate_idx) = [];
    end
    u_comm_sorted_idx = u_comm_before_idx(max_overlap_idx);
    u_comm_sorted_idx(source_duplicate) = [1:numel(source_duplicate)] + max(u_comm_before_idx);
    for u_i = 1:numel(u_comm_after_idx)
        comm_sorted_idx(comm_idx(:,run_i) == u_comm_after_idx(u_i), run_i) = u_comm_sorted_idx(u_i);
    end
end
N_all.multi_module_consensus_sorted = comm_sorted_idx;

end
