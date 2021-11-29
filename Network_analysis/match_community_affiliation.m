function Ci_sorted = match_community_affiliation(Ci)

% This function matches multilayer community affiliation vectors.
% In a nutshell, two contiguous affiliation vectors are compared with each
% other, and the second vector is matched to the first vector to maximize
% the the extent of overlap.
%
%
% :Usage:
% ::
%     Ci_sorted = match_community_affiliation(Ci)
%
%
% :Input:
% ::
%   - Ci                 Community affiliation vectors. (nodes X time or category)
%
%
% :Output:
% ::   
%   - Ci_sorted          Matched community affiliation vectors.
%
%
% :Example:
% ::
%
%
%     Author and copyright information:
%
%     Copyright (C) May 2021  Jae-Joong Lee
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

Ci_sorted = Ci;
n_t = size(Ci, 2);
overlap_mat = [];

for t = 2:n_t
    u_comm_before_idx = unique(Ci_sorted(:,t-1));
    u_comm_before_idx(isnan(u_comm_before_idx)) = [];
    u_comm_after_idx = unique(Ci(:,t));
    u_comm_after_idx(isnan(u_comm_after_idx)) = [];
    for u_i = 1:numel(u_comm_after_idx)
        for u_j = 1:numel(u_comm_before_idx)
            overlap_mat(u_j, u_i) = sum(Ci_sorted(:,t-1) == u_comm_before_idx(u_j) & Ci(:,t) == u_comm_after_idx(u_i));
        end
    end
    [max_overlap_val, max_overlap_idx] = max(overlap_mat);
    u_comm_sorted_idx = u_comm_before_idx(max_overlap_idx);
    
    target_duplicate = find(histcounts(max_overlap_idx) > 1);
    if ~isempty(target_duplicate)
        fprintf('\nDivergence found between time (or category) %d and %d ... \n', t-1, t);
        source_duplicate = [];
        for dup_i = 1:numel(target_duplicate)
            source_duplicate{dup_i} = find(max_overlap_idx == target_duplicate(dup_i));
            [~, max_duplicate_idx] = max(max_overlap_val(source_duplicate{dup_i}));
            source_duplicate{dup_i}(max_duplicate_idx) = [];
        end
        source_duplicate = cat(2, source_duplicate{:});
        u_comm_sorted_idx(source_duplicate) = [1:numel(source_duplicate)] + max(u_comm_before_idx);
    end
    
    for u_i = 1:numel(u_comm_after_idx)
        Ci_sorted(Ci(:,t) == u_comm_after_idx(u_i), t) = u_comm_sorted_idx(u_i);
    end
end

end