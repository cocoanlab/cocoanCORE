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
    u_comm_after_idx = unique(Ci_sorted(:,t));
    u_comm_after_idx(isnan(u_comm_after_idx)) = [];
    for u_i = 1:numel(u_comm_after_idx)
        for u_j = 1:numel(u_comm_before_idx)
            overlap_mat(u_j, u_i) = sum(Ci_sorted(:,t-1) == u_comm_before_idx(u_j) & Ci_sorted(:,t) == u_comm_after_idx(u_i));
        end
    end
    [sort_overlap_val, sort_overlap_idx] = sort(overlap_mat(:), 'descend');
    [sort_overlap_rowidx, sort_overlap_colidx] = ind2sub(size(overlap_mat), sort_overlap_idx);
    
    u_comm_sorted_idx = NaN(size(u_comm_after_idx));
    wh_count = true(size(sort_overlap_val));
    wh_count(sort_overlap_val == 0) = false;
    while true
        count_idx = find(wh_count, 1);
        if isempty(count_idx); break; end
        u_comm_sorted_idx(sort_overlap_colidx(count_idx)) = sort_overlap_rowidx(count_idx);
        wh_count(sort_overlap_rowidx == sort_overlap_rowidx(count_idx)) = false;
        wh_count(sort_overlap_colidx == sort_overlap_colidx(count_idx)) = false;
    end
    u_comm_sorted_idx(isnan(u_comm_sorted_idx)) = [1:sum(isnan(u_comm_sorted_idx))] + max(u_comm_before_idx);

    for u_i = 1:numel(u_comm_after_idx)
        Ci_sorted(Ci(:,t) == u_comm_after_idx(u_i), t) = u_comm_sorted_idx(u_i);
    end
end

end