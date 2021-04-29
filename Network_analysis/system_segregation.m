function [S_all, S_each] = system_segregation(A, Ci) 

% Calculate systematic segregation, which is defined as the degree to which
% within-community connectivity is higher than between-community connectivity.
% Modified from
% https://github.com/mychan24/system-segregation-and-graph-tools/blob/master/matlab_scripts/segregation.m.
% (See also Chan et al., 2014, PNAS)
%
% :Usage:
% ::
%     [S_all, S_each] = system_segregation(A, Ci) 
%
%
% :Input:
% ::
%   - A                  Input adjacency matrix (N X N).
%                        Z-transformed matrix is preferred.
%
%
% :Output:
% ::  
%   - S_all              System segregation across all communities.
%   - S_each             System segregation for each communitiy.
%
%
% :Examples:
% ::
%
%    [S_all, S_each] = system_segregation(A, Ci) 
%
%
%
%     Author and copyright information:
%
%     Copyright (C) Apr 2021  Jae-Joong Lee
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

n_node = length(A);

wh_within = Ci == Ci';
wh_within(1:n_node+1:end) = false;
wh_between = Ci ~= Ci';

A_within = A(wh_within);
A_between = A(wh_between);

S_all = (mean(A_within) - mean(A_between)) / mean(A_within);

u_comm = unique(Ci);

for comm_i = 1:numel(u_comm) % loop through communities
    
    wh_each = Ci == u_comm(comm_i); % find index for within communitiy edges
    wh_within_each = wh_within;
    wh_within_each(~wh_each, :) = false;
    wh_between_each = wh_between;
    wh_between_each(~wh_each, :) = false;
    
    A_within_each = A(wh_within_each);
    A_between_each = A(wh_between_each);
    
    S_each(comm_i,1) = (mean(A_within_each) - mean(A_between_each)) / mean(A_within_each);

end

end