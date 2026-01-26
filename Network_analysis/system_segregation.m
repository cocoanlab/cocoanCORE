function [S_all, S_each, S_eachtype] = system_segregation(A, Ci, Ti) 

% Calculate systematic segregation, which is defined as the degree to which
% within-community connectivity is higher than between-community connectivity.
% Modified from
% https://github.com/mychan24/system-segregation-and-graph-tools/blob/master/matlab_scripts/segregation.m.
% (See also Chan et al., 2014, PNAS)
%
% :Usage:
% ::
%     [S_all, S_each, S_type] = system_segregation(A, Ci, Ti) 
%
%
% :Input:
% ::
%   - A                  Input adjacency matrix (N X N).
%                        Z-transformed matrix is preferred.
%   - Ci                 Community index (N X 1).
%   - Ti                 Community type index (N_C X 1).
%
%
% :Output:
% ::  
%   - S_all              System segregation across all communities.
%   - S_each             System segregation for each communitiy.
%   - S_eachtype         System segregation for each communitiy type.
%
%
% :Examples:
% ::
%
%    [S_all, S_each, S_eachtype] = system_segregation(A, Ci, Ti)
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

mA_within = mean(A(wh_within));
mA_between = mean(A(wh_between));

S_all = (mA_within - mA_between) ./ mA_within;

u_comm = unique(Ci);

for comm_i = 1:numel(u_comm) % loop through communities
    
    wh_each = Ci == u_comm(comm_i); % find index for within communitiy edges
    wh_within_each = wh_within;
    wh_within_each(~wh_each, :) = false;
    wh_between_each = wh_between;
    wh_between_each(~wh_each, :) = false;
    
    mA_within_each(comm_i,1) = mean(A(wh_within_each));
    mA_between_each(comm_i,1) = mean(A(wh_between_each));

end

S_each = (mA_within_each - mA_between_each) ./ mA_within_each;

if nargin > 2

    assert(numel(u_comm)==numel(Ti));
    u_commtype = unique(Ti(Ti~=0 & ~isnan(Ti)));
    
    for commtype_i = 1:numel(u_commtype) % loop through community types
        wh_eachtype = Ti == u_commtype(commtype_i); % find index for communitiy type
        mA_within_eachtype(commtype_i,1) = mean(mA_within_each(wh_eachtype));
        mA_between_eachtype(commtype_i,1) = mean(mA_between_each(wh_eachtype));
    end
    
    S_eachtype = (mA_within_eachtype - mA_between_eachtype) ./ mA_within_eachtype;

end

end