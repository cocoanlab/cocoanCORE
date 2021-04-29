function [A, thresh_perm, rng_status] = calc_allegiance(C, varargin)

% This function calculates allegiance matrix.
%
%
% :Usage:
% ::
%     A = calc_allegiance(C, varargin)
%
%
% :Input:
% ::
%   - C                  Community index matrix.
%                        Each column indiactes community index.
%                        If the number of column > 1, then output will
%                        be the sum of allegiance for each column.
%
%
% :Optional Input:
% ::
%   - verbose            Print current progress out.
%   - threshold          threshold using permutation of nodes.
%                        'max': use maximum allegiance for thresholding
%                        'mean': use mean allegiance for thresholding
%                        'median': use median allegiance for thresholding
%                        'prctile01', 'prctile05', 'prctile10', 'prctile20'
%                         : use top 1%, 5%, 10%, 20% for thresholding
%   - rndseed            random generator seed.
%
%
% :Output:
% ::
%   - A                  Allegiance matrix.
%                        If the number of column > 1, then output will
%                        be the sum of allegiance for each column.
%                        Data type can be variable depending on the number
%                        of columns.
%   - thresh_perm        Maximum allegiance of permuted community. This
%                        value is used for thresholding allegiance matrix.
%
%
% :Example:
% ::
%
%
% ..
%     Author and copyright information:
%
%     Copyright (C) Oct 2019  Jae-Joong Lee
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

doverbose = 0;
dothreshold = 0;
rndseed = NaN;
thresh_type = 'max';
for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            case {'verbose'}
                doverbose = 1;
            case {'threshold'}
                dothreshold = 1;
                thresh_type = varargin{i+1};
            case {'rndseed'}
                rndseed = varargin{i+1};
        end
    end
end

if isnan(rndseed)
    rng_status = rng('shuffle');
else
    rng_status = rng(rndseed);
end

[n_node, n_rep] = size(C);

if n_rep <= intmax('uint8') 
    A = zeros(n_node, n_node, 'uint8');
elseif n_rep > intmax('uint8') && n_rep <= intmax('uint16')
    A = zeros(n_node, n_node, 'uint16');
elseif n_rep > intmax('uint16') && n_rep <= intmax('uint32')
    A = zeros(n_node, n_node, 'uint32');
elseif n_rep > intmax('uint32')
    A = zeros(n_node, n_node, 'uint64');
end

if dothreshold
    A_perm = A;
    for rep_i = 1:n_rep
        if doverbose; fprintf('Working on Permutation: repetition %d ... \n', rep_i); end
        Ci = C(randperm(n_node),rep_i);
        u_c = unique(Ci);
        u_c(isnan(u_c)) = [];
        for comm_i = 1:numel(u_c)
            wh_comm = Ci == u_c(comm_i);
            A_perm(wh_comm, wh_comm) = A_perm(wh_comm, wh_comm) + 1;
        end
    end
    A_perm = A_perm(triu(true(n_node,n_node), 1));
    
    switch thresh_type
        case 'max'
            thresh_perm = max(A_perm);
        case 'mean'
            thresh_perm = mean(A_perm);
        case 'median'
            thresh_perm = median(A_perm);
        case 'prctile01'
            sort_perm = sort(A_perm, 'descend');
            wh_prctile = round(numel(sort_perm) * 0.01);
            thresh_perm = sort_perm(wh_prctile);
        case 'prctile05'
            sort_perm = sort(A_perm, 'descend');
            wh_prctile = round(numel(sort_perm) * 0.05);
            thresh_perm = sort_perm(wh_prctile);
        case 'prctile10'
            sort_perm = sort(A_perm, 'descend');
            wh_prctile = round(numel(sort_perm) * 0.10);
            thresh_perm = sort_perm(wh_prctile);
        case 'prctile20'
            sort_perm = sort(A_perm, 'descend');
            wh_prctile = round(numel(sort_perm) * 0.20);
            thresh_perm = sort_perm(wh_prctile);
        case 'prctile50'
            sort_perm = sort(A_perm, 'descend');
            wh_prctile = round(numel(sort_perm) * 0.50);
            thresh_perm = sort_perm(wh_prctile);
        case 'half'
            thresh_perm = n_rep / 2;
    end
            
    A_perm = [];
end

for rep_i = 1:n_rep
    if doverbose; fprintf('Working on Allegiance: repetition %d ... \n', rep_i); end
    Ci = C(:,rep_i);
    u_c = unique(Ci);
    u_c(isnan(u_c)) = [];
    for comm_i = 1:numel(u_c)
        wh_comm = Ci == u_c(comm_i);
        A(wh_comm, wh_comm) = A(wh_comm, wh_comm) + 1;
    end
end
A(1:n_node+1:end) = 0;
A = double(A) ./ n_rep;

if dothreshold
    thresh_perm = double(thresh_perm) ./ n_rep;
    A = A .* double(A > thresh_perm);
end

end