function r_out = reformat_r_new(r, varargin)

% Do some reformatting or calculations on adjacency matrix or correlation 
% matrix, e.g., flatten, reconstruct, remove diagonal, etc. 
%
% :Usage:
% ::
%     r = reformat_r_new(r, varargin)
%
% ..
%     Author and copyright information:
%
%     Copyright (C) 2018  Wani Woo
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
%
% :Inputs:
%
%   **r:**
%        This can be a flattened or full correlation OR adjacency matrix.
%
% :Optional Inputs:
%
%
%   **'include_diag':**
%        This preserves diagonal values when using other options in this
%        function.
%
%   **'flatten':**
%        This flatten (i.e., vectorize) the input matrix (removing diagonal
%        values as default)
%
%   **'reconstruct':** 
%        This option reconstuct the full matrix from the input vector 
%        (with diagonal = 0)
%
%   **'remove_diag':**
%        This option make diagonals zeros
%
%   **'one_diag' or 'ones_diag':**
%        This option make diagonals ones
%
%   **'upper_triangle':**
%        This option keep upper triangle and make lower triangle and
%        diagonals into zeros.
%
%   **'lower_triangle':**
%        This option keep lower triangle and make upper triangle and
%        diagonals into zeros.
%
%   **'symmetric_avg':**
%        This option does (r + r')/2
%
%   **'symmetric_sum':**
%        This option does (r + r')
%        Good for upper triangle or lower triangle matrix
%
%   **'fisherz' or 'r2z':**
%        Fisher's r-to-z transformation
%
%   **'z2r':**
%        Reverse fisher transformation: z-to-r; r = tanh(z)
%
%
% :Outputs:
%
%
%   **r:**
%        newly calculated/reformatted r
%
% :Examples:
% ::
%
%    r = reformat_r_new(r, 'flatten')
%    r = reformat_r_new(r, 'reconstruct')
%
% :See also:
% vis_corr, vis_network
%
% ..
%    Programmers' notes:
%    Created 8/28/18 by wani woo
%
%    8/29/18 : J.J. - add 'upper_triangle' and 'lower_triangle' option
%    (used for 'vis_network' function, because the input should not be
%    duplicate matrix)
%
%    10/9/19 : J.J. - add 'include_diag' option
%                   - add 'group_mean' and 'group_sum' option
%
%    5/20/20 : J.J. - add 'threshold' option
% ..

k = 1;

do_include_diag = false;
do_flatten = false;
do_reconstruct = false;
do_remove_diag = false;
do_upper_triangle = false;
do_lower_triangle = false;
do_symmetric_avg = false;
do_symmetric_sum = false;
do_fisherz = false;
do_ones_diag = false;
do_z2r = false;
do_group = false;
do_threshold = false;

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            % functional commands
            case {'include_diag'}
                do_include_diag = true;
                k = 0;
            case {'flatten'}
                do_flatten = true;
            case {'reconstruct'}
                do_reconstruct = true;
            case {'remove_diag'}
                do_remove_diag = true;
            case {'one_diag', 'ones_diag'}
                do_ones_diag = true;
            case {'upper_triangle'}
                do_upper_triangle = true;
            case {'lower_triangle'}
                do_lower_triangle = true;
            case {'symmetric_avg'}
                do_symmetric_avg = true;
            case {'symmetric_sum'}
                do_symmetric_sum = true;
            case {'fisherz', 'r2z'}
                do_fisherz = true;
            case {'z2r'}
                do_z2r = true;
            case {'group_mean'}
                do_group = true;
                group_fun = @(x) mean(x(:));
                group_idx = varargin{i+1};
            case {'group_sum'}
                do_group = true;
                group_fun = @(x) sum(x(:));
                group_idx = varargin{i+1};
            case {'threshold'}
                do_threshold = true;
                thresh_lvl = varargin{i+1};
        end
    end
end

%%

% get n
if size(r,1) == 1 || size(r,2) == 1 
    n = 1/2+sqrt(max(size(r,1), size(r,2))*2+1/4);
    if do_include_diag
        n = n - 1;
    end
else
    if size(r,1) ~= size(r,2)
        warning('r should be symmetric. Check your input.');
    end
    n = size(r,1);
end


if do_flatten, r_out = r(triu(true(n,n), k)); end

if do_reconstruct
    r_out = zeros(n,n);
    r_out(triu(true(n,n), k)) = r;
    r_out = r_out + r_out';
    if do_include_diag
        r_out(1:n+1:end) = r_out(1:n+1:end) ./ 2;
    end
end
    
if do_remove_diag, r_out = r; r_out(1:n+1:end) = 0; end
    
if do_ones_diag, r_out = r; r_out(1:n+1:end) = 1; end

if do_upper_triangle, r_out = r .* triu(true(n,n), k); end
    
if do_lower_triangle, r_out = r .* tril(true(n,n), -k); end

if do_symmetric_avg, r_out = (r + r')./2; end

if do_symmetric_sum, r_out = (r + r'); end
    
if do_fisherz, r_out = .5 * log( (1+r) ./ (1-r) ); end

if do_z2r, r_out = tanh(r); end

if do_group
    u_group_idx = unique(group_idx);
    n_group = numel(u_group_idx);
    r_out = zeros(n_group,n_group);
    for i = 1:numel(u_group_idx)
        for j = i:numel(u_group_idx)
            r_out(i,j) = group_fun(r(group_idx == u_group_idx(i), group_idx == u_group_idx(j)));
        end
    end
    r_out = r_out + r_out';
    r_out(1:n_group+1:end) = r_out(1:n_group+1:end) ./ 2;
end

if do_threshold
    if size(r,1) > 1 && size(r,2) > 1 % square matrix
        thresh_val = prctile(reformat_r_new(r, 'flatten'), 100 * (1 - thresh_lvl));
    else
        thresh_val = prctile(r, 100 * (1 - thresh_lvl));
    end
    r_out = r .* double(r > thresh_val);
end