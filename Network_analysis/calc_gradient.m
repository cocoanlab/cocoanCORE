function [grad_embedding, grad_lambda] = calc_gradient(A, varargin)

% This function calculates connectivity gradient with an efficient algorithm,
% heavily based on Brainspace toolbox ('https://github.com/MICA-MNI/BrainSpace').
%
%
% :Usage:
% ::
%     [grad_embedding, grad_lambda] = calc_gradient(A, varargin)
%
%
% :Input:
% ::
%   - A                  Adjacency matrix in Perason's r value.
%
%
% :Optional Input:
% ::
%   - sparsity           Level of sparsity. Top (1-sparsity) remains.
%                        (default: 0.90).
%   - kernel             Type of kernel.
%                        'CS': Cosine similarity. (default)
%                        'NA': Normalized angle.
%   - algorithm          Type of algorithm.
%                        If a cell containing list of algorithms is
%                        provided, this function iteratively generates
%                        outputs for each algorithm.
%                        'DE': Diffusion Embedding. (default)
%                        'LE': Laplacien Eigenmap.
%                        'PCA': principal component analysis.
%   - n_components       Number of components.
%                        (default: 10).
%   - delA               Delete A variable in the base workspace.
%                        This is useful for saving memory, but please
%                        be very careful about using this option.
%                        If A was not saved on your base workspace (e.g.,
%                        by calling corr and r2z function), this function
%                        will not have any effect.
%
%
% :Output:
% ::
%   - grad_embedding     Gradient vectors.
%   - grad_lambda        Explainced variance of each gradient vector.
%
%
% :Example:
% ::
%   A = corr(randn(500,1000); % generate random adjacency matrix
%   [gm_embedding, gm_lambda] = calc_gradient(A, 'sparsity', 0.9, 'kernel', 'CS', 'algorithm', 'DE', 'n_components', 10, 'delA');
%
%
% ..
%     Author and copyright information:
%
%     Copyright (C) June 2021  Jae-Joong Lee
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

fprintf('\n\n Connectivity gradient analysis ... \n\n');

grad_sparsity = 0.90;
grad_kernel = 'CS';
grad_algorithm = 'DE';
n_components = 10;

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            case {'sparsity'}
                grad_sparsity = varargin{i+1};
            case {'kernel'}
                grad_kernel = varargin{i+1};
                grad_kernel = upper(grad_kernel);
            case {'algorithm'}
                grad_algorithm = varargin{i+1};
                grad_algorithm = upper(grad_algorithm);
                if ~iscell(grad_algorithm)
                    grad_algorithm = {grad_algorithm};
                end
            case {'n_components'}
                n_components = varargin{i+1};
            case {'delA'}
                A_varname = inputname(1);
                if ~isempty(A_varname)
                    evalin('base', sprintf('clear %s', A_varname));
                end
        end
    end
end

%% Basic setting and Sanity check

n_node = size(A,1);

if size(A,1) ~= size(A,2)
    error('Not rectangular form!');
end
if any(isnan(A), 1:2)
    error('NaN corr found!');
end
if max(A, [], 1:2) > 1
    error('corr > 1 found!');
end
        
%% Thresholding

fprintf('Gradient: Thresholding, sparsity %.2f ... \n', grad_sparsity);

sort_A = sort(A, 'descend');
thr_A = sort_A(ceil(n_node * (1-grad_sparsity)),:);
clear sort_A;
if min(thr_A) < 0
    fprintf('thesholded corr < 0 found. these will be set to zero. \n')
    thr_A(thr_A < 0) = 0;
end
A(bsxfun(@lt, A, thr_A)) = 0;
clear thr_A;

fprintf('\n');

%% Kernel

fprintf('Gradient: Kernel, method %s ... \n', upper(grad_kernel));

switch grad_kernel
    
    case {'CS', 'NA'}
        
        % Cos(A,B) = A dot B / |A||B|
        
        A_dot = zeros(n_node, n_node, 'double');
        A_dot = A.' * A;
        clear A;
        l2norm_A = sqrt(diag(A_dot));
        
        kernel_A = zeros(n_node, n_node, 'double');
        kernel_A = l2norm_A * l2norm_A.';
        clear l2norm_A;
        kernel_A = A_dot ./ kernel_A;
        clear A_dot;
        kernel_A(kernel_A > 1) = 1; % precision issue
        
        % Sanity check
        
        wh_neg = kernel_A < 0;
        if any(wh_neg, 1:2)
            fprintf('cosine similarity < 0 found. these will be set to zero. \n')
            kernel_A(wh_neg) = 0;
        end
        clear wh_neg;
        
        if any(all(kernel_A == 0))
            error('all-zero column of cosine similarity found!');
        end
        
        if ~issymmetric(kernel_A)
            error('asymmetry found!');
        end
        
        % Normalized angle if specified
        
        if strcmp(grad_kernel, 'NA')
            kernel_A = 1 - acos(kernel_A)./pi;
        end
        
    otherwise
        
        error('Unknown kernel method');
        
end

fprintf('\n');

%% Algorithnm

for alg_i = 1:numel(grad_algorithm)
    
    grad_algorithm_each = grad_algorithm{alg_i};

    fprintf('Gradient: Calculation, algorithm %s ... \n', grad_algorithm{alg_i});
    
    switch grad_algorithm_each
        
        case 'PCA'
            
            [U, S, ~] = svds(kernel_A - mean(kernel_A), n_components);
            embedding = U * S;
            lambda = diag(S) .^ 2;
            clear U S;
            lambda = lambda ./ sum(lambda);
            
        case 'LE'
            
            % Construct eigenmaps (solve Ly = lambda*Dy)
            d = sum(kernel_A, 2);
%             [embedding, lambda] = eigs(diag(d)-kernel_A, diag(d), n_components+1, 'smallestabs'); % only need bottom (no_dims + 1) eigenvectors
            [embedding, lambda] = eigs(bsxfun(@times, d.^-1, diag(d)-kernel_A), n_components+1, 'smallestreal'); % only need bottom (no_dims + 1) eigenvectors
            clear d;
            
            % Sort eigenvectors in ascending order
            lambda = diag(lambda);
            [lambda, ind] = sort(lambda, 'ascend');
            lambda = lambda(2:n_components+1);
            
            % Final embedding
            embedding = embedding(:,ind(2:n_components+1));
            
        case 'DE'
            
            % default parameters
            alpha = 0.5;
            
            % Parameter for later use.
            L = zeros(n_node, n_node, 'double');
            d = sum(kernel_A,2) .^ -alpha;
            L = kernel_A .* (d*d.');
            
            M = zeros(n_node, n_node, 'double');
            d2 = sum(L,2) .^ -1;
            M = bsxfun(@times, L, d2);
            clear L d d2;
            
            % Get the eigenvectors and eigenvalues
            [eigvec,eigval] = eigs(M, n_components+1, 'largestabs');
            clear M;
            eigval = diag(eigval);
            
            % Scale eigenvectors by the largest eigenvector.
            psi = bsxfun(@rdivide, eigvec, eigvec(:,1));
            
            % Automatically determines the diffusion time and scales the eigenvalues.
            scaled_eigval = eigval(2:end) ./ (1 - eigval(2:end));
            
            % Calculate embedding and bring the data towards output format.
            embedding = bsxfun(@times, psi(:,2:(n_components+1)), scaled_eigval(1:n_components).');
            lambda = scaled_eigval;
            clear psi scaled_eigval;
            
    end
    
    grad_embedding{alg_i} = embedding;
    grad_lambda{alg_i} = lambda;
    
    clear embedding lambda;
    
end

fprintf('\n');

end