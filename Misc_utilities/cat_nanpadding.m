function out = cat_nanpadding(varargin)

% Concatenate vectors with different numbers of rows. It will be padding
% the extra rows with NaNs. 
%
% :Usage:
% ::
%
%    out = cat_nanpadding(x1, x2, x3, ...)
%
%
% :Inputs:
%
%   **x1, x2,...:**
%           x's should be a column vector.
%

out = nan(max(cellfun(@numel, varargin)), numel(varargin));

for i = 1:numel(varargin)
    
    out(1:numel(varargin{i}),i) = varargin{i};
    
end


end