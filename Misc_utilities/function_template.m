function [out] = function_template(input)

% Explanation...
%
% :Usage:
% ::
%
%    [out] = function_template(input)
%
% :Inputs:
%
%   **x:**
%        blahblah
%
%   **y:**
%        blah
%
% :Optional Inputs: Enter keyword followed by variable with values
%
%   **'xx':**
%        blah blah
%
%   **'yy':**
%        blah
%
% :Some advanced options:
%
%   **'xxx':**
%        xxx
%
%   **'yyy':**
%        yyy
%
% :Output:
%
%   **h:**
%        graphic handles for a bar plot
%
% :Examples:
% ::
%
%    % data
%    y = [-0.6518   -0.6934   -0.5417   -0.6496   -0.5946   -0.3839
%        1.1511    0.9090    1.1681    1.2892    0.9346    1.1383];
% 
%    col =  [0    0.1157    0.2686];
%
%    % run
%    [out] = function_template(input)
%    set(gca, 'ytickLabel', num2str(get(gca, 'ytick')'));
%    set(gcf, 'position', [1   531   399   169]);
% 
%    savename = 'example_barwani.pdf';
% 
%    try
%        pagesetup(gcf);
%        saveas(gcf, savename);
%    catch
%        pagesetup(gcf);
%        saveas(gcf, savename);   
%    end
%
% ..
% Copyright (C) 2018  Choong-Wan Woo, Cocoan Lab
%
% ..

% Programmers' notes:
%   03/01/2018: Wani modified ...
%   03/02/2018: Wani modified ...


for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            % functional commands
            case {'xx'}
                xx = varargin{i+1};
            case {'yy'}
                yy = varargin{i+1};
        end
    end
end



end
