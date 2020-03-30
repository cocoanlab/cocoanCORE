function y = time_varying_multivariate_estimate(data, window_width, varargin)

% Performs U-centered (bias corrected) distance correlation (default) between 
% X and Y on symmetric moving average of data
%
% This is a modified version of time_varying_estimate.m, which was in 
% https://github.com/canlab/CanlabCore/tree/master/CanlabCore/Statistics_tools
%
% :Usage:
% ::
%
%     y = time_varying_multivariate_estimate(meth, data, [window width], [function handle])
%
%
% :Inputs:
%
%   **data:** 
%   - each cell should have multivariate data 
%   - row (#common dimension, e.g., observation, time, conditions) 
%       x column (#features, e.g., #voxels)
%   - e.g., data{1} has 2000 (# time) x 200 (# voxels)
%           data{2} has 2000 (# time) x 100 (# voxels)
%
% :Optional Inputs: 
%
%   **method:** Window-type options, either:
%   - 'gaussian' 
%   - 'tukey'
% 
%     * Note: These kernels do not work with dCor method because it modified
%            the euclidean distances across conditions. 
%
%   - 'none' (default)
%     * This is the default because the dCor method is default
%
%   **stepby:**
%
%   - Sometimes input data is very high resolution, and it would take too much
%     time to work element by element across the inputs.  You can enter an
%     option here to compute estimates at every n-th lag.
%
%
% :Examples:
% ::
%
%    y =  time_varying_multivariate_estimate(data,20);
%    
%    % generate random data
%    data{1} = rand(2000,200); % time x voxel
%    data{2} = rand(2000,100);
%
%    y = time_varying_multivariate_estimate(data, 20)
%
% ..
%    By Wani Woo 
%    Last updated: Mar 2020
% ..

stepby = 1;                 % step size for shift
meth = 'none';              % default method
fhandle = @(x,y) dcor(x,y); % default fhandle: dCor
% center_local_data = false;   % mean-center local data: good for corr, bad for moving average

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            % functional commands
            case {'meth', 'method', 'methods'}
                meth = varargin{i+1};
            case {'fhandle'}
                fhandle = varargin{i+1};
            case {'stepby'}
                stepby = varargin{i+1};
        end
    end
end


% set up the kernel
% -------------------------------------------
switch meth
    
    case 'gaussian'
        
         ntrials = window_width;
         kern = normpdf(-3:6/ntrials:3); 
         kern = kern./max(kern);            % norm to max of 1
         
         mymax = find(kern == max(kern));
         nshift = mymax - 1;  % kernel shifts by n points; adjust
         
         kern = kern';  % column

    case 'tukey'
        % Window length is the zero-influence to zero-influence time
        kern = tukeywin(window_width);
        nshift = round(window_width ./ 2);
        
        kern = [kern; 0];
        
    case 'none'
        ntrials = window_width;
        nshift = round(window_width ./ 2);
        kern = ones(ntrials,1);
        kern = [kern; 0];
        
    otherwise
        
        error('Unknown method.')
        
end


% set up data
% -------------------------------------------

if numel(data) > 2
    error('Check the size of the input. Currently, this function can take only two sets of multivariate data.');
end

for i = 1:numel(data)
    [nobs_all(i,1),ncols_all(i,1)] = size(data{i});
end

nobs = unique(nobs_all);
% ncols = unique(ncols_all);
% 
if numel(nobs)>1
    error('Check the size of the data. They should have the same numbers of observations');
end

% replicate kernel for each column
% kern = repmat(kern,1,ncols);

if nobs < nshift, error('Not enough observations to support kernel.'); end

y = zeros(nobs,1);


% pad data at ends to avoid edge artifacts
% -------------------------------------------
for i = 1:numel(data)
    paddat{i} = data{i}(end:-1:end-nshift,:); % this will create biased 
                                              % distance for the padded data (for dCor)
    data{i} = [data{i}; paddat{i}];
    paddat{i} = data{i}(nshift:-1:1,:);
    data{i} = [paddat{i}; data{i}];
end


% execute
% -------------------------------------------
for i = [(nshift+1):stepby:(nobs + nshift) (nobs + nshift)]
    % start at nshift to avoid ends; data is padded
    
    % set up windowed data
    
    for j = 1:numel(data)
        dati{j} = data{j}(i - nshift : i + nshift,:);
        
%         if center_local_data
%             dati{j} = scale(dati{j},1) .* repmat(kern,1,ncols_all(j));   % center and multiply by kern so data taper towards mean
%         else
        dati{j} = dati{j} .* repmat(kern,1,ncols_all(j));
%         end
        
    end
    
    y(i, :) = fhandle(dati{1}, dati{2});  
    
end

if stepby == 1
    y = y( (nshift+1):(nobs + nshift) );
    
else
    indx = [(nshift+1):stepby:(nobs + nshift) (nobs + nshift)];
    y = interp1(indx, y(indx), (nshift+1):(nobs + nshift));

end

y([1:nshift (end-nshift+1):end]) = NaN; % make the first and last values 
                                        % within the window width NaNs to 
                                        % take into account the bias

end