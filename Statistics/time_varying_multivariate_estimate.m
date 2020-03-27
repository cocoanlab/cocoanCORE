function y = time_varying_multivariate_estimate(meth,data,varargin)

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
%   **meth:**
%   Window-type options, either:
%   - 'gaussian'
%   - 'tukey'
%
% :Optional Inputs: 
%
%   **stepby:**
%
%   Sometimes input data is very high resolution, and it would take too much
%   time to work element by element across the inputs.  You can enter an
%   option here to compute estimates at every n-th lag.
%
% :Examples:
% ::
%
%    y =  time_varying_multivariate_estimate('gaussian',data,20);
%    
%    % generate random data
%    data{1} = rand(2000,100); % time x voxel
%    data{2} = rand(2000,100);
%
%    y = time_varying_multivariate_estimate('gaussian',data, 20)
%
% ..
%    By Wani Woo 
%    Last updated: Mar 2020
% ..

nshift = 0;
stepby = 1;                 % step size for shift
center_local_data = false;   % mean-center local data: good for corr, bad for moving average

% set up the kernel
% -------------------------------------------
switch meth
    
    case 'gaussian'
        
         ntrials = varargin{1};
         kern = normpdf(-3:6/ntrials:3); 
         kern = kern./max(kern);            % norm to max of 1
         
         mymax = find(kern == max(kern));
         nshift = mymax - 1;  % kernel shifts by n points; adjust
         
         kern = kern';  % column

    case 'tukey'
        % Window length is the zero-influence to zero-influence time
        kern = tukeywin(varargin{1});
        nshift = round(varargin{1} ./ 2);
        
        % kludgey adjust in case we need extra element
        %if length(i - shift : i + nshift) > length(kern)
            kern = [kern; 0];
        %end
        
    otherwise error('Unknown method.')
        
end

% set up function handle
% -------------------------------------------
if length(varargin) > 1
    fhandle = varargin{2};
else
    fhandle = @(x,y) dcor(x,y);
end

if length(varargin) > 2
    stepby = varargin{3};
end

% set up data
% -------------------------------------------

if numel(data) > 2
    error('Check the size of the input. This function can take only two sets of multivariate data.');
end

for i = 1:numel(data)
    [nobs_all(i,1),ncols_all(i,1)] = size(data{i});
end

nobs = unique(nobs_all);
ncols = unique(ncols_all);

if numel(nobs)>1 || numel(ncols)>1 
    error('Check the size of the data. They should have the same numbers of rows and columns');
end

% replicate kernel for each column
kern = repmat(kern,1,ncols);

if nobs < nshift, error('Not enough observations to support kernel.'); end

y = zeros(nobs,1);


% pad data at ends to avoid edge artifacts
% -------------------------------------------
for i = 1:numel(data)
    paddat{i} = data{i}(end:-1:end-nshift,:);
    
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
        
        if center_local_data
            dati{j} = scale(dati{j},1) .* kern;   % center and multiply by kern so data taper towards mean
        else
            dati{j} = dati{j} .* kern;
        end
        
    end
    
    y(i, :) = fhandle(dati{1}, dati{2});  %r(1,2);
    
    %y(:,i) = tmpy(nshift+1:nshift+nobs);
    
end

if stepby == 1
    y = y( (nshift+1):(nobs + nshift) );

else
    indx = [(nshift+1):stepby:(nobs + nshift) (nobs + nshift)];
    y = interp1(indx, y(indx), (nshift+1):(nobs + nshift));

end

end