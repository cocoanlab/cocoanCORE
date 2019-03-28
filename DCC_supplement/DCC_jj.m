function DCC_mat = DCC_jj(roi_values, varargin)
% function DCC_mat = DCC_jj(roi_values, varargin)
%
% Estimate a multivariate GARCH model using the DCC estimator of Engle and Sheppard
% 
% INPUTS:
%
%     roi_values      Timeseries of ROI values. T * P matrix (Time * Parcel)
%
%
% Optional Inputs:
%
%     'whiten'
%                     Use ARMA(1, 1) model to remove autocorrelation
%
%     'simple'
%                     Estimate parameters from the whole connections.
%                     Faster than connection-wise parameter estimation.
%                     However, assumption may not be accurate.
%
%     'residualize'
%                     Followed by model matrix X, [time points x predictors], to remove
%                     prior to whitening and dynamic correlation estimation
%
%     'parellel'
%                     For parellel computing.
%
%     'dosaveload' 
%                     Only without 'simple' option.
%                     Followed by savedir, to save and load data from savedir
%                     which contains and will contain temporalily saved file.
%
%     'doverbose' 
%                     Print out current estimating process.
%                     
%
% OUTPUTS:
%
%     DCC_mat             P by P by T array of conditional correlations
%
%
% EXAMPLES :
%
%     DCC_mat = DCC_jj(roi_values, 'whiten', 'simple', 'doverbose');
%     DCC_mat = DCC_jj(roi_values, 'whiten', 'doverbose', 'dosaveload', savedir);
%
% Modified by J.J. Lee
% 2017.07.14


%% Get optional variables

dowhiten = 0;
onlywhiten = 0;
doresid = 0;
dopar = 0;
dosaveload = 0;
dosimple = 0;
doflat = 1;
doverbose = 0;

% optional inputs with default values
for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}

            case {'white', 'whiten'}
                dowhiten = 1;
                
            case {'onlywhiten'}
                dowhiten = 1;
                onlywhiten = 1;
                
            case {'simple'}
                dosimple = 1;
                
            case {'resid', 'residualize', 'X', 'model'}
                doresid = 1; 
                X = varargin{i + 1}; varargin{i + 1} = [];
                
            case {'parellel'}
                dopar = 1;
                
            case {'dosaveload'}
                dosaveload = 1;
                savedir = varargin{i + 1}; varargin{i + 1} = [];
                
            case {'noflat'}
                doflat = 0;

            case {'doverbose'}
                doverbose = 1;
                
            otherwise, warning(['Unknown input string option:' varargin{i}]);
        end
    end
end


%% Application of optional variables

[t1,k1]=size(roi_values);

if doresid
   
    % break up for efficiency - otherwise Matlab slows down with large
    % matrices...or at least it used to...
    px = pinv(X);
    pxy = px * roi_values;
    xpxy = X * pxy;
    roi_values = roi_values - xpxy; % X * pinv(X) * y;
    
end

if dowhiten
    
    Mdl = arima(1,0,1);
    
    opts = optimoptions('fmincon');
    opts.Algorithm = 'sqp';
    deltaX = 10^(-5);
    opts.TolX = deltaX;
    opts.TolCon = 10^(-6);
    opts.MaxFunEvals = 1000;
    
    indx = 0;
    fprintf('ARMA(1,1) on %d variables: %05d', k1, indx);
    
    if dopar
        parfor k = 1:k1
            EstMdl = estimate(Mdl, roi_values(:,k), 'Display', 'off', 'options', opts);
            [res, v] = infer(EstMdl,roi_values(:,k));
            roi_values(:,k) = res./sqrt(v);
        end
    else
        for k = 1:k1
            EstMdl = estimate(Mdl, roi_values(:,k), 'Display', 'off', 'options', opts);

            [res, v] = infer(EstMdl,roi_values(:,k));
            roi_values(:,k) = res./sqrt(v);

            indx = indx + 1;
            fprintf('\b\b\b\b\b%05d', indx);
        end
    end
    
    fprintf('\n');
    
end

%% Main function

if onlywhiten
    
    DCC_mat = roi_values;
    
elseif dosimple
    
    DCC_mat = DCC_CAPS2_garch(roi_values, doverbose);
    if doflat
        DCC_mat = shiftdim(DCC_mat, 2);
        DCC_mat = DCC_mat(:, triu(true(k1,k1),1))';
    end
    
elseif ~dosimple
    
    DCC_mat = zeros(k1*(k1-1)/2, t1);
    dat = zeros(t1, 2);
    tempdat = zeros(2, 2, t1);
    
    [~, ~, ~, ~, c_taketime] = sec2hms(0);
    count_i = 0;
    fprintf('total %.6d, estimated time %s : %.6d', k1*(k1-1)/2, c_taketime, count_i);
    
    
    if dosaveload && exist(savedir, 'file')
        load(savedir);
        startpoint = find(DCC_mat==0, 1);
    else
        startpoint = 1;
    end
    
    
    
    if dopar
        parfor c_i = 1:k1*(k1-1)/2
            c2_i = ceil(sqrt(2*c_i+0.25)-0.5); % n(n-1)/2 < x < n(n+1)/2.  sqrt(2x+1/4)-1/2 <= n < sqrt(2x+1/4)+1/2.
            c1_i = c_i - c2_i*(c2_i-1)/2;
            c2_i = c2_i + 1; % because of triu
            dat = [roi_values(:,c1_i) roi_values(:,c2_i)];
            tempdat = DCC_CAPS2_garch(dat, doverbose);
            DCC_mat(c_i,:) = squeeze(tempdat(1,2,:));
        end
    else
        tic;
        for c2_i = 2:k1
            for c1_i = 1:(c2_i-1)
                
                count_i = count_i + 1;
                if count_i >= startpoint
                    dat = [roi_values(:,c1_i) roi_values(:,c2_i)];
                    tempdat = DCC_CAPS2_garch(dat, doverbose);
                    DCC_mat(count_i,:) = squeeze(tempdat(1,2,:));
                    
                    if mod(count_i, 10) ~= 1
                        fprintf(repmat('\b', 1, 6));
                        fprintf('%.6d', count_i);
                    elseif mod(count_i, 10) == 1
                        c_time = toc;
                        [~, ~, ~, ~, c_taketime] = sec2hms(c_time * (k1*(k1-1)/2) / count_i);
                        fprintf(repmat('\b', 1, 32));
                        fprintf('%s : %.6d', c_taketime, count_i);
                    end
                    
                    if dosaveload && mod(count_i, 1800) == 0 % almost 1 hr
                        save(savedir, 'out');
                    end
                end
                
            end
        end
    end
    
    if ~doflat
        out_temp_3d = zeros(k1,k1,t1);
        for t_i = 1:t1
            out_temp_2d = zeros(k1,k1);
            out_temp_2d(triu(true(k1,k1),1)) = DCC_mat(:,t_i);
            out_temp_3d(:,:,t_i) = out_temp_2d + out_temp_2d';
        end
        DCC_mat = out_temp_3d;
    end
end

fprintf('    done.\n');
    


end
