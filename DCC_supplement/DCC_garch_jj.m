function Ct = DCC_garch_jj(dat, doverbose)

% function Ct = DCC_garch_jj(dat, doverbose)
%
% Modified version of dcc_mvgarch.m.
%
% 
% INPUTS:
%
%     dat             Timeseries of ROI values. T * P matrix (Time * Parcel)
%
%
% Optional Inputs:
%
%     'doverbose' 
%                     Print out current estimating process.
%
%
% OUTPUTS:
%
%     Ct              P*(P-1)/2 (=Pchoose2) by T matrix of dynamic conditional correlations
%
%
% Difference from the previous code (dcc_mvgarch.m):
%
%     - Fix dccP, dccQ, archP, garchQ as 1
%     - DCC values are the only output (output variables such as loglikelihood, Qt, stdresid, likelihoods, stderrors, A, B, and jointscores are ignored)
%     - Can deal with NaN connectivity values (this happens when some of the ROI values were zero across all timepoints).
%     - If there is a error (mostly, failure to convergence), error message is printed out.
%     - Simplify code (see below):
%       
%
%         *** Original code (dcc_mvgarch.m) ***
%
%         ...
%         %We now have Ht and the likelihood
%         [loglikelihood, Rt, likelihoods, Qt]=dcc_mvgarch_full_likelihood(parameters, data, archP,garchQ,dccP,dccQ);
%         likelihoods=-likelihoods;
%         loglikelihood=-loglikelihood;
%         Ht=zeros(k,k,t);
%         stdresid=zeros(t,k);
%         Hstd=H.^(0.5);
%         for i=1:t
%             Ht(:,:,i)=diag(Hstd(i,:))*Rt(:,:,i)*diag(Hstd(i,:));
%             stdresid(i,:)=data(i,:)*Ht(:,:,i)^(-0.5);
%         end
%         ...
%
%         and the Ht output was additionally processed in DCCsimple.m as follows:
%
%         ...
%         % Compute dynamic correlation matrix Ct
%         Ct = Ht;
%         for i=1:T,
%             Ct(:,:,i) = Ht(:,:,i)./sqrt(diag(Ht(:,:,i))*diag(Ht(:,:,i))');
%         end
%         ...
%
%         However, this process makes 'Ct' mathematically same as 'Rt'.
%         Since we are only interested in 'Ct', we remove the codes below the likelihood estimation process.
%
%
%         *** New code (DCC_garch_jj.m) ***
%
%         ...
%         [~, Ct, ~, ~] = dcc_mvgarch_full_likelihood(parameters, dat, archP, garchQ, dccP, dccQ);
%         ...
%
%         This can avoid redundant computation process and allow fast speed.
%         You can also compare this function with the original version.
%         
%         X = randn(100,20);
%         DCC_orig = DCCsimple(X); % uses 'dcc_mvgarch.m'
%         DCC_new = DCC_jj(X, 'simple', 'noflat'); % uses 'DCC_garch_jj.m'
% 
%         % correlation coefficient
%         corr(DCC_orig(:), DCC_new(:))
%
%         % histogram of differences
%         histogram(DCC_orig(:) - DCC_new(:))
%
%
%
% Created by J.J. Lee (modified from dcc_mvgarch.m)
% 2017.07.14


%% BASIC setting : Find zero timeseries to replace it with NaN

[t2, k2_all] = size(dat);

miscol = [];
for i = 1:k2_all
    miscol = [miscol isequal(dat(:,i), repmat(0,t2,1))];
end
miscol = logical(miscol);

if any(miscol)
    doNaN = 1;
    dat(:,miscol) = [];
    wh_val = ones(k2_all, k2_all);
    wh_val(miscol,:) = 0;
    wh_val(:,miscol) = 0;
    wh_val = logical(wh_val);
else
    doNaN = 0;
end

%% BASIC setting : Parameter

[t2, k2] = size(dat);

dccP = 1;
dccQ = 1;
archP = repmat(1, 1, k2);
garchQ = repmat(1, 1, k2);


%% Main Function
        
% Now let's do the univariate garching using fattailed_garch as it's faster then garchpq

stdresid = dat;
options = optimset('fmincon');
options = optimset(options, 'Display', 'off', 'Diagnostics', 'off', 'MaxFunEvals', 1000*max(archP+garchQ+1), ...
    'MaxIter', 1000*max(archP+garchQ+1), 'LargeScale', 'off', 'MaxSQPIter', 1000);
for i=1:k2
    if doverbose, fprintf(1,'Estimating GARCH model for Series %d\n',i); end
    try
        [univariate{i}.parameters, ~, ~, ~, univariate{i}.ht, ~] = ...
            fattailed_garch(dat(:,i), archP(i), garchQ(i), 'NORMAL', [], options);
    catch
        error('DCC does not converge: "fattailed_garch"');
    end
    stdresid(:,i) = dat(:,i)./sqrt(univariate{i}.ht);
end


dccstarting = [ones(1,dccP)*.01/dccP ones(1,dccQ)*.97/dccQ];
if doverbose, fprintf(1,'\n\nEstimating the DCC model\n'); end

try
    [dccparameters,~,~,~,~,~] = fmincon('dcc_mvgarch_likelihood', dccstarting, ones(size(dccstarting)), ...
        [1-2*options.TolCon], [], [], zeros(size(dccstarting))+2*options.TolCon, [], [], ...
        options, stdresid, dccP, dccQ);
catch
    error('DCC does not converge: "dcc_mvgarch_likelihood"');
end



% We now have all of the estimated parameters
parameters = [];
H = zeros(t2,k2);
for i= 1:k2
    parameters = [parameters; univariate{i}.parameters];
    H(:,i) = univariate{i}.ht;
end
parameters = [parameters; dccparameters'];



% We now have Ht and the likelihood

if ~doNaN
    [~, Ct, ~, ~] = dcc_mvgarch_full_likelihood(parameters, dat, archP, garchQ, dccP, dccQ);
elseif doNaN
    [~, Ct_temp, ~, ~] = dcc_mvgarch_full_likelihood(parameters, dat, archP, garchQ, dccP, dccQ);
    Ct = NaN(k2_all,k2_all,t2);
    Ct(repmat(wh_val, 1, 1, t2)) = Ct_temp;
end

% Mathematically same as...
%
% [~, Rt, ~, ~] = dcc_mvgarch_full_likelihood(parameters, dat, archP, garchQ, dccP, dccQ);
% Ht = zeros(k2,k2,t2);
% if ~doNaN
%     Ct = zeros(k2,k2,t2);
%     Hstd = H.^(0.5);
%     
%     for i= 1:t2
%         Ht(:,:,i) = diag(Hstd(i,:)) * Rt(:,:,i) * diag(Hstd(i,:));
%         Ct(:,:,i) = Ht(:,:,i) ./ sqrt(diag(Ht(:,:,i))*diag(Ht(:,:,i))');
%     end
%     
% elseif doNaN
%     Ct_temp = zeros(k2,k2,t2);
%     Hstd = H.^(0.5);
%     
%     for i= 1:t2
%         Ht(:,:,i) = diag(Hstd(i,:)) * Rt(:,:,i) * diag(Hstd(i,:));
%         Ct_temp(:,:,i) = Ht(:,:,i) ./ sqrt(diag(Ht(:,:,i))*diag(Ht(:,:,i))');
%     end
%     
%     Ct = NaN(k2_all,k2_all,t2);
%     Ct(repmat(wh_val, 1, 1, t2)) = Ct_temp;
% end

