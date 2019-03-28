function Ct = DCC_garch_jj(dat, doverbose)
% Modified version of dcc_mvgarch.m
% Simplify output

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
        
% Now lest do the univariate garching using fattailed_garch as it's faster then garchpq

stdresid = dat;
options = optimset('fmincon');
options = optimset(options, 'Display', 'off', 'Diagnostics', 'off', 'MaxFunEvals', 1000*max(archP+garchQ+1), ...
    'MaxIter', 1000*max(archP+garchQ+1), 'LargeScale', 'off', 'MaxSQPIter', 1000);
for i=1:k2
    if doverbose, fprintf(1,'Estimating GARCH model for Series %d\n',i); end
    [univariate{i}.parameters, ~, ~, ~, univariate{i}.ht, ~, dcc_exit] = ...
        fattailed_garch_CAPS2(dat(:,i), archP(i), garchQ(i), 'NORMAL', [], options);
    stdresid(:,i) = dat(:,i)./sqrt(univariate{i}.ht);
    if dcc_exit < 1
        error('DCC does not converge. : fattailed_garch_CAPS2');
    end
end


dccstarting = [ones(1,dccP)*.01/dccP ones(1,dccQ)*.97/dccQ];
if doverbose, fprintf(1,'\n\nEstimating the DCC model\n'); end

[dccparameters,~,dcc_exit,~,~,~] = fmincon('dcc_mvgarch_likelihood', dccstarting, ones(size(dccstarting)), ...
    [1-2*options.TolCon], [], [], zeros(size(dccstarting))+2*options.TolCon, [], [], ...
    options, stdresid, dccP, dccQ);

if dcc_exit < 1
    error('DCC does not converge. : dcc_mvgarch_liklihood');
end



% We now have all of the estimated parameters
parameters = [];
H = zeros(t2,k2);
for i= 1:k2
    parameters = [parameters; univariate{i}.parameters];
    H(:,i) = univariate{i}.ht;
end
parameters = [parameters; dccparameters'];



%We now have Ht and the likelihood
if ~doNaN
    [~, Ct, ~, ~] = dcc_mvgarch_full_likelihood(parameters, dat, archP, garchQ, dccP, dccQ);
elseif doNaN
    [~, Ct_temp, ~, ~] = dcc_mvgarch_full_likelihood(parameters, dat, archP, garchQ, dccP, dccQ);
    Ct = NaN(k2_all,k2_all,t2);
    Ct(repmat(wh_val, 1, 1, t2)) = Ct_temp;
end
    
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

