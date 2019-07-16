function [offset, stats] = hex_glm1(y, angle, varargin)

% [offset, stats] = hex_glm1(y, angle, optional input)
%
% angle should be in degree
% 

folds = 6;

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            % functional commands
            case {'folds', 'fold'}
                folds = varargin{i+1};
        end
    end
end

% design matrix
x = zeros(length(angle), 3);
x(:,1) = ones(size(x,1),1);
x(:,2) = cos(deg2rad(folds*angle));
x(:,3) = sin(deg2rad(folds*angle));

% multiple regression
%[b, ~, stats] = glmfit(x, y);
b = (pinv(x)*y)';

% calulate model fit
% SS_resid = sum((stats.resid).^2);
mean_y = repmat(mean(y), size(y,1),1);

SS_total = sum((y-mean_y).^2);
SS_reg = sum((x*b' - mean_y).^2);

% SS_reg = SS_total - SS_resid;
df_reg = size(x,2)-1;
df_resid = size(x,1)-size(x,2);

MSR = SS_reg./df_reg;
MSE = (SS_total-SS_reg)./df_resid;

stats.F = MSR./MSE;

stats.f2p = 1-fcdf(stats.F,df_reg, df_resid);

stats.f2z = f2z(stats.F, df_reg, df_resid);

stats.r_squared = 1 - (SS_total-SS_reg)./SS_total;

% calulate offset

offset_init = atan2d(b(:,3), b(:,2));
offset_init(offset_init<0) = offset_init(offset_init<0)+360;

offset = (offset_init./folds)';


end
