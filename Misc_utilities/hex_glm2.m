function stats = hex_glm2(y, angle, offset, varargin)

% stats = hex_glm2(y, angle, offset, varargin)
%
% y: #trial x #roi
% angle: #trial x 1
% offset: 1 x #roi
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
n_roi = size(y,2);
x = cell(n_roi,1);

for i = 1:n_roi
    x{i} = zeros(length(angle), 2);
    x{i}(:,1) = ones(size(x{i},1),1);
    x{i}(:,2) = cos(deg2rad(folds*(angle-offset(i)))); 
end


for i = 1:n_roi
    
    % multiple regression for 6 folds
    y_i = y(:,i);
    stats.b{i} = (pinv(x{i})*y_i)';
    
    % calulate model fit
    % SS_resid = sum((stats.resid).^2);
    mean_y = repmat(mean(y_i), size(y_i,1),1);
    
    SS_total = sum((y_i-mean_y).^2);
    SS_reg = sum((x{i}*stats.b{i}' - mean_y).^2);
    
    % SS_reg = SS_total - SS_resid;
    df_reg = size(x{i},2)-1;
    df_resid = size(x{i},1)-size(x{i},2);
    
    MSR = SS_reg./df_reg;
    MSE = (SS_total-SS_reg)./df_resid;
    
    stats.F{i} = MSR./MSE;
    
    stats.f2p{i} = 1-fcdf(stats.F{i},df_reg, df_resid);
    
    stats.f2z{i} = f2z(stats.F{i}, df_reg, df_resid);
    
    stats.r_squared{i} = 1 - (SS_total-SS_reg)./SS_total;
    
    
    % build fold regressors (e.g., 12 regressors for 6 folds)
    
    interval = 360/(folds*2);
    deg = -(interval/2):interval:(360-(interval/2));
    angle_offset = angle-offset(i);
    angle_offset(angle_offset<-15)=angle_offset(angle_offset<-15)+360;
    
    clear xx1; 
    
    for j = 1:numel(deg)-1
        xx1(:,j) = double(deg(j)<(angle_offset) & (angle_offset)<=deg(j+1));
    end
    
    xx1(:,end+1) = ones(size(xx1,1),1);
    
    stats.b_align{i} = (pinv(xx1)*y(:,i))';
    
    stats.b_align{i} = stats.b_align{i}(1:end-1);
    
end

end
