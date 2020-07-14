function out = plot_y_yfit(yval, yfit, varargin)

% out = plot_y_yfit(yval, yfit, varargin)
%
% case {'xlim'}
%     xlim = varargin{i+1};
% case {'ylim'}
%     ylim = varargin{i+1};

position = [1   747   218   208];

prop = 0.2;
if iscell(yval)
    if ~any(isrow(yval)-isrow(yval{1}))
        xlim = [min(min(cell2mat(yval)))-abs(min(min(cell2mat(yval))))*prop max(max(cell2mat(yval)))+max(max(cell2mat(yval)))*prop];
        ylim = [min(min(cell2mat(yfit)))-abs(min(min(cell2mat(yfit))))*prop max(max(cell2mat(yfit)))+max(max(cell2mat(yfit)))*prop];
    else
        xlim = [min(min(cell2mat(yval')))-abs(min(min(cell2mat(yval'))))*prop max(max(cell2mat(yval')))+max(max(cell2mat(yval')))*prop];
        ylim = [min(min(cell2mat(yfit')))-abs(min(min(cell2mat(yfit'))))*prop max(max(cell2mat(yfit')))+max(max(cell2mat(yfit')))*prop];
    end
else
    xlim = [min(min(yval))-abs(min(min(yval)))*prop max(max(yval))+max(max(yval))*prop];
    ylim = [min(min(yfit))-abs(min(min(yfit)))*prop max(max(yfit))+max(max(yfit))*prop];
end

data_alpha  = 1;
line_alpha  = 1;
dotsize = 40;
do_random = 0;
do_xyline = false;

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            % functional commands
            case {'xlim'}
                xlim = varargin{i+1};
            case {'ylim'}
                ylim = varargin{i+1};
            case {'position'}
                position = varargin{i+1};
            case {'data_alpha'}
                data_alpha = varargin{i+1};
            case {'line_alpha'}
                line_alpha = varargin{i+1};
            case {'dotsize'}
                dotsize = varargin{i+1};
            case {'random_col'}
                do_random = 1;
            case {'xyline'}
                do_xyline = true;
        end
    end
end

%% plotting
if do_random
    %     colors = [0.9098    0.4902    0.4471
    %         0.7451    0.6039    0.2000
    %         0.4275    0.6863    0.2039
    %         0.3373    0.7373    0.5922
    %         0.3059    0.7098    0.9020
    %         0.6235    0.5608    0.9725
    %         0.9137    0.4353    0.8235];
    colors = [166,206,227
        31,120,180
        178,223,138
        51,160,44
        251,154,153
        227,26,28
        253,191,111
        255,127,0
        202,178,214
        106,61,154
        177,89,40]./255;
else
    colors = [255,237,160
        254,217,118
        254,178,76
        253,141,60
        252,78,42
        227,26,28
        189,0,38]./255;
end

create_figure('y_yfit_plot');

clear test_por;
if iscell(yval)
    for i = 1:numel(yval)
        x = yval{i};
        y = yfit{i};
        b = glmfit(x,y);
        test_por(i) = corr(x,y);
    end
else
    for i = 1:size(yval,2)
        x = yval(:,i);
        y = yfit(:,i);
        b = glmfit(x,y);
        test_por(i) = corr(x,y);
    end
end

if do_random
    colors = repmat(colors, size(yval,2), 1);
    k = randperm(size(yval,2));
else
    dif = 1/size(colors,1);
    
    k = zeros(size(test_por));
    for i = 1:size(colors,1)
        idx = test_por <= (dif*i+.0001) & test_por >= dif*(i-1);
        k(idx) = i;
    end
end

%%
marker_shapes = repmat('osd^v><', 1, 40);

if iscell(yval)
    aa = numel(yval);
else
    aa = size(yval,2);
end

for i = 1:aa
    hold on;
    
    if iscell(yval)
        x = yval{i};
        y = yfit{i};
    else
        x = yval(:,i);
        y = yfit(:,i);
    end
    b = glmfit(x,y);
    out.b(i,1) = b(2);
    out.r(i,1) = corr(x,y);
    try
        line_h(i) = line(xlim, b'*[ones(1,2); xlim], 'linewidth', 1.5, 'color', colors(k(i),:)); % cmap(round(i*1.5),:));
        line_h(i).Color(4) = line_alpha;
        h = scatter(x, y, dotsize, colors(k(i),:), 'filled', 'markerfacealpha', data_alpha, 'marker', marker_shapes(i));
    catch
        line_h(i) = line(xlim, b'*[ones(1,2); xlim], 'linewidth', 1.5, 'color', [0.1961    0.5333    0.7412]); % cmap(round(i*1.5),:));
        line_h(i).Color(4) = line_alpha;
        h = scatter(x, y, dotsize, [0.1961    0.5333    0.7412], 'filled', 'markerfacealpha', data_alpha, 'marker', marker_shapes(i));
    end 
end

if do_xyline
    line(xlim, xlim, 'linewidth', 4, 'linestyle', ':', 'color', [.5 .5 .5]);
end

set(gcf, 'position', position);
set(gca, 'tickdir', 'out', 'TickLength', [.03 .03], 'linewidth', 1.5, 'xlim', xlim, 'ylim', ylim, 'fontsize', 18);

end