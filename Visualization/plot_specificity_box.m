function out = plot_specificity_box(yfit1, yfit2, varargin)

xlim = [-.1 1.2];
ylim = [-.1 .4];

cols = [0.8353    0.2431    0.3098
    0.9922    0.6824    0.3804];

savefig = false;

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            % functional commands
            case {'xlim'}
                xlim = varargin{i+1};
            case {'ylim'}
                ylim = varargin{i+1};
            case {'colors', 'color'}
                cols = varargin{i+1};
        end
    end
end

create_figure('plot');

data = [yfit1 yfit2];

close all;
boxplot_wani_2016(data, 'color', cols, 'linewidth', 2, 'boxlinewidth', 1, 'mediancolor', 'k', 'violin');

xdot{1} = ones(size(data,1),1)*1+.32;
xdot{2} = ones(size(data,1),1)*2-.32;

wh = data(:,1) > data(:,2);

out.h1 = line([xdot{1}(wh) xdot{2}(wh)]', data(wh,1:2)', 'color', [227,26,28]./255, 'linewidth', 1.5);
out.h2 = line([xdot{1}(~wh) xdot{2}(~wh)]', data(~wh,1:2)', 'color', [0.1961    0.5333    0.7412], 'linewidth', 1.5);

scatter(xdot{1}, data(:,1), 20, cols(1,:), 'filled'); 
scatter(xdot{2}, data(:,2), 20, cols(2,:), 'filled'); 

set(gcf, 'position', [1   518   204   187]);
set(gca, 'fontsize', 18, 'xlim', [0.5 2.5], 'linewidth', 1.5, 'ticklength', [.03 .03], 'xticklabel', '');

end