function [gcf_position] = setfig(varargin)
% change the figure handle for
% size, width/length ratio, and color of the figure window
% Esp, 'square' option will make the window as square. 
%
% :Optional inputs: (with their default values)
% ::
%   - 'gcf'                 change the gcf 
%   - 'gca'                 change the gca 
%   - 'square'              make the figure square
%   - 'rectangle'           make the figure rectangle
%   - 'color'               change the window color, default is white 
%   - position              position you want to change (1 x 4)
%   - window_color          color of the background want to change (string or 1 x 3)
%   - ratio                 width/length of window want to change 
%
% :Output:
% ::
%   - 'gcf_position'        new gcf position
%
% :Example:
% ::
%     setfig('gcf', gcf_position, 'color', window_color);
%     setfig('square');
%     setfig('rectangle', ratio);
%     setfig('gcf', get(gcf,'position'));
%
%
% Copyright Byeol Kim, 2020
%
set_gcf = false;
set_gca = false;
change_ratio = false;
window_color = 'w'; % default is white
gcf_position = get(gcf, 'Position');
gca_position = get(gca, 'Position');

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            % functional commands
            case {'gcf'}
                set_gcf = true;
                gcf_position = varargin{i+1};
            case {'gca'}
                set_gca = true;
                gca_position = varargin{i+1};
            case {'square'}
                change_ratio = true;
                ratio = 1;
            case {'rectangle'}
                change_ratio = true;
                ratio = varargin{i+1};
            case {'color'}
                window_color = varargin{i+1};
        end
    end
end

if change_ratio
    gcf_position(3) = gcf_position(4) * ratio;
    if ratio == 1
        gca_position = [min(gca_position(1:2)), min(gca_position(1:2)), max(gca_position(3:4)), max(gca_position(3:4))];
    end
end

set(gcf, 'Position', gcf_position, 'color', window_color)
set(gca, 'Position', gca_position, 'color', window_color)

end
