function angle = get_angle(xy)

% Calculate trajectory angles in a 2D space using xy coordinates
%
% :Usage:
% ::
%
%    angle = get_degree(xy)
%
% :Inputs:
%
%   **xy:**
%        xy coordinates, (n x 2) matrix, where n is the number of trials in
%        a sequential order. xy can be a cell array. 
%
% :Output:
%
%   **angle:**
%        angle (in degree) of the movement trajectory from t to t+1 trials
%        on a 2d space. The reference (0 degree) is the (1,0) unit vector.
%
% 
% ..
%    Copyright (C) 2019  Jiwon Jeon, Mijin Kim, Suhye Kim, & Wani Woo
% ..

% Programmers' notes:
% 

if iscell(xy)
    for i = 1:numel(xy)
        xy_diff = xy{i}(2:end,:)-xy{i}(1:end-1,:);
        angle{i} = atan2d(xy_diff(:,2), xy_diff(:,1));
        angle{i}(angle{i}<0)=angle{i}(angle{i}<0)+360;
    end
else
    xy_diff = xy(2:end,:)-xy(1:end-1,:);
    angle = atan2d(xy_diff(:,2), xy_diff(:,1));
    angle(angle<0)=angle(angle<0)+360;
end

end