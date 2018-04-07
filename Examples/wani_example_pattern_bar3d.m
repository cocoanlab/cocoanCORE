dat = fmri_data(which('weights_NSF_grouppred_cvpcr.img'));

x = [-2 2];
y = [-12 28];
z = [18 52];
dat = cut_dat(dat, x, y, z);
info = roi_contour_map(dat);
Z = info{1}.vZ;

close all;
Z = Z-min(min(Z))+.005;
Z = fliplr(Z);
maxZ = max(max(Z));

b = bar3(Z,1);

for k = 1:length(b)
    zdata = get(b(k), 'zdata');
    set(b(k), 'CData', zdata);
    set(b(k), 'FaceColor', 'interp');
    set(b(k), 'LineWidth', 0.5, 'edgecolor', [.2 .2 .2]);
end

% set(gca, 'zlim', [.005 maxZ]);

col = [0 0 0
    0.0461    0.3833    0.5912
    0.2461    0.5833    0.7912
    0.4000    0.7608    0.6471
    0.6706    0.8667    0.6431
    0.9020    0.9608    0.5961
    1.0000    1.0000    0.7490
    0.9961    0.8784    0.5451
    0.9922    0.6824    0.3804
    0.9569    0.4275    0.2627
    0.8353    0.2431    0.3098
    0.6196    0.0039    0.2588];

colormap(col);
view(-14, 78);
axis off;

set(gcf, 'color', 'w', 'position', [289     5   847   701]);

% take a snapshot
%% S2
dat = fmri_data(which('weights_NSF_grouppred_cvpcr.img'));

x = [34 66];
y = [-32 -2];
z = [15 17];
dat = cut_dat(dat, x, y, z);
info = roi_contour_map(dat);
Z = info{1}.vZ;

close all;
Z = Z-min(min(Z))+.005;
% Z = fliplr(Z);
maxZ = max(max(Z));
b = bar3(Z,1);

for k = 1:length(b)
    zdata = get(b(k), 'zdata');
    set(b(k), 'CData', zdata);
    set(b(k), 'FaceColor', 'interp');
    set(b(k), 'LineWidth', 0.5, 'edgecolor', [.2 .2 .2]);
end
% set(gca, 'zlim', [.005 maxZ]);
col = [0 0 0
    0.0461    0.3833    0.5912
    0.2461    0.5833    0.7912
    0.4000    0.7608    0.6471
    0.6706    0.8667    0.6431
    0.9020    0.9608    0.5961
    1.0000    1.0000    0.7490
    0.9961    0.8784    0.5451
    0.9922    0.6824    0.3804
    0.9569    0.4275    0.2627
    0.8353    0.2431    0.3098
    0.6196    0.0039    0.2588];

colormap(col);
view(21, 70);
axis off;

set(gcf, 'color', 'w', 'position', [289     5   847   701]);

%% dpINS
dat = fmri_data(which('weights_NSF_grouppred_cvpcr.img'));

x = [40 42];
y = [-23 16];
z = [-13 20];
dat = cut_dat(dat, x, y, z);
info = roi_contour_map(dat);
Z = info{1}.vZ;

close all;
Z = Z-min(min(Z))+.001;
Z = fliplr(Z);
maxZ = max(max(Z));
b = bar3(Z,1);

for k = 1:length(b)
    zdata = get(b(k), 'zdata');
    set(b(k), 'CData', zdata);
    set(b(k), 'FaceColor', 'interp');
    set(b(k), 'LineWidth', 0.5, 'edgecolor', [.2 .2 .2]);
end
% set(gca, 'zlim', [.005 maxZ]);
col = [0 0 0
    0.0461    0.3833    0.5912
    0.2461    0.5833    0.7912
    0.4000    0.7608    0.6471
    0.6706    0.8667    0.6431
    0.9020    0.9608    0.5961
    1.0000    1.0000    0.7490
    0.9961    0.8784    0.5451
    0.9922    0.6824    0.3804
    0.9569    0.4275    0.2627
    0.8353    0.2431    0.3098
    0.6196    0.0039    0.2588];

colormap(col);
view(-29, 78);
axis off;

set(gcf, 'color', 'w', 'position', [289     5   847   701]);