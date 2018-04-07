function dat = cut_dat(dat, x, y, z)

% function dat = cut_dat(dat, x, y, z)
% 
% % example
% x = [-2 2];
% y = [-6 32];
% z = [30 60];
% dat = cut_dat(dat, x, y, z);


XYZ1 = [x(2) y(1) z(1)];
XYZ2 = [x(1) y(2) z(2)];

XYZout1 = mm2voxel(XYZ1, dat.volInfo.mat);
XYZout2 = mm2voxel(XYZ2, dat.volInfo.mat);

xyzlist = dat.volInfo.xyzlist;
idx = ones(size(xyzlist));

for i = 1:3
    idx(dat.volInfo.xyzlist(:,i) < XYZout1(i) | dat.volInfo.xyzlist(:,i) > XYZout2(i),i) = 0;
end

idx = sum(idx,2)==3;
dat.dat = dat.dat .* idx;
try 
    dat.sig = dat.sig .* idx;
catch
end

end