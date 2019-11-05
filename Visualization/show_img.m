function show_img(base_img, overlay_img, cmap, clim)

n_img = numel(overlay_img);
spm_check_registration(repmat(base_img, n_img, 1));
myfig = findobj('Tag','Graphics'); % colormap
colormap(myfig, cmap);
for i = 1:n_img
    overlay_spm = spm_vol(overlay_img{i});
    overlay_dat = spm_read_vols(overlay_spm);
    M = overlay_spm.mat;
    [xyz, ~, k] = find(overlay_dat(:));
    [x, y, z] = ind2sub(overlay_spm.dim, xyz);
    spm_orthviews('AddBlobs', i, [x y z]', k, M, 'img');
    spm_orthviews('redraw');
end

global st
for i = 1:n_img
    st.vols{i}.blobs{1}.min = clim{i}(1);
    st.vols{i}.blobs{1}.max = clim{i}(2);
    spm_orthviews('redraw');
end

end