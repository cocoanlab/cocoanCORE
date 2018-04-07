function pagesetup(handle)
% function pagesetup(handle)
% setup pager size as figure size
    sz = get(handle, 'Position'); sz = sz(3:4);
    set(handle,'Paperposition', [0 0 sz], 'PaperUnits', 'points', 'PaperSize', sz);
end