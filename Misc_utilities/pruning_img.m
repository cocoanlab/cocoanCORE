function [dat, r] = pruning_img(stat_img, p, k)

% [dat, r] = pruning_img(stat_img, p, k)
%
% stat_img: statistical image
% p: multiple p values [.05 .01 .001]
% k: multiple k values [1 1 10]

stat_img.dat(stat_img.dat > 0) = 1;
stat_img.dat(stat_img.dat < 0) = -1;

% ================= pruning  ====================================== 
stat_img = threshold(stat_img, p(1), 'unc', 'k', k(1));
final_dat = stat_img.dat .* stat_img.sig;

if length(p) > 1
    for i = 2:length(p)
        stat_img = replace_empty(stat_img);
        stat_img = threshold(stat_img, p(i), 'unc', 'k', k(i));
        final_dat((stat_img.dat .* stat_img.sig) == 1) = i;
        final_dat((stat_img.dat .* stat_img.sig) == -1) = -i;
    end
end

% ================= prepare the frame ===================================== 
stat_img = replace_empty(stat_img);
stat_img = threshold(stat_img, p(1), 'unc', 'k', k(1));
stat_img.dat = final_dat .* stat_img.sig; 


% ================= remove regions that do not have 3 or -3================
r = region(stat_img);

k = [];
for i = 1:numel(r)
    if any(abs(r(i).val) == 3)
        k = [k i];
    end
end

fprintf('\ndeleting... remaining # of regions = %d', numel(k));

r = r(k);
for i = 1:numel(r)
    r(i).Z = r(i).val';
end

dat = region2imagevec(r);

end