function histo_wani(dat)

clf;
N = size(dat.dat,2);
for i=1:size(dat.dat,2)
    dattmp = dat.dat(:,i);
    for j=1:100
        if N > (j+1)*j
            j=j+1;
        else
            k=j;
            break
        end
    end
    subplot(k,k+1,i);
    [h, x] = hist(dattmp, 100);
    han = bar(x, h);
    set(han, 'FaceColor', [.3 .3 .3], 'EdgeColor', 'none');
    % axis([-1 1 0 12000]);
    % axis([-1 1]);
    
    xlabel('Values'); ylabel('Frequency');
    hist_title = ['Histogram of values ' num2str(i)];
    title(hist_title);
    drawnow
        
    clear dattmp
    
    globalmean = nanmean(dat.dat);  % global mean of each obs
    globalstd = nanstd(dat.dat);  % global mean of each obs
    nobs = length(globalmean);
    sz = rescale_range(globalstd, [1 6]); % marker size related to global std
    %sz(sz < .5) = .5;
end

end

function rx = rescale_range(x, y)
% re-scale x to range of y
m = range(y)./range(x);

if isinf(m)
    % no range/do not rescale
    rx = x;
else
    x = x - min(x);
    rx = y(1) + x * ((y(2) - y(1)) ./ max(x));
end
end