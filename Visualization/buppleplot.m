function h = buppleplot(x, cols)

xyr_all = NaN(round(sum(x)/10),3);
n_filled = 0;
frames = 25;
n = round(sum(x)/10);

while n_filled < n
    frames = frames + 1;
    
    for i = 1:20
        iter_n = frames^2;
        xyr = [rand(iter_n,2)*frames ones(iter_n,1)];
        dist_xy = tril(squareform(pdist(xyr(:,1:2))));
        
        dist_xy(triu(true(size(dist_xy)))) = inf;
        isOK = all(dist_xy >= 2*2, 2);
    end
    
    xyr_all(n_filled+1:n_filled+sum(isOK),:) = xyr(isOK,:);
    n_filled = n_filled + sum(isOK);
end

scatter(xyr_all(:,1), xyr_all(:,2), xyr_all(:,3)*100)


% & sqrt(xyr(size(circData,1)+1:end,1).^2+xyr(size(circData,1)+1:end,2).^2)<20;

end