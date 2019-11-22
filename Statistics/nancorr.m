function C = nancorr( X,Y )
% NANCORR calculates the sample correlation coefficient
%    for the series with NaNs expected.
%    X is the one series, Y is another.
% X=X(:);
% Y=Y(:);
L1=size(X,1);
L2=size(Y,1);

if L1 ~= L2
    error('The samples must be of the same length')
end

for i=1:L1
    if any(isnan(X(i,:)))
        Y(i,:)=NaN;
    end
    if any(isnan(Y(i,:)))
        X(i,:)=NaN;
    end
end
        
Xm=nanmean(X);
Ym=nanmean(Y);
C=nansum((X-Xm).*(Y-Ym))./sqrt((nansum((X-Xm).^2)).*(nansum((Y-Ym).^2)));
end