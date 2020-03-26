function [R,T,p,df] = dcor(X,Y)

% [R,T,p,df] = dcor(X,Y)
% Computes the U-centered (bias corrected) distance correlation between X and Y. Also outputs the
% associated t-value (T) and p-value (p). Rows represent the examples, and columns the variables.
%
% Based on: http://www.mathworks.com/matlabcentral/fileexchange/49968-dcorr--x--y--
% and, the R package energy and papers by Szekely and Rizzo (2007; 2013 and 2014). 
%
% Author: Linda Geerligs (lindageerligs@gmail.com), Date: 22-10-2015
% 
% Obtained from https://github.com/MRC-CBU/riksneurotools/blob/master/Conn/dcor_uc.m

a = pdist2(X, X);
b = pdist2(Y, Y);
n=size(X,1);

A = Ucenter(a,n);
B = Ucenter(b,n);

dcovXY = sum(sum(A.*B)) ./ (n*(n-3));
dvarX = sum(sum(A.*A)) ./ (n*(n-3));
dvarY = sum(sum(B.*B)) ./ (n*(n-3));

R=dcovXY / sqrt(dvarX * dvarY);

df=(n*(n-3))/2 -1;
T=sqrt(df).*(R./sqrt(1-R.^2));
p = tcdf(T,df,'upper');

if R<0
    R=0;
else
    R=sqrt(R);
end

    function A = Ucenter(a,n)
        m = sum(a,2);
        M = sum(m)/((n - 1) * (n - 2));
        m = m./(n-2);
        A = a - repmat(m,[1 n]);
        A = A - repmat(m,[1 n])';
        A = A+M;
        A(eye(size(A))==1)=0;
    end

end