function [xe,xo,m] = LVEvenOddSymmZero(x,n)
% [xe,xo,m] = LVEvenOddSymmZero([1 2 3 4 5],[0 1 2 3 4])
% [xe,xo,m] = LVEvenOddSymmZero([1 2 3 4 5],[-1 0 1 2 3])
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
m = -fliplr(n); m1=min([m,n]); m2 = max([m,n]); m = m1:m2;
nm = n(1)-m(1); n1 = 1:1:length(n); x1 = zeros(1,length(m));
x1(n1+nm) = x; xe = 0.5*(x1 + fliplr(x1));
xo = 0.5*(x1 - fliplr(x1));