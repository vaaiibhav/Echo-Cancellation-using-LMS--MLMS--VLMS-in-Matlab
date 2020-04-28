function [y,nOut] = LVMultSeqs(y1,n1,y2,n2)
% y1 is first sequence and n1 is its sample index or position vector, 
% y2 is first sequence and n2 is its sample index or position vector.
% [y,nOut] = LVMultSeqs([1 2 3 4],[-1 0 1 2],[6 5 4 3 2 1],[-3 -2 -1 0 1 2])
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
nOut = [min(min(n1),min(n2)):1:max(max(n1),max(n2))];
mnfv = min(nOut); mxfv = max(nOut);
y = [zeros(1,min(n1)-mnfv) y1 zeros(1,mxfv-max(n1))].* ...
    [zeros(1,min(n2)-mnfv) y2 zeros(1,mxfv-max(n2))];