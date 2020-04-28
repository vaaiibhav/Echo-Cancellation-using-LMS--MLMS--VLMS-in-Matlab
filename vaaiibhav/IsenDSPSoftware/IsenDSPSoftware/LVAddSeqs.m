function [y, nOut] = LVAddSeqs(y1,n1,y2,n2)
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
nOut = [min(min(n1),min(n2)):1:max(max(n1),max(n2))];
mnfv = min(nOut); mxfv = max(nOut);
y = [zeros(1,min(n1)-mnfv) y1 zeros(1,mxfv-max(n1))] + ...
[zeros(1,min(n2)-mnfv) y2 zeros(1,mxfv-max(n2))];