function [dn] = LVBin2DecDotProd(bn)
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
% [dn] = LVBin2DecDotProd([1 0 1 0 1 0 1 0])
dn = 2.^(length(bn)-1:-1:0)*bn';