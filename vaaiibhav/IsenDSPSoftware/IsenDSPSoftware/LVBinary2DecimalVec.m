function [dn] = LVBinary2DecimalVec(bn)
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
%
%  [dn] = LVBinary2DecimalVec([1 0 1 0 1 0 1 0])
dn = sum(2.^(length(bn)-1:-1:0).*bn);