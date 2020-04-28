function [dn] = LVBinary2Decimal(bn)
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
%  [dn] = LVBinary2Decimal([1 0 1 0 1 0 1 0])
Lenbn = length(bn); dn = 0; for ctr = Lenbn:-1:1; 
dn = dn + 2^(ctr-1)*bn(Lenbn-ctr+1); end;