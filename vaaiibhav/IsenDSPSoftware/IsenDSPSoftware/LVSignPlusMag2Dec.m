function [dn] = LVSignPlusMag2Dec(bn)
%  [dn] = LVSignPlusMag2Dec([1 0 1 0 1 0 1 0])
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
pbn = fliplr(bn(2:length(bn)));
dn = (-1)^(bn(1))*sum(2.^(0:1:length(pbn)-1).*pbn);