function [CorC] = LVCorrCosinesZerothLag(k1,k2,N)
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
n = 0:1:N-1; 
CorC = sum(cos(2*pi*n*k1/N).*cos(2*pi*n*k2/N));
