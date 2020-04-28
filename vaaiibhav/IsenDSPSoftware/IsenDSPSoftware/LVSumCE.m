function s = LVSumCE(k,N)
% s = LVSumCE(2,32)
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
n = 0:1:N-1; s = sum(exp(j*2*pi*n*k/N));