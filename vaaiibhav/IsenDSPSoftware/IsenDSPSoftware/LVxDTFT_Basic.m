function d = LVxDTFT_Basic(x,M,R)
% LVxDTFT_Basic([1 0 1],300,1)
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
N = length(x); W = exp(-j*R*pi/M); k = 0:1:M-1;
n = 0:1:N-1; dMat = W.^(n'*k); d = x*dMat;
