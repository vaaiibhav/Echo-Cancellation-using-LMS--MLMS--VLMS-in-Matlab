function b = LVIdealLPFImpResp(wc,L)
% LVIdealLPFImpResp(0.25*pi,61)
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
M = (L-1)/2; n = 0:1:L-1;
b = sin(wc*(n - M + eps))./(pi*(n - M + eps));
figure(55); stem(n,b)