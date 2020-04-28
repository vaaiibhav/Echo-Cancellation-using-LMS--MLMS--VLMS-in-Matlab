function LVButterFRViaPoly(N,OmegaC,HiFreqLim)
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
%
% LVButterFRViaPoly(8,1,2)
k = 0:1:2*N-1; 
P = OmegaC*exp(j*(pi/(2*N))*(2*k+N+1))
NetP = P(find(real(P)<0)); b = 1; a = poly(NetP); 
H = LVsFreqResp(b,a,HiFreqLim,12);