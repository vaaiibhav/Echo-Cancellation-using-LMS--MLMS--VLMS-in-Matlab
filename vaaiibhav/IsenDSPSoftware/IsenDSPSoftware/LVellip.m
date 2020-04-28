function [Z,P,K] = LVellip(N,Rp,As,OmgC)
% [Z,P,K] = LVellip(5,0.2,40,2)
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
[z,p,k] = ellipap(N,Rp,As);
Z = OmgC*z, P = OmgC*p,
nrmGainFrZero = prod(abs(z))/prod(abs(p)); 
nrmK = k*nrmGainFrZero;
UnrmGnFrZero = prod(abs(Z))/prod(abs(P)); 
K = nrmK/UnrmGnFrZero;
H = LVsFreqResp(K*poly(Z),poly(P),3*OmgC,18);

