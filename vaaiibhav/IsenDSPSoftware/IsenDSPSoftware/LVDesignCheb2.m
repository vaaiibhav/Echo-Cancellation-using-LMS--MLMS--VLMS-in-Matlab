function [Z,P,K] = LVDesignCheb2(Rp,As,OmgP,OmgS)
% [Z,P,K] = LVDesignCheb2(0.2,40,0.9,1)
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
OmgC = 1; OmgP = OmgP/OmgS;
E = 1/sqrt(10^(As/10)-1); G = 10^(-Rp/20);
N = acosh(G/(E*sqrt(1-G^2)))/acosh(1/OmgP);
N = ceil(abs(N)); V0 = asinh(1/E)/N;
k = -(N-1):2:(N-1); r = k*pi/(2*N);
Ch1P = -sinh(V0)*cos(r) + j*cosh(V0)*sin(r);
P = OmgS*(1./Ch1P); div = sin(k*pi/(2*N)); 
Zdiv = find(div==0); NZerDiv = div(find(~(div==0))); 
Z = OmgS*(j./NZerDiv); s=j*[0:0.001:2*OmgS];
Hs = abs(polyval(poly(Z),s)./polyval(poly(P),s)); 
K = 1/max(abs(Hs));