function [Z,P,K] = LVDesignCheby1Filter(Rp,As,OmgP,OmgS)
% [Z,P,K] = LVDesignCheby1Filter(0.5,40,0.5,0.65)
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
E = sqrt(10^(Rp/10)-1); A = 10^(As/20); OmgC = OmgP; 
OmgT = OmgS/OmgP; g = sqrt((A^2-1)/(E^2));
N = ceil(log10(g + sqrt(g^2 - 1))/log10(OmgT + sqrt(OmgT^2-1)))
k = -(N-1):2:N-1; r = k*pi/(2*N); v0 = asinh(1/E)/N;
P = OmgC*(-sinh(v0)*cos(r) + j*cosh(v0)*sin(r));
K = prod(abs(P)); Z = []; if rem(N,2)==0 % N even
K = 1/sqrt(1 + E^2)*K; end
