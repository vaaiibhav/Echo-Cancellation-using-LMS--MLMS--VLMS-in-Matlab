function LVButterPolesAndSysFcn(N,OmC)
% LVButterPolesAndSysFcn(3,3)
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
%
% N = 3; OmC = 3; 
[z,p,k] = buttap(N),
P = OmC*p, K=k*OmC^N,
b = real(poly(z)), a = real(poly(P)),
[Bbq,Abq,Gain] = LVDirToCascadeClassIIR(b,a,K)