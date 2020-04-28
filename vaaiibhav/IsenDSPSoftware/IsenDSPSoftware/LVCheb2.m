function LVCheb2(N,As,OmgC)
% LVCheb2(3,40,2)
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
[z,p,k] = cheb2ap(N,As)
z = OmgC*z, p = OmgC*p,
K = prod(abs(p))./prod(abs(z))
a = poly(p); b = K*poly(z);
H = LVsFreqResp(b,a,2*OmgC,16);
[Bbq,Abq,Gain] = LVDirToCascadeClassIIR(b/K,a,K)