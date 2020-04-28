function LVCheby1PolesAndSysFcn(N,OmC,Rp)
% LVCheby1PolesAndSysFcn(3,4,0.2)
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
[z,p,k] = cheb1ap(N,Rp),
P = OmC*p, K = k*OmC^N,
b = real(poly(z)), a = real(poly(P)),
H = LVsFreqResp(K,poly(P),2*OmC,19);
[Bbq,Abq,Gain] = LVDirToCascadeClassIIR(b,a,K)    
[b,a,k]=LVCas2DirClassIIR(Bbq,Abq,Gain)
    
