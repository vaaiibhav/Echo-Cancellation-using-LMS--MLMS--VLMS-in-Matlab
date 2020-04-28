function LVs2zViaBilinearEx(N,Rp,OmegaC,Fs,CheborButter) 
% LVs2zViaBilinearEx(3,1,1,1,1)
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
if CheborButter==1
[Z,P,K] = cheb1ap(N,Rp);
else; [Z,P,K] = buttap(N); end
K = K*OmegaC^N; P = OmegaC*P; 
sB = K; sA = poly(P);
Num = 1; Den = 1; D = [1 1];
for Ctr = 1:1:length(P)
    N = [(2*Fs - P(Ctr)), -(2*Fs + P(Ctr))];
    Num = conv(Num,N); Den = conv(Den,D); 
end; 
zB = real(Den); zA = real(Num);
LVsFRzFr(sB,sA,OmegaC,K*zB,zA,37)


