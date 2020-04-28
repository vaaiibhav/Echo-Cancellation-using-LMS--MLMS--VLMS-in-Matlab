function [Z,P,K] = LVCheb1Lpf2Hpf(N,Rp,WpLP,WpHP)
% [Z,P,K] = LVCheb1Lpf2Hpf(7,1,3,3)
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
%
[z,p,k] = cheb1ap(N,Rp);
P = WpLP*p; G = k*WpLP^(length(P));
b1 = G; a1 = poly(P);
Num = 1; Den = 1; D = [1 0];
c = WpLP*WpHP; N = [0 c];
for Ctr = 1:1:length(P)
Num = conv(Num,N-P(Ctr)*D);
Den = conv(Den,D); end
a = real(Num); Scale = 1/a(1);
b2 = real(Den); a2 = Scale*a;
K = Scale*G; Z = roots(b2); P = roots(a2); 
LVsFreqRespDouble(b1,a1,2*WpHP,29,K*b2,a2)