function [Z,P,K] = LVDesignCheby1Notch(OmP1,OmS1,OmS2,OmP2,Rp,As)
% [Z,P,K] = LVDesignCheby1Notch(4,5,8,10,1,40) 
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
Om0Sq = OmP1*OmP2; OmPlp = OmP1/(Om0Sq-OmP1^2);
OmSlp1 = OmS2/(OmS2^2 - Om0Sq);
OmSlp2 = OmS1/(Om0Sq - OmS1^2);   
OmSlp = min([OmSlp1,OmSlp2]); 
[Z1,P1,K1] = LVDesignCheby1Filter(Rp,As,OmPlp,OmSlp);
H = LVsFreqResp(K1*poly(Z1),poly(P1),2*OmSlp,19);  
Num = 1; Den = 1; D = [1 0 Om0Sq]; N = [0 1 0];
for Ctr = 1:1:length(P1)
Num = conv(Num,[N-P1(Ctr)*D]);
Den = conv(Den,D); end; 
a = real(Num); b = real(Den); s=j*[0:0.001:2*OmS2]; 
Hs = abs(polyval(b,s)./polyval(a,s)); K = 1/max(Hs); 
Z = roots(b); P = roots(a); H = LVsFreqResp(K*b,a,2*OmS2,31);
[NetAs,NetRp1,NetRp2] = LVxRealizedFiltParamNotch(H,OmP1,OmS1,OmS2,OmP2,2*OmS2)