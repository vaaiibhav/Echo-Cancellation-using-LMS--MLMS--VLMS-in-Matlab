function [Z,P,K] = LVDesignCheby1BPF(OmS1,OmP1,OmP2,OmS2,Rp,As)
% [Z,P,K] = LVDesignCheby1BPF(4,5,8,10,1,40) 
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
    Om0Sq = OmP1*OmP2;
    OmPlp = (OmP2^2 - Om0Sq)/OmP2
    OmSlp1 = (OmS2^2 - Om0Sq)/OmS2;
    OmSlp2 = (Om0Sq - OmS1^2)/OmS1;
    OmSlp = min([OmSlp1,OmSlp2])
    [Z1,P1,K1] = LVDesignCheby1Filter(Rp,As,OmPlp,OmSlp);
    Num = 1; Den = 1; 
    N = [1 0 Om0Sq]; D = [0 1 0];
    for Ctr = 1:1:length(P1)
    Num = conv(Num,[N-P1(Ctr)*D]);
    Den = conv(Den,D); end
    a = real(Num);  b = real(Den); s=j*[0:0.001:2*OmS2];
    Hs = abs(polyval(b,s)./polyval(a,s)); 
    K = 1/max(Hs);   Z = roots(b); P = roots(a); 
    LVsFreqRespDouble(K1*poly(Z1),poly(P1),2*OmS2,30,K*b,a)