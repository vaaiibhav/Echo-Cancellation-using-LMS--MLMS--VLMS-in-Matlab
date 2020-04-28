function N = LVxOrdEllipBPF(Rp,As,ws1,wp1,wp2,ws2)
% Computes the needed order to use in a call to the built-in function ellip
% when designing a bandpass filter meeting the specifications of Rp, As,
% and the band edges ws1,wp1,wp2, and ws2.
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
Fs = 1;
T = 1/Fs; 
OmS1 = (2/T)*tan(ws1*pi/2);
OmP1 = (2/T)*tan(wp1*pi/2);
OmP2 = (2/T)*tan(wp2*pi/2);
OmS2 = (2/T)*tan(ws2*pi/2);
Om0Sq = OmP1*OmP2;
OmPlp = (OmP2^2 - Om0Sq)/OmP2;
OmSlp1 = (OmS2^2 - Om0Sq)/OmS2;
OmSlp2 = (Om0Sq - OmS1^2)/OmS1;
OmSlp = min([OmSlp1,OmSlp2]);
OmgP =  OmPlp;
OmgS = OmSlp;
E = (10^(Rp/10)-1)^0.5;
A=10^(As/20); OmgC = OmgP;
k = OmgP/OmgS; k1 = E/sqrt(A^2-1);
[K1k, K2k] = ellipke([k, (1-k^2)].^2);
[K1k1, K2k1] = ellipke([k1, (1-k1^2)].^2);
N = ceil(K1k(1)*K1k1(2)/(K1k1(1)*K1k(2))); 