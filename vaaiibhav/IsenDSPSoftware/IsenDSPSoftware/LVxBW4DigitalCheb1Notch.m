function [rS1,rS2,Rp1,Rp2] = LVxBW4DigitalCheb1Notch(H,OmP1,OmP2,As)
% OmP1 and OmP2 are the passband edges, in normalized frequency (i.e., in multiples of pi radians) 
% of a notch filter designed by the function cheby1, rS1 and rS2 are stopband edges at which a 
% desired value of stopband attenuation As is realized
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
HiFreqLim = 1; % i.e. 1*pi radians;
H = abs(H); H = H/max(H); LenH = length(H);
HdB = -20*log10(H+eps); y = find(HdB >= As);
indS1 = min(y); indS2 = max(y);
rS1 = (indS1/LenH)*HiFreqLim;
rS2 = (indS2/LenH)*HiFreqLim;
indP1 = round((OmP1/HiFreqLim)*LenH);
indP2 = round((OmP2/HiFreqLim)*LenH);
Rp1 = 20*log10(1/H(indP1));
Rp2 = 20*log10(1/H(indP2));
