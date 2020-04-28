function [actRp,actAs,WF] = LVxDesignEquirippHPF(Rp,As,ws,wp) 
% function [actRp,actAs,WF] = LVxDesignEquirippHPF(Rp,As,ws,wp)
% Rp and As are the maximum passband ripple and minimum stopband attenuation, respectively.
% ws and wp are the stopband and passband edges, specified in normalized frequency, i.e., units of pi.
% actRp and actAs are the realized values of Rp and As from the filter design.
% Test calls:
% LVxDesignEquirippHPF(0.2,60,0.45,0.55)
% LVxDesignEquirippHPF(0.02,60,0.4,0.55)
% LVxDesignEquirippHPF(0.5,70,0.15,0.25)
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool

Rfac= 10^(-Rp/20); DeltaP = (1-Rfac)/(1+Rfac);
DeltaS = (1+DeltaP)*10^(-As/20);
deltaF = abs(ws - wp)/2;
compL = (-20*log10(sqrt(DeltaP*DeltaS)) - 13)/(14.6*deltaF) + 1;
L = round(compL), Ord = L - 2; LenFFT = 2^15;
SBAtten = 0; PBR = 100;
% make Ord even so filter length is constrained to be odd
if ~(rem(Ord,2)==0)
    Ord = Ord + 1;
end
while (SBAtten < As)|(PBR > Rp); Ord = Ord + 2;
WF = firpm(Ord,[0,ws,wp,1],[0,0,1,1],[DeltaS/DeltaP,1]);
fr = abs(fft(WF,LenFFT));
fr = fr(1,1:LenFFT/2+1)/(max(abs(fr)));
Lfr = length(fr); SB = fr(1,1:round(ws*(Lfr-1)));
PB = fr(1,round(wp*(Lfr-1))+1:Lfr);
PBR = -20*log10(min(PB) + eps)
SBAtten = -20*log10(max(SB) + eps)
end; 
Fin_L = Ord + 1;

actRp = PBR;
actAs = SBAtten; 

figure(77);
plot([0:1:LenFFT/2]/(LenFFT/2), 20*log10(fr+eps)); 
xlabel('Normalized Frequency (Units of \pi)')
ylabel('Magnitude, dB')
axis([0,1,-(As+20),5])