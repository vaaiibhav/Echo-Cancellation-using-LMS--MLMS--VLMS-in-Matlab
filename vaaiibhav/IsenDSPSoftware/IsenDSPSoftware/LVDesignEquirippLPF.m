function [actRp,actAs,WF] = LVDesignEquirippLPF(Rp,As,wp,ws)    
% LVDesignEquirippLPF(0.2,60,0.45,0.55)
% wp = 0.45; ws = 0.55; As = 60; Rp = 0.2;
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
Rfac= 10^(-Rp/20); DeltaP = (1-Rfac)/(1+Rfac);
DeltaS = (1+DeltaP)*10^(-As/20);
deltaF = abs(ws - wp)/2;
compL = (-20*log10(sqrt(DeltaP*DeltaS)) - 13)/(14.6*deltaF) + 1;
L = round(compL), Ord = L - 2; LenFFT = 2^15;
SBAtten = 0; PBR = 100;
while (SBAtten < As)|(PBR > Rp); Ord = Ord + 1
WF = firpm(Ord,[0,wp,ws,1],[1,1,0,0],[DeltaS/DeltaP,1]);
fr = abs(fft(WF,LenFFT));
fr = fr(1,1:LenFFT/2+1)/(max(abs(fr)));
Lfr = length(fr); PB = fr(1,1:round(wp*(Lfr-1)));
SB = fr(1,round(ws*(Lfr-1))+1:Lfr);
PBR = -20*log10(min(PB) + eps)
SBAtten = -20*log10(max(SB) + eps)
end; Fin_L = Ord + 1, actRp = PBR;
actAs = SBAtten; figure(77);
plot([0:1:LenFFT/2]/(LenFFT/2), 20*log10(fr+eps)); 
xlabel('Normalized Frequency (Units of \pi)')
ylabel('Magnitude, dB')
axis([0,1,-(As+20),5])