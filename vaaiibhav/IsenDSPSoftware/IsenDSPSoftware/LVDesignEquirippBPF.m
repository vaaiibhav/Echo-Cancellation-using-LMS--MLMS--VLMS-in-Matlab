function LVDesignEquirippBPF(Rp,As,ws1,wp1,wp2,ws2) 
% LVDesignEquirippBPF(0.2,60,0.4,0.45,0.65,0.7)
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
BndLm = [0,ws1,wp1,wp2,ws2,1];
Rfac= 10^(-Rp/20); DeltaP = (1-Rfac)/(1+Rfac);
DeltaS = (1+DeltaP)*10^(-As/20);
deltaF = abs(wp1 - ws1)/2;
compL = (-20*log10(sqrt(DeltaP*DeltaS)) - 13)/(14.6*deltaF) + 1;
L = round(compL), Ord = L - 2; LenFFT = 2^15;
SBAtten = 0; PBR = 100;
while (SBAtten < As)|(PBR > Rp); Ord = Ord + 1;
WF = firpm(Ord,BndLm,[0,0,1,1,0,0],[1,DeltaS/DeltaP,1]);
fr = abs(fft(WF,LenFFT)); 
fr = fr(1,1:LenFFT/2+1)/max(fr);
Lfr = length(fr); 
PB = fr(1,round(wp1*Lfr):round(wp2*Lfr)); 
PBR = -20*log10(min(PB)+eps), 
SB1 = fr(1,1:round(ws1*Lfr));
SB2 = fr(1,round(ws2*Lfr):Lfr); 
SBAtten1 = -20*log10(max(SB1)+eps);
SBAtten2 = -20*log10(max(SB2)+eps);
SBAtten = min([SBAtten1,SBAtten2]),
end; Fin_L = Ord + 1, Fin_Rp = PBR,
Fin_As = SBAtten, figure(8); clf;
plot([0:1:LenFFT/2]/(LenFFT/2), 20*log10(fr+eps)); 
xlabel('Normalized Frequency (Units of \pi)')
ylabel('Magnitude, dB'); axis([0,1,-(As+20),5])
