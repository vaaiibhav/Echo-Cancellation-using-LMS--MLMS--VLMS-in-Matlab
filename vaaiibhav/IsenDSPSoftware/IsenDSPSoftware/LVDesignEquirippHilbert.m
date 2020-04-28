function LVDesignEquirippHilbert(Rp,wp1,wp2) 
%
% LVDesignEquirippHilbert(0.2,0.1,0.9)
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
PBR = 100; Ord = 10; LenFFT = 2^13;
figure(99); while (PBR > Rp); Ord = Ord + 2;
b = remez(Ord,[0.1,0.9],[1,1],'Hilbert');
y = abs(fft(b,LenFFT)); 
y = y(1,1:LenFFT/2+1)/(max(abs(y)));
LenGrid = LenFFT/2;
PB = y(1,round(wp1*LenGrid)+1:round(wp2*LenGrid)+1);
PBR = -20*log10(min(PB)+eps);
plot([0:1:LenGrid]/LenGrid,20*log10(y+eps))
xlabel(['Normalized Frequency (Units of \pi)'])
ylabel('Magnitude, dB'); L = Ord + 1,
axis([0,1,-40,5]); pause(0.3); end
Final_Rp = PBR, Final_L = Ord + 1