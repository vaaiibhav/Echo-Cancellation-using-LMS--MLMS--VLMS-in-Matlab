function LVLPFSincBlackman(wp,ws,Rp,As) 
% LVLPFSincBlackman(0.25*pi,0.32*pi,0.1,74)
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
wc = (wp+ws)/2; wt = ws - wp; limL = ceil(11*pi/wt), 
startL = fix(0.9*limL); figure(59); LenFFT = 8192;
for L = startL:1:ceil(1.1*limL); M = (L-1)/2;
n = 0:1:L-1; b = sin(wc*(n - M + eps))./(pi*(n - M + eps)); 
b = b.*(blackman(L)'); fr = abs(fft(b,LenFFT)); 
fr=fr(1,1:(LenFFT/2+1)); Lfr = length(fr); 
PB = fr(1,1:round((wp/pi)*Lfr)); SB = fr(1,round((ws/pi)*Lfr):Lfr); 
PBR = -20*log10(min(PB)),  SBAtten = -20*log10(max(SB)), 
plot([0:1:LenFFT/2]/(LenFFT/2), 20*log10(fr+eps)); 
xlabel('Frequency, Units of \pi');
ylabel(['Mag, dB']); axis([0 1 -100 5])
text(0.6,-10,['est L = ', num2str(limL)]);
text(0.6,-20,['L = ', num2str(L)]);
text(0.6,-30,['design Rp = ', num2str(Rp)]);
text(0.6,-40,['actual Rp = ', num2str(PBR)]);
text(0.6,-50,['design As = ', num2str(As)]);
text(0.6,-60,['actual As = ', num2str(SBAtten)]);
if SBAtten>=As; break; end; pause(0.25); end