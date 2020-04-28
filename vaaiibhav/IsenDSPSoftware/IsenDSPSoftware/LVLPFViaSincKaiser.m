function LVLPFViaSincKaiser(wc,L)
% LVLPFViaSincKaiser(0.5*pi,50)
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
M = (L-1)/2; n = 0:1:L-1;
Imp = sin(wc*(n - M + eps))./(pi*(n - M + eps)); 
win = kaiser(L,10)'; fr = abs(fft(Imp.*win, 2048));
fr = fr(1,1:fix(length(fr)/2+1));
xvec = [0:1:length(fr)-1]/length(fr);
figure(6); plot(xvec,20*log10(fr+eps))
xlabel(['Frequency, Units of \pi'])
ylabel(['Mag, dB']); axis([0 inf -110 10])