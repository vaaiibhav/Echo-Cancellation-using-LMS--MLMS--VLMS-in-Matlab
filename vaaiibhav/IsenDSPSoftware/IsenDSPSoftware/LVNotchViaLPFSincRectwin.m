function LVNotchViaLPFSincRectwin(wc1,wc2,L)
% LVNotchViaLPFSincRectwin(0.45*pi,0.75*pi,71)
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
M = (L-1)/2; n = 0:1:L-1;
ImpLo2 = sin(wc2*(n - M + eps))./(pi*(n - M + eps));
ImpLo1 = sin(wc1*(n - M + eps))./(pi*(n - M + eps));
ImpBand = ImpLo2 - ImpLo1; 
ImpWide = sin(pi*(n - M + eps))./(pi*(n - M + eps));
ImpStop = ImpWide - ImpBand;
frImpStop = abs(fft(ImpStop,1024));
figure(58); subplot(211); stem(n,ImpStop); 
xlabel('(a) Sample'); ylabel('Amplitude'); 
subplot(212); plot([0:1:512]/512, frImpStop(1,1:513));
xlabel(['(b) Frequency, Units of \pi']); 
ylabel('Magnitude'); axis([0 1 0 1.2])