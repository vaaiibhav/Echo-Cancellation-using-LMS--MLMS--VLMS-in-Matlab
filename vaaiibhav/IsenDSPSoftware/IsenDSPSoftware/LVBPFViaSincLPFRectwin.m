function LVBPFViaSincLPFRectwin(wc1,wc2,L)
%
% LVBPFViaSincLPFRectwin(0.25*pi,0.55*pi,79)
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
M = (L-1)/2; n = 0:1:L-1; 
ImpLoWide = sin(wc2*(n - M + eps))./(pi*(n - M + eps));
ImpLoNarrow = sin(wc1*(n - M + eps))./(pi*(n - M + eps));
ImpBand = ImpLoWide - ImpLoNarrow; 
frImpBand = abs(fft(ImpBand,1024));
figure(57); subplot(211); stem(n,ImpBand); 
xlabel('(a) Sample'); ylabel('Amplitude')
axis([0 length(ImpBand) 1.2*min(ImpBand) 1.1*max(ImpBand)])
subplot(212); plot([0:1:512]/512,frImpBand(1,1:513));
xlabel(['(b) Frequency, Units of \pi']); 
ylabel('Magnitude'); axis([0 1 0 1.2])