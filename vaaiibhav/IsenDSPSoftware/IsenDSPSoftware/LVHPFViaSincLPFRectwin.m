function LVHPFViaSincLPFRectwin(wc,L)
% LVHPFViaSincLPFRectwin(0.3*pi,51)
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
    M = (L-1)/2; n = 0:1:L-1;
    ImpLo = sin(wc*(n - M + eps))./(pi*(n - M + eps));
    ap = sin(pi*(n - M + eps))./(pi*(n - M + eps));
    ImpHi = ap - ImpLo; frImpHi = abs(fft(ImpHi,1024));
    frImpLo = abs(fft(ImpLo,1024)); figure(56); 
    set(56,'color',[1,1,1]); subplot(221); stem(n,ImpLo); 
    xlabel('(a) Sample'); ylabel('Amplitude')
    axis([0 length(ImpLo) 1.2*min(ImpLo) 1.1*max(ImpLo)])
    subplot(222); plot([0:1:512]/512, frImpLo(1,1:513));
    xlabel(['(b) Freq, Units of \pi']); ylabel('Magnitude')
    axis([0 1 0 1.1]); subplot(223); stem(n,ImpHi); 
    xlabel('(c) Sample'); ylabel('Amplitude')
    axis([0 length(ImpHi) 1.2*min(ImpHi) 1.1*max(ImpHi)])
    subplot(224); plot([0:1:512]/512,frImpHi(1,1:513)); 
    xlabel(['(d) Freq, Units of \pi']); ylabel('Magnitude')
    axis([0 1 0 1.1])