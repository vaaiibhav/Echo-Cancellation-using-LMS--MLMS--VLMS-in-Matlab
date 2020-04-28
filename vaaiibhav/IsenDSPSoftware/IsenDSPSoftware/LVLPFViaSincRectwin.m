function LVLPFViaSincRectwin(wp,ws,L)
% LVLPFViaSincRectwin(0.4*pi,0.5*pi,17)
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
    wc = (wp + ws)/2; M = (L-1)/2; n = 0:1:L-1; 
    b = sin(wc*(n - M + eps))./(pi*(n - M + eps)); 
    LenFFT = 8192; fr = abs(fft(b,LenFFT)); fr=fr(1,1:LenFFT/2+1); 
    Lfr = length(fr); PB = fr(1,1:round((wp/pi)*Lfr)); 
    SB = fr(1,round((ws/pi)*Lfr):Lfr); 
    PBR = -20*log10(min(PB)), SBAtten = -20*log10(max(SB)), 
    figure(59); plot([0:1:LenFFT/2]/(LenFFT/2), 20*log10(fr+eps)); 
    xlabel('Frequency, Units of \pi'); 
    text(0.1,-30,['actual Rp = ',num2str(PBR,3),' dB'])
    text(0.1,-45,['actual As = ',num2str(SBAtten,3),' dB'])
    ylabel(['Mag, dB']); axis([0 1 -inf inf])